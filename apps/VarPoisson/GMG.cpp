/***************************************************************************
 *  ThunderEgg, a library for solving Poisson's equation on adaptively
 *  refined block-structured Cartesian grids
 *
 *  Copyright (C) 2019  ThunderEgg Developers. See AUTHORS.md file at the
 *  top-level directory.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 ***************************************************************************/

#include "Init.h"
#include "Writers/ClawWriter.h"
#include <ThunderEgg/BiCGStab.h>
#include <ThunderEgg/BiCGStabPatchSolver.h>
#include <ThunderEgg/BiLinearGhostFiller.h>
#include <ThunderEgg/Domain.h>
#include <ThunderEgg/DomainTools.h>
#include <ThunderEgg/DomainWrapOp.h>
#include <ThunderEgg/Experimental/DomGen.h>
#include <ThunderEgg/GMG/CycleFactory.h>
#include <ThunderEgg/GMG/DrctIntp.h>
#include <ThunderEgg/GMG/LinearRestrictor.h>
#include <ThunderEgg/MatWrapper.h>
#include <ThunderEgg/PetscShellCreator.h>
#include <ThunderEgg/SchwarzPrec.h>
#include <ThunderEgg/Timer.h>
#include <ThunderEgg/VarPoisson/StarPatchOperator.h>
#ifdef HAVE_VTK
#include "Writers/VtkWriter2d.h"
#endif
#ifdef HAVE_P4EST
#include "TreeToP4est.h"
#include <ThunderEgg/P4estDomainGenerator.h>
#endif
#include "CLI11.hpp"
#include <cmath>
#include <iostream>
#include <memory>
#include <mpi.h>

// =========== //
// main driver //
// =========== //

using namespace std;
using namespace ThunderEgg;
using namespace ThunderEgg::Experimental;
using namespace ThunderEgg::VarPoisson;
int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);
	PetscInitialize(nullptr, nullptr, nullptr, nullptr);

	// parse input
	CLI::App app{"ThunderEgg 3d poisson solver example"};

	app.set_config("--config", "", "Read an ini file", false);
	// program options
	int n;
	app.add_option("-n,--num_cells", n, "Number of cells in each direction, on each patch")
	->required();

	int loop_count = 1;
	app.add_option("-l", loop_count, "Number of times to run program");

	string matrix_type = "wrap";
	app.add_set_ignore_case("--matrix_type", matrix_type, {"wrap", "crs"},
	                        "Which type of matrix operator to use");

	int div = 0;
	app.add_option("--divide", div, "Number of levels to add to quadtree");

	bool no_zero_rhs_avg = false;
	app.add_flag("--nozerof", no_zero_rhs_avg,
	             "Make the average of the rhs on neumann problems zero");

	double tolerance = 1e-12;
	app.add_option("-t,--tolerance", tolerance, "Tolerance of Krylov solver");

	bool neumann;
	app.add_flag("--neumann", neumann, "Use neumann boundary conditions");

	bool solve_schur = false;
	app.add_flag("--schur", solve_schur, "Solve the schur compliment system");

	string mesh_filename = "";
	app.add_option("--mesh", mesh_filename, "Filename of mesh to use")
	->required()
	->check(CLI::ExistingFile);

	string problem = "trig";
	app.add_set_ignore_case("--problem", problem, {"trig", "gauss", "zero", "circle"},
	                        "Which problem to solve");

	string solver_type = "thunderegg";
	app.add_set_ignore_case("--solver", solver_type, {"petsc", "thunderegg"},
	                        "Which Solver to use");

	string preconditioner = "";
	app.add_set_ignore_case("--prec", preconditioner, {"GMG", "Schwarz"},
	                        "Which Preconditoner to use");

	string patch_solver = "fftw";
	app.add_set_ignore_case("--patch_solver", patch_solver, {"fftw", "bcgs", "dft"},
	                        "Which patch solver to use");

	int ps_max_it = 1000;
	app.add_option("--patch_solver_max_it", ps_max_it, "max iterations for bcgs patch solver");
	double ps_tol = 1e-12;
	app.add_option("--patch_solver_tol", ps_tol, "tolerance for bcgs patch solver");

	bool setrow = false;
	app.add_flag("--setrow", setrow, "Set row in matrix");

	string petsc_opts = "";
	app.add_option("--petsc_opts", petsc_opts, "petsc options");

	// GMG options

	auto gmg = app.add_subcommand("GMG", "GMG solver options");

	GMG::CycleOpts copts;

	gmg->add_option("--max_levels", copts.max_levels,
	                "The max number of levels in GMG cycle. 0 means no limit.");

	gmg->add_option(
	"--patches_per_proc", copts.patches_per_proc,
	"Lowest level is guaranteed to have at least this number of patches per processor.");

	gmg->add_option("--pre_sweeps", copts.pre_sweeps, "Number of sweeps on down cycle");

	gmg->add_option("--post_sweeps", copts.post_sweeps, "Number of sweeps on up cycle");

	gmg->add_option("--mid_sweeps", copts.mid_sweeps,
	                "Number of sweeps inbetween up and down cycle");

	gmg->add_option("--coarse_sweeps", copts.coarse_sweeps, "Number of sweeps on coarse level");

	gmg->add_option("--cycle_type", copts.cycle_type, "Cycle type");

	// output options

	string claw_filename = "";
	app.add_option("--out_claw", claw_filename, "Filename of clawpack output");

#ifdef HAVE_VTK
	string vtk_filename = "";
	app.add_option("--out_vtk", vtk_filename, "Filename of vtk output");
#endif

	string matrix_filename = "";
	app.add_option("--out_matrix", matrix_filename, "Filename of matrix output");
	string solution_filename = "";
	app.add_option("--out_solution", solution_filename, "Filename of solution output");
	string resid_filename = "";
	app.add_option("--out_resid", resid_filename, "Filename of residual output");
	string error_filename = "";
	app.add_option("--out_error", error_filename, "Filename of error output");
	string rhs_filename = "";
	app.add_option("--out_rhs", rhs_filename, "Filename of rhs output");
	string gamma_filename = "";
	app.add_option("--out_gamma", gamma_filename, "Filename of gamma output");

#ifdef HAVE_P4EST
	bool use_p4est = false;
	app.add_flag("--p4est", use_p4est, "use p4est");
#endif

	string config_out_filename = "";
	auto   out_config_opt
	= app.add_option("--output_config", config_out_filename, "Save CLI options to config file");

	CLI11_PARSE(app, argc, argv);

	if (config_out_filename != "") {
		app.remove_option(out_config_opt);
		ofstream file_out(config_out_filename);
		file_out << app.config_to_str(true, true);
		file_out.close();
	}

	PetscOptionsInsertString(nullptr, petsc_opts.c_str());

	int num_procs;
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	int my_global_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_global_rank);

	// Set the number of discretization points in the x and y direction.
	std::array<int, 2> ns;
	ns.fill(n);

	///////////////
	// Create Mesh
	///////////////
	shared_ptr<Domain<2>> domain;
	Tree<2>               t;
	t = Tree<2>(mesh_filename);
	for (int i = 0; i < div; i++) {
		t.refineLeaves();
	}

	shared_ptr<DomainGenerator<2>> dcg;
#ifdef HAVE_P4EST
	if (use_p4est) {
		TreeToP4est ttp(t);

		auto bmf = [](int block_no, double unit_x, double unit_y, double &x, double &y) {
			x = unit_x;
			y = unit_y;
		};
		auto inf
		= [=](Side<2> s, const array<double, 2> &, const array<double, 2> &) { return neumann; };

		dcg.reset(new P4estDomainGenerator(ttp.p4est, ns, 1, inf, bmf));
#else
	if (false) {
#endif
	} else {
		dcg.reset(new DomGen<2>(t, ns, 1, neumann));
	}

	domain = dcg->getFinestDomain();

	// the functions that we are using
	function<double(const std::array<double, 2> &)> ffun;
	function<double(const std::array<double, 2> &)> gfun;
	function<double(double, double)>                nfunx;
	function<double(double, double)>                nfuny;
	function<double(const std::array<double, 2> &)> hfun;

	if (false) {
		ffun = [](const std::array<double, 2> &coord) {
			double x = coord[0];
			double y = coord[1];
			return exp(x * y)
			       * (pow(x, 5) * (-1 + y) * pow(y, 2) + y * (-2 + (5 - 3 * y) * y)
			          - pow(x, 4) * y * (4 + (-7 + y) * y)
			          - 2 * x * (1 + 2 * (-2 + y) * (-1 + y) * pow(y, 2))
			          - pow(x, 2) * (-5 + 8 * y + (-6 + y) * (-1 + y) * pow(y, 3))
			          + pow(x, 3) * (-3 + y * (12 - 6 * y - pow(y, 3) + pow(y, 4))));
		};
		gfun = [](const std::array<double, 2> &coord) {
			double x = coord[0];
			double y = coord[1];
			return (1 - x) * x * (1 - y) * y * exp(x * y);
		};
		hfun = [](const std::array<double, 2> &coord) {
			double x = coord[0];
			double y = coord[1];
			// return 1 + x * y;
			return 1;
		};
	} else {
		ffun = [](const std::array<double, 2> &coord) {
			double x = coord[0];
			double y = coord[1];
			return -5 * M_PI * M_PI * sinl(M_PI * y) * cosl(2 * M_PI * x);
		};
		gfun = [](const std::array<double, 2> &coord) {
			double x = coord[0];
			double y = coord[1];
			return sinl(M_PI * y) * cosl(2 * M_PI * x);
		};
		hfun = [](const std::array<double, 2> &coord) { return 1; };
		/*
		ffun = [](const std::array<double, 2> &coord) {
		    return 0;
		};
		gfun = [](const std::array<double, 2> &coord) {
		    return 1;
		};
		hfun = [](const std::array<double, 2> &coord) { return 1; };
		*/
	}

	std::shared_ptr<Timer> timer = make_shared<Timer>(MPI_COMM_WORLD);
	for (int loop = 0; loop < loop_count; loop++) {
		timer->start("Domain Initialization");

		/*
		if (f_outim) {
		    IfaceMatrixHelper imh(dc);
		    PW<Mat>           IA = imh.formCRSMatrix();
		    PetscViewer       viewer;
		    PetscViewerBinaryOpen(PETSC_COMM_WORLD, args::get(f_outim).c_str(), FILE_MODE_WRITE,
		                          &viewer);
		    MatView(IA, viewer);
		    PetscViewerDestroy(&viewer);
		}
		*/

		shared_ptr<ValVector<2>>   u     = ValVector<2>::GetNewVector(domain);
		shared_ptr<ValVector<2>>   exact = ValVector<2>::GetNewVector(domain);
		shared_ptr<ValVector<2>>   f     = ValVector<2>::GetNewVector(domain);
		shared_ptr<ValVector<2>>   au    = ValVector<2>::GetNewVector(domain);
		shared_ptr<ValVector<2>>   h     = ValVector<2>::GetNewVector(domain);
		shared_ptr<PetscVector<1>> h_bc  = PetscVector<2>::GetNewBCVector(domain);

		DomainTools<2>::setValues(domain, f, ffun);
		DomainTools<2>::setValues(domain, exact, gfun);
		DomainTools<2>::setValuesWithGhost(domain, h, hfun);
		DomainTools<2>::setBCValues(domain, h_bc, hfun);
		if (neumann) {
			// Init::initNeumann2d(*domain, f->vec, exact->vec, ffun, gfun, nfunx, nfuny);
		} else {
			// DomainTools<2>::addDirichlet(domain, f, gfun);
		}

		timer->stop("Domain Initialization");

		// patch operator
		shared_ptr<BiLinearGhostFiller>  gf(new BiLinearGhostFiller(domain));
		shared_ptr<StarPatchOperator<2>> p_operator(new StarPatchOperator<2>(h, domain, gf));
		// shared_ptr<fivePoint> p_operator(new fivePoint(domain, gf));
		StarPatchOperator<2>::addDrichletBCToRHS(domain, f, gfun, hfun);

		// set the patch solver
		auto p_solver
		= make_shared<BiCGStabPatchSolver<2>>(domain, gf, p_operator, ps_tol, ps_max_it);
		p_solver->setTimer(timer);

		if (neumann && !no_zero_rhs_avg) {
			double fdiff = domain->integrate(f) / domain->volume();
			if (my_global_rank == 0)
				cout << "Fdiff: " << fdiff << endl;
			f->shift(-fdiff);
		}

		std::shared_ptr<Operator<2>> A = p_operator;
		PW<Mat>                      A_petsc;
		std::shared_ptr<Operator<2>> M;
		///////////////////
		// setup start
		///////////////////
		timer->start("Linear System Setup");

		timer->start("Matrix Formation");
		timer->stop("Matrix Formation");
		// preconditoners
		timer->start("Preconditioner Setup");
		if (preconditioner == "Schwarz") {
			M = p_solver;
			//	M.reset(new SchwarzPrec<2>(sch, p_solver, p_interp));
		} else if (preconditioner == "GMG") {
			timer->start("GMG Setup");

			BiLinearGhostFiller::Generator  filler_gen(gf);
			StarPatchOperator<2>::Generator op_gen(p_operator, filler_gen);
			// fivePoint::Generator     op_gen(p_operator, filler_gen);
			BiCGStabPatchSolver<2>::Generator   smooth_gen(p_solver, filler_gen, op_gen);
			GMG::LinearRestrictor<2>::Generator restrictor_gen;
			GMG::DrctIntp<2>::Generator         interpolator_gen;
			M = GMG::CycleFactory<2>::getCycle(copts, dcg, restrictor_gen, interpolator_gen,
			                                   smooth_gen, op_gen);

			timer->stop("GMG Setup");
		} else if (preconditioner == "cheb") {
			throw 3;
		}
		timer->stop("Preconditioner Setup");

		// PW<KSP> solver;
		// setup petsc if needed
		if (solver_type == "petsc") {
			timer->start("Petsc Setup");

			// KSPCreate(MPI_COMM_WORLD, &solver);
			// KSPSetFromOptions(solver);
			// KSPSetOperators(solver, A_petsc, A_petsc);
			if (M != nullptr) {
				PC M_petsc;
				// KSPGetPC(solver, &M_petsc);
				PetscShellCreator::getPCShell(M_petsc, M, domain);
			}
			// KSPSetUp(solver);
			// KSPSetTolerances(solver, tolerance, PETSC_DEFAULT, PETSC_DEFAULT, 5000);

			timer->stop("Petsc Setup");
		}
		///////////////////
		// setup end
		///////////////////
		timer->stop("Linear System Setup");

		timer->start("Linear Solve");
		if (solver_type == "petsc") {
			// KSPSolve(solver, f->vec, u->vec);
			int its;
			// KSPGetIterationNumber(solver, &its);
			if (my_global_rank == 0) {
				cout << "Iterations: " << its << endl;
			}
		} else {
			std::shared_ptr<VectorGenerator<2>> vg(new ValVectorGenerator<2>(domain));

			u->set(0);
			int its = BiCGStab<2>::solve(vg, A, u, f, M, 1000, 1e-12, timer);
			if (my_global_rank == 0) {
				cout << "Iterations: " << its << endl;
			}
		}
		timer->stop("Linear Solve");

		A->apply(u, au);

		// residual
		shared_ptr<ValVector<2>> resid = ValVector<2>::GetNewVector(domain);
		// VecAXPBYPCZ(resid->vec, -1.0, 1.0, 0.0, au->vec, f->vec);
		double residual = resid->twoNorm();
		double fnorm    = f->twoNorm();

		// error
		shared_ptr<ValVector<2>> error = ValVector<2>::GetNewVector(domain);
		// VecAXPBYPCZ(error->vec, -1.0, 1.0, 0.0, exact->vec, u->vec);
		if (neumann) {
			double uavg = domain->integrate(u) / domain->volume();
			double eavg = domain->integrate(exact) / domain->volume();

			if (my_global_rank == 0) {
				cout << "Average of computed solution: " << uavg << endl;
				cout << "Average of exact solution: " << eavg << endl;
			}

			error->shift(eavg - uavg);
		}
		double error_norm = error->twoNorm();
		double exact_norm = exact->twoNorm();

		double ausum = domain->integrate(au);
		double fsum  = domain->integrate(f);
		if (my_global_rank == 0) {
			std::cout << std::scientific;
			std::cout.precision(13);
			std::cout << "Error: " << error_norm / exact_norm << endl;
			std::cout << "Error-inf: " << error->infNorm() << endl;
			std::cout << "Residual: " << residual / fnorm << endl;
			std::cout << u8"ΣAu-Σf: " << ausum - fsum << endl;
			cout.unsetf(std::ios_base::floatfield);
			int total_cells = domain->getNumGlobalCells();
			cout << "Total cells: " << total_cells << endl;
		}

		// output
		if (claw_filename != "") {
			ClawWriter writer(domain);
			writer.addVector(u);
			writer.addVector(f);
			writer.addVector(error);
			writer.write();
		}
#ifdef HAVE_VTK
		if (vtk_filename != "") {
			VtkWriter2d writer(domain, vtk_filename);
			writer.add(u, "Solution");
			writer.add(error, "Error");
			writer.add(resid, "Residual");
			writer.add(h, "Coeff");
			writer.add(f, "RHS");
			writer.add(exact, "Exact");
			writer.write();
		}
#endif
	}
	cout.unsetf(std::ios_base::floatfield);
	if (my_global_rank == 0) {
		cout << *timer;
	}
	PetscFinalize();
	return 0;
}
