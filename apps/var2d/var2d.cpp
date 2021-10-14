/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured 
 *  Cartesian grids.
 *
 *  Copyright (c) 2018-2021 Scott Aiton
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

#include "Writers/ClawWriter.h"
#include <ThunderEgg.h>
#ifdef HAVE_VTK
#include "Writers/VtkWriter2d.h"
#endif
#include "CLI11.hpp"
#include "GetP4estDomainGenerator.h"
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <unistd.h>
//#include "IfaceMatrixHelper.h"

// =========== //
// main driver //
// =========== //

using namespace std;
using namespace ThunderEgg;
using namespace ThunderEgg::VarPoisson;
int main(int argc, char *argv[])
{
	PetscInitialize(nullptr, nullptr, nullptr, nullptr);
	Communicator comm(MPI_COMM_WORLD);

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

	// Set the number of discretization points in the x and y direction.
	std::array<int, 2> ns;
	ns.fill(n);

	///////////////
	// Create Mesh
	///////////////
	auto bmf = [](int block_no, double unit_x, double unit_y, double &x, double &y) {
		x = unit_x;
		y = unit_y;
	};

	p4est_connectivity * conn = p4est_connectivity_new_brick(1, 1, false, false);
	P4estDomainGenerator dcg  = getP4estDomainGenerator(conn, mesh_filename, ns, 1, bmf);

	Domain<2> domain = dcg.getFinestDomain();

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
			return 1 + x * y;
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
	}

	std::shared_ptr<Timer> timer = make_shared<Timer>(comm);
	for (int loop = 0; loop < loop_count; loop++) {
		timer->start("Domain Initialization");

		Vector<2> u(domain, 1);
		Vector<2> exact(domain, 1);
		Vector<2> f(domain, 1);
		Vector<2> au(domain, 1);
		Vector<2> h(domain, 1);

		DomainTools::SetValues<2>(domain, f, ffun);
		DomainTools::SetValues<2>(domain, exact, gfun);
		DomainTools::SetValuesWithGhost<2>(domain, h, hfun);

		timer->stop("Domain Initialization");

		// patch operator
		BiLinearGhostFiller gf(domain, GhostFillingType::Faces);
		auto                p_operator = make_shared<StarPatchOperator<2>>(h, domain, gf);
		p_operator->addDrichletBCToRHS(f, gfun, hfun);

		// set the patch solver
		Iterative::BiCGStab<2> p_bcgs;
		p_bcgs.setTolerance(ps_tol);
		p_bcgs.setMaxIterations(ps_max_it);
		auto p_solver = make_shared<Iterative::PatchSolver<2>>(p_bcgs, *p_operator);

		std::shared_ptr<Operator<2>> A = p_operator;
		std::shared_ptr<Operator<2>> M;
		///////////////////
		// setup start
		///////////////////
		timer->start("Linear System Setup");

		timer->start("Matrix Formation");
		timer->stop("Matrix Formation");
		// preconditoners
		timer->start("Preconditioner Setup");
		if (preconditioner == "GMG") {
			timer->start("GMG Setup");

			Domain<2> curr_domain = domain;

			int domain_level = 0;
			curr_domain.setTimer(timer);
			domain_level++;

			Domain<2>                next_domain = dcg.getCoarserDomain();
			GMG::LinearRestrictor<2> restrictor(curr_domain, next_domain, true);

			GMG::CycleBuilder<2> builder(copts);
			builder.addFinestLevel(*p_operator, *p_solver, restrictor);

			auto prev_coeffs = h;
			auto prev_domain = curr_domain;
			curr_domain      = next_domain;
			while (dcg.hasCoarserDomain()) {
				curr_domain.setTimer(timer);
				domain_level++;

				auto                next_domain = dcg.getCoarserDomain();
				BiLinearGhostFiller new_gf(curr_domain, GhostFillingType::Faces);
				Vector<2>           new_coeffs(curr_domain, 1);

				DomainTools::SetValuesWithGhost<2>(curr_domain, new_coeffs, hfun);

				StarPatchOperator<2> new_p_operator(new_coeffs, curr_domain, new_gf);

				Iterative::PatchSolver<2> new_p_solver(p_bcgs, new_p_operator);

				GMG::DirectInterpolator<2> interpolator(curr_domain, prev_domain);
				restrictor = GMG::LinearRestrictor<2>(curr_domain, next_domain);

				builder.addIntermediateLevel(new_p_operator, new_p_solver, restrictor, interpolator);
				prev_domain = curr_domain;
				curr_domain = next_domain;
			}
			curr_domain.setTimer(timer);

			GMG::DirectInterpolator<2> interpolator(curr_domain, prev_domain);
			BiLinearGhostFiller        coarse_gf(curr_domain, GhostFillingType::Faces);
			Vector<2>                  coarse_coeffs(curr_domain, 1);
			DomainTools::SetValuesWithGhost<2>(curr_domain, coarse_coeffs, hfun);

			StarPatchOperator<2> coarse_p_operator(coarse_coeffs, curr_domain, coarse_gf);

			Iterative::PatchSolver<2> coarse_p_solver(p_bcgs, coarse_p_operator);
			builder.addCoarsestLevel(coarse_p_operator, coarse_p_solver, interpolator);

			M = builder.getCycle();

			timer->stop("GMG Setup");
		}
		timer->stop("Preconditioner Setup");

		timer->stop("Linear System Setup");

		timer->start("Linear Solve");
		u.set(0);
		Iterative::BiCGStab<2> solver;
		solver.setMaxIterations(1000);
		solver.setTolerance(1e-12);
		solver.setTimer(timer);
		int its = solver.solve(*A, u, f, M.get());
		if (comm.getRank() == 0) {
			cout << "Iterations: " << its << endl;
		}
		timer->stop("Linear Solve");

		A->apply(u, au);

		// residual
		Vector<2> resid(domain, 1);
		resid.scaleThenAddScaled(0, -1, au, 1, f);
		double residual = resid.twoNorm();
		double fnorm    = f.twoNorm();

		// error
		Vector<2> error(domain, 1);
		error.scaleThenAddScaled(0, -1, exact, 1, u);
		if (neumann) {
			double uavg = DomainTools::Integrate<2>(domain, u) / domain.volume();
			double eavg = DomainTools::Integrate<2>(domain, exact) / domain.volume();

			if (comm.getRank() == 0) {
				cout << "Average of computed solution: " << uavg << endl;
				cout << "Average of exact solution: " << eavg << endl;
			}

			error.shift(eavg - uavg);
		}
		double error_norm = error.twoNorm();
		double exact_norm = exact.twoNorm();

		double ausum = DomainTools::Integrate<2>(domain, au);
		double fsum  = DomainTools::Integrate<2>(domain, f);
		if (comm.getRank() == 0) {
			std::cout << std::scientific;
			std::cout.precision(13);
			std::cout << "Error: " << error_norm / exact_norm << endl;
			std::cout << "Error-inf: " << error.infNorm() << endl;
			std::cout << "Residual: " << residual / fnorm << endl;
			std::cout << u8"ΣAu-Σf: " << ausum - fsum << endl;
			cout.unsetf(std::ios_base::floatfield);
			int total_cells = domain.getNumGlobalCells();
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
	cout << *timer;
	timer->saveToFile("timings.json");
	PetscFinalize();
	return 0;
}
