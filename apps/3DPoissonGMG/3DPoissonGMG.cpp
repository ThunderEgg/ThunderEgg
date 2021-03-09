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

#include "CLI11.hpp"
#include <ThunderEgg.h>
#ifdef HAVE_VTK
#include "Writers/VtkWriter.h"
#endif
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <unistd.h>

// =========== //
// main driver //
// =========== //

using namespace std;
using namespace ThunderEgg;
using namespace ThunderEgg::Experimental;
using namespace ThunderEgg::Poisson;

int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);

	// parse input
	CLI::App app{"ThunderEgg 3d poisson solver example"};

	app.set_config("--config", "", "Read an ini file", false);
	// program options
	int n;
	app.add_option("-n,--num_cells", n, "Number of cells in each direction, on each patch")
	->required();

	int div = 0;
	app.add_option("--divide", div, "Number of levels to add to oct-tree");

	bool no_zero_rhs_avg = false;
	app.add_flag("--nozerof", no_zero_rhs_avg,
	             "Make the average of the rhs on neumann problems zero");

	double tolerance = 1e-12;
	app.add_option("-t,--tolerance", tolerance, "Tolerance of Krylov solver");

	bool neumann;
	app.add_flag("--neumann", neumann, "Use neumann boundary conditions");

	string mesh_filename = "";
	app.add_option("--mesh", mesh_filename, "Filename of mesh to use")
	->required()
	->check(CLI::ExistingFile);

	string problem = "trig";
	app.add_set_ignore_case("--problem", problem, {"trig", "gauss", "zero"},
	                        "Which problem to solve");

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

	string patch_solver = "fftw";
	gmg->add_set_ignore_case("--patch_solver", patch_solver, {"fftw", "dft", "bcgs"},
	                         "Which patch solver to use");
	// output options

#ifdef HAVE_VTK
	string vtk_filename = "";
	app.add_option("--out_vtk", vtk_filename, "Filename of vtk output");
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

	int num_procs;
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	int my_global_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_global_rank);

	// Set the number of discretization points in the x and y direction.
	std::array<int, 3> ns;
	ns.fill(n);

	///////////////
	// Create Mesh
	///////////////
	shared_ptr<Domain<3>> domain;
	Tree<3>               t;
	t = Tree<3>(mesh_filename);
	for (int i = 0; i < div; i++) {
		t.refineLeaves();
	}

	// the functions that we are using
	function<double(const std::array<double, 3> &)> ffun;
	function<double(const std::array<double, 3> &)> gfun;
	function<double(const std::array<double, 3> &)> nfunx;
	function<double(const std::array<double, 3> &)> nfuny;
	function<double(const std::array<double, 3> &)> nfunz;

	if (problem == "zero") {
		ffun  = [](const std::array<double, 3> &) { return 0; };
		gfun  = [](const std::array<double, 3> &) { return 0; };
		nfunx = [](const std::array<double, 3> &) { return 0; };
		nfuny = [](const std::array<double, 3> &) { return 0; };
		nfunz = [](const std::array<double, 3> &) { return 0; };
	} else if (problem == "gauss") {
		gfun = [](const std::array<double, 3> &coord) {
			double x = coord[0];
			double y = coord[1];
			double z = coord[2];
			return exp(cos(10 * M_PI * x)) - exp(cos(11 * M_PI * y)) + exp(cos(12 * M_PI * z));
		};
		ffun = [](const std::array<double, 3> &coord) {
			double x = coord[0];
			double y = coord[1];
			double z = coord[2];
			return -M_PI * M_PI
			       * (100 * exp(cos(10 * M_PI * x)) * cos(10 * M_PI * x)
			          - 100 * exp(cos(10 * M_PI * x)) * pow(sin(10 * M_PI * x), 2)
			          - 121 * exp(cos(11 * M_PI * y)) * cos(11 * M_PI * y)
			          + 121 * exp(cos(11 * M_PI * y)) * pow(sin(11 * M_PI * y), 2)
			          + 144 * exp(cos(12 * M_PI * z)) * cos(12 * M_PI * z)
			          - 144 * exp(cos(12 * M_PI * z)) * pow(sin(12 * M_PI * z), 2));
		};
		nfunx = [](const std::array<double, 3> &coord) {
			double x = coord[0];
			return -10 * M_PI * sin(10 * M_PI * x) * exp(cos(10 * M_PI * x));
		};
		nfuny = [](const std::array<double, 3> &coord) {
			double y = coord[1];
			return 11 * M_PI * sin(11 * M_PI * y) * exp(cos(11 * M_PI * y));
		};
		nfunz = [](const std::array<double, 3> &coord) {
			double z = coord[2];
			return -12 * M_PI * sin(12 * M_PI * z) * exp(cos(12 * M_PI * z));
		};
	} else {
		ffun = [](const std::array<double, 3> &coord) {
			double x = coord[0];
			double y = coord[1];
			double z = coord[2];
			x += .3;
			y += .3;
			z += .3;
			return -77.0 / 36 * M_PI * M_PI * sin(M_PI * x) * cos(2.0 / 3 * M_PI * y)
			       * sin(5.0 / 6 * M_PI * z);
		};
		gfun = [](const std::array<double, 3> &coord) {
			double x = coord[0];
			double y = coord[1];
			double z = coord[2];
			x += .3;
			y += .3;
			z += .3;
			return sin(M_PI * x) * cos(2.0 / 3 * M_PI * y) * sin(5.0 / 6 * M_PI * z);
		};
		nfunx = [](const std::array<double, 3> &coord) {
			double x = coord[0];
			double y = coord[1];
			double z = coord[2];
			x += .3;
			y += .3;
			z += .3;
			return M_PI * cos(M_PI * x) * cos(2.0 / 3 * M_PI * y) * sin(5.0 / 6 * M_PI * z);
		};
		nfuny = [](const std::array<double, 3> &coord) {
			double x = coord[0];
			double y = coord[1];
			double z = coord[2];
			x += .3;
			y += .3;
			z += .3;
			return -2.0 / 3 * M_PI * sin(M_PI * x) * sin(2.0 / 3 * M_PI * y)
			       * sin(5.0 / 6 * M_PI * z);
		};
		nfunz = [](const std::array<double, 3> &coord) {
			double x = coord[0];
			double y = coord[1];
			double z = coord[2];
			x += .3;
			y += .3;
			z += .3;
			return 5.0 / 6 * M_PI * sin(M_PI * x) * cos(2.0 / 3 * M_PI * y)
			       * cos(5.0 / 6 * M_PI * z);
		};
	}

	shared_ptr<DomainGenerator<3>> dcg(new DomGen<3>(t, ns, 1, neumann));

	domain = dcg->getFinestDomain();

	shared_ptr<GhostFiller<3>> ghost_filler = make_shared<TriLinearGhostFiller>(domain);
	// patch operator
	shared_ptr<StarPatchOperator<3>> patch_operator
	= make_shared<StarPatchOperator<3>>(domain, ghost_filler, neumann);

	std::shared_ptr<Timer> timer = make_shared<Timer>(MPI_COMM_WORLD);
	timer->start("Domain Initialization");

	// Initialize Vectors
	shared_ptr<ValVector<3>> u     = ValVector<3>::GetNewVector(domain, 1);
	shared_ptr<ValVector<3>> exact = ValVector<3>::GetNewVector(domain, 1);
	shared_ptr<ValVector<3>> f     = ValVector<3>::GetNewVector(domain, 1);
	shared_ptr<ValVector<3>> au    = ValVector<3>::GetNewVector(domain, 1);

	DomainTools::SetValues<3>(domain, f, ffun);
	DomainTools::SetValues<3>(domain, exact, gfun);
	if (neumann) {
		patch_operator->addDrichletBCToRHS(f, gfun);
	} else {
		patch_operator->addNeumannBCToRHS(f, gfun, {nfunx, nfuny, nfunz});
	}

	timer->stop("Domain Initialization");

	if (neumann && !no_zero_rhs_avg) {
		double fdiff = domain->integrate(f) / domain->volume();
		if (my_global_rank == 0)
			cout << "Fdiff: " << fdiff << endl;
		f->shift(-fdiff);
	}

	std::shared_ptr<Operator<3>> A;
	std::shared_ptr<Operator<3>> M;
	///////////////////
	// setup start
	///////////////////
	timer->start("Linear System Setup");

	A = patch_operator;
	// preconditoners
	timer->start("Preconditioner Setup");

	auto curr_domain = domain;

	int domain_level = 0;
	curr_domain->setId(domain_level);
	curr_domain->setTimer(timer);
	domain_level++;

	auto next_domain = dcg->getCoarserDomain();
	auto restrictor  = make_shared<GMG::LinearRestrictor<3>>(curr_domain, next_domain, 1, true);

	// set the patch solver
	auto p_bcgs = make_shared<Iterative::BiCGStab<3>>();
	// p_bcgs->setTolerance(ps_tol);
	// p_bcgs->setMaxIterations(ps_max_it);

	// set the patch solver
	shared_ptr<PatchSolver<3>> p_solver;
	if (patch_solver == "dft") {
		p_solver = make_shared<DFTPatchSolver<3>>(patch_operator);
	} else if (patch_solver == "fftw") {
		p_solver = make_shared<FFTWPatchSolver<3>>(patch_operator);
	} else {
		p_solver = make_shared<Iterative::PatchSolver<3>>(p_bcgs, patch_operator);
	}

	GMG::CycleBuilder<3>                builder(copts);
	std::shared_ptr<VectorGenerator<3>> vg(new ValVectorGenerator<3>(domain, 1));
	builder.addFinestLevel(patch_operator, p_solver, restrictor, vg);

	auto prev_domain = curr_domain;
	curr_domain      = next_domain;
	while (dcg->hasCoarserDomain()) {
		curr_domain->setId(domain_level);
		curr_domain->setTimer(timer);
		domain_level++;

		auto next_domain = dcg->getCoarserDomain();
		auto new_vg      = make_shared<ValVectorGenerator<3>>(curr_domain, 1);
		auto new_gf      = make_shared<TriLinearGhostFiller>(curr_domain);
		auto new_coeffs  = new_vg->getNewVector();

		auto new_p_operator = make_shared<StarPatchOperator<3>>(curr_domain, new_gf);

		shared_ptr<PatchSolver<3>> new_p_solver;
		if (patch_solver == "dft") {
			new_p_solver = make_shared<DFTPatchSolver<3>>(new_p_operator);
		} else if (patch_solver == "fftw") {
			new_p_solver = make_shared<FFTWPatchSolver<3>>(new_p_operator);
		} else {
			new_p_solver = make_shared<Iterative::PatchSolver<3>>(p_bcgs, new_p_operator);
		}

		auto interpolator = make_shared<GMG::DirectInterpolator<3>>(curr_domain, prev_domain, 1);
		restrictor        = make_shared<GMG::LinearRestrictor<3>>(curr_domain, next_domain, 1);

		builder.addIntermediateLevel(new_p_operator, new_p_solver, restrictor, interpolator,
		                             new_vg);
		prev_domain = curr_domain;
		curr_domain = next_domain;
	}
	curr_domain->setId(domain_level);
	curr_domain->setTimer(timer);

	auto interpolator = make_shared<GMG::DirectInterpolator<3>>(curr_domain, prev_domain, 1);
	auto coarse_vg    = make_shared<ValVectorGenerator<3>>(curr_domain, 1);
	auto coarse_gf    = make_shared<TriLinearGhostFiller>(curr_domain);

	auto coarse_p_operator = make_shared<StarPatchOperator<3>>(curr_domain, coarse_gf);

	shared_ptr<PatchSolver<3>> coarse_p_solver;
	if (patch_solver == "dft") {
		coarse_p_solver = make_shared<DFTPatchSolver<3>>(coarse_p_operator);
	} else if (patch_solver == "fftw") {
		coarse_p_solver = make_shared<FFTWPatchSolver<3>>(coarse_p_operator);
	} else {
		coarse_p_solver = make_shared<Iterative::PatchSolver<3>>(p_bcgs, coarse_p_operator);
	}
	builder.addCoarsestLevel(coarse_p_operator, coarse_p_solver, interpolator, coarse_vg);

	M = builder.getCycle();

	timer->stop("Preconditioner Setup");

	timer->stop("Linear System Setup");

	timer->start("Linear Solve");

	Iterative::BiCGStab<3> solver;
	solver.setTimer(timer);
	int its = solver.solve(vg, A, u, f, M, true);
	if (my_global_rank == 0) {
		cout << "Iterations: " << its << endl;
	}
	timer->stop("Linear Solve");

	A->apply(u, au);

	shared_ptr<ValVector<3>> resid = ValVector<3>::GetNewVector(domain, 1);
	resid->addScaled(-1, au, 1, f);

	double residual = resid->twoNorm();
	double fnorm    = f->twoNorm();

	// error
	shared_ptr<ValVector<3>> error = ValVector<3>::GetNewVector(domain, 1);
	error->addScaled(-1, exact, 1, u);
	if (neumann) {
		double uavg = domain->integrate(u) / domain->volume();
		double eavg = domain->integrate(exact) / domain->volume();

		if (my_global_rank == 0) {
			cout << "Average of computed solution: " << uavg << endl;
			cout << "Average of exact solution: " << eavg << endl;
		}

		error->shift(eavg - uavg);
	}
	double error_norm     = error->twoNorm();
	double error_norm_inf = error->infNorm();
	double exact_norm     = exact->twoNorm();

	double ausum = domain->integrate(au);
	double fsum  = domain->integrate(f);
	if (my_global_rank == 0) {
		std::cout << std::scientific;
		std::cout.precision(13);
		std::cout << "Error (2-norm):   " << error_norm / exact_norm << endl;
		std::cout << "Error (inf-norm): " << error_norm_inf << endl;
		std::cout << "Residual: " << residual / fnorm << endl;
		std::cout << u8"ΣAu-Σf: " << ausum - fsum << endl;
		cout.unsetf(std::ios_base::floatfield);
		int total_cells = domain->getNumGlobalCells();
		cout << "Total cells: " << total_cells << endl;
		cout << "Cores: " << num_procs << endl;
	}

	// output
#ifdef HAVE_VTK
	if (vtk_filename != "") {
		VtkWriter writer(dc, vtk_filename);
		writer.add(u->vec, "Solution");
		writer.add(error->vec, "Error");
		writer.add(exact->vec, "Exact");
		writer.add(resid->vec, "Residual");
		writer.add(f->vec, "RHS");
		writer.write();
	}
#endif
	cout.unsetf(std::ios_base::floatfield);

	if (my_global_rank == 0) {
		cout << *timer;
	}
	MPI_Finalize();
	return 0;
}
