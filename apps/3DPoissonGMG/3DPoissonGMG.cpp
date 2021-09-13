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
#include "TreeToP8est.h"
#include <cmath>

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

	// command line options
	CLI::App app{"ThunderEgg 3d poisson solver example"};

	app.set_config("--config", "", "Read an ini file", false);

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

	// the functions that we are using
	function<double(const std::array<double, 3> &)> ffun; // rhs
	function<double(const std::array<double, 3> &)> gfun; // lhs
	// functions needed for neumann boundary conditions
	function<double(const std::array<double, 3> &)> nfunx; // lhs derivative in x direction
	function<double(const std::array<double, 3> &)> nfuny; // lhs derivative in y direction
	function<double(const std::array<double, 3> &)> nfunz; // lhs derivative in z direction

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

	// create a Timer object to keep track of timings
	std::shared_ptr<Timer> timer = make_shared<Timer>(MPI_COMM_WORLD);

	// read in a tree from file
	Tree<3> t;
	t = Tree<3>(mesh_filename);
	for (int i = 0; i < div; i++) {
		t.refineLeaves();
	}
	TreeToP8est ttp(t);

	auto bmf = [](int block_no, double unit_x, double unit_y, double unit_z, double &x, double &y, double &z) {
		x = unit_x;
		y = unit_y;
		z = unit_z;
	};

	// Set the number of cells in the x, y and z direction.
	std::array<int, 3> ns;
	ns[0] = n;
	ns[1] = n;
	ns[2] = n;

	// A DomainGenerator will create domains for the Multigrid algorithm from a tree
	int num_ghost_cells = 1; // the poission operator needs 1 row/column of ghost cells on the edges of a patch
	//
	shared_ptr<DomainGenerator<3>> domain_generator = make_shared<P8estDomainGenerator>(ttp.p8est, ns, 1, bmf);

	// Get the finest domain from the tree
	shared_ptr<Domain<3>> domain = domain_generator->getFinestDomain();

	// A patch operator needs a GhostFiller object to define how to fill ghost cells for the patches.
	// This one will use a tri-linear interpolation scheme at the refinement boundarys of the domain
	shared_ptr<GhostFiller<3>> ghost_filler = make_shared<TriLinearGhostFiller>(domain, GhostFillingType::Faces);

	// create patch operator that uses a typical 2nd order 7 point poisson stencil
	shared_ptr<StarPatchOperator<3>> patch_operator = make_shared<StarPatchOperator<3>>(domain, ghost_filler, neumann);

	timer->start("Domain Initialization");

	// Create some new vectors for the domain
	int num_components = 1; //the poisson operator just has one value in each cell

	shared_ptr<ValVector<3>> u     = ValVector<3>::GetNewVector(domain, num_components);
	shared_ptr<ValVector<3>> exact = ValVector<3>::GetNewVector(domain, num_components);
	shared_ptr<ValVector<3>> f     = ValVector<3>::GetNewVector(domain, num_components);

	// fill the vectors with some values
	DomainTools::SetValues<3>(domain, f, ffun); // fill the f vector with the associated domain using the ffun function
	DomainTools::SetValues<3>(domain, exact, gfun);

	// modify the rhs to set the boundary conditions
	// the patch operator contains some helper functions for this
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

	///////////////////
	// setup start
	///////////////////
	timer->start("Linear System Setup");

	std::shared_ptr<Operator<3>> A = patch_operator; // a PatchOperator can be passed the BiCGStab to use as an operator

	// preconditoners
	timer->start("Preconditioner Setup");

	int domain_level = 0;
	domain->setTimer(timer);
	domain_level++;

	// set the patch solver
	auto patch_bcgs = make_shared<Iterative::BiCGStab<3>>();
	// p_bcgs->setTolerance(ps_tol);
	// p_bcgs->setMaxIterations(ps_max_it);

	// create a CycleBuilder and set the options
	GMG::CycleBuilder<3> builder(copts);

	// the GMG cycle needs a vector generator for domain of each level to generate temporary work vectors
	std::shared_ptr<VectorGenerator<3>> vg(new ValVectorGenerator<3>(domain, num_components));

	// Create a smoother for the finest level
	shared_ptr<GMG::Smoother<3>> finest_smoother;
	bitset<6>                    neumann_bitset = neumann ? 0xFF : 0x0;
	if (patch_solver == "dft") {
		finest_smoother = make_shared<DFTPatchSolver<3>>(patch_operator, neumann_bitset);
	} else if (patch_solver == "fftw") {
		finest_smoother = make_shared<FFTWPatchSolver<3>>(patch_operator, neumann_bitset);
	} else {
		finest_smoother = make_shared<Iterative::PatchSolver<3>>(patch_bcgs, patch_operator);
	}

	// the next coarser domain is needed for the restrictor
	shared_ptr<Domain<3>> coarser_domain = domain_generator->getCoarserDomain();

	shared_ptr<GMG::Restrictor<3>> finest_restrictor = make_shared<GMG::LinearRestrictor<3>>(domain, coarser_domain, num_components);

	//add the finest level
	builder.addFinestLevel(patch_operator, finest_smoother, finest_restrictor, vg);

	shared_ptr<Domain<3>> finer_domain   = domain;
	shared_ptr<Domain<3>> current_domain = coarser_domain;

	//generate each of the middle levels
	while (domain_generator->hasCoarserDomain()) {
		current_domain->setTimer(timer);
		domain_level++;

		// get the coarser domain
		coarser_domain = domain_generator->getCoarserDomain();

		// vector generator
		auto middle_vector_generator = make_shared<ValVectorGenerator<3>>(current_domain, num_components);

		// create operator for middle domain
		auto middle_ghost_filler   = make_shared<TriLinearGhostFiller>(current_domain, GhostFillingType::Faces);
		auto middle_patch_operator = make_shared<StarPatchOperator<3>>(current_domain, middle_ghost_filler);

		// smoother
		shared_ptr<GMG::Smoother<3>> middle_smoother;
		if (patch_solver == "dft") {
			middle_smoother = make_shared<DFTPatchSolver<3>>(middle_patch_operator, neumann_bitset);
		} else if (patch_solver == "fftw") {
			middle_smoother = make_shared<FFTWPatchSolver<3>>(middle_patch_operator, neumann_bitset);
		} else {
			middle_smoother = make_shared<Iterative::PatchSolver<3>>(patch_bcgs, middle_patch_operator);
		}

		// restrictor and interpolator
		shared_ptr<GMG::Interpolator<3>> interpolator = make_shared<GMG::DirectInterpolator<3>>(current_domain, finer_domain, num_components);
		shared_ptr<GMG::Restrictor<3>>   restrictor   = make_shared<GMG::LinearRestrictor<3>>(current_domain, coarser_domain, num_components);

		// add the middle level
		builder.addIntermediateLevel(middle_patch_operator, middle_smoother, restrictor, interpolator, middle_vector_generator);

		finer_domain   = current_domain;
		current_domain = coarser_domain;
	}
	current_domain->setTimer(timer);

	//add the coarsest level to the builder

	// vector generator
	shared_ptr<VectorGenerator<3>> coarsest_vector_generator = make_shared<ValVectorGenerator<3>>(current_domain, 1);

	// patch operator
	shared_ptr<GhostFiller<3>>   coarsest_ghost_filler   = make_shared<TriLinearGhostFiller>(current_domain, GhostFillingType::Faces);
	shared_ptr<PatchOperator<3>> coarsest_patch_operator = make_shared<StarPatchOperator<3>>(current_domain, coarsest_ghost_filler);

	// smoother
	shared_ptr<GMG::Smoother<3>> coarsest_smoother;
	if (patch_solver == "dft") {
		coarsest_smoother = make_shared<DFTPatchSolver<3>>(coarsest_patch_operator, neumann_bitset);
	} else if (patch_solver == "fftw") {
		coarsest_smoother = make_shared<FFTWPatchSolver<3>>(coarsest_patch_operator, neumann_bitset);
	} else {
		coarsest_smoother = make_shared<Iterative::PatchSolver<3>>(patch_bcgs, coarsest_patch_operator);
	}

	// coarsets level only needs an interpolator
	shared_ptr<GMG::Interpolator<3>> interpolator = make_shared<GMG::DirectInterpolator<3>>(current_domain, finer_domain, 1);

	// add the coarsest level
	builder.addCoarsestLevel(coarsest_patch_operator, coarsest_smoother, interpolator, coarsest_vector_generator);

	// get the preconditioner operator
	std::shared_ptr<Operator<3>> M = builder.getCycle();

	timer->stop("Preconditioner Setup");

	timer->stop("Linear System Setup");

	timer->start("Linear Solve");

	// solve the system using bcgs with GMG as the preconditioner

	Iterative::BiCGStab<3> solver;
	solver.setTimer(timer);
	solver.setTolerance(tolerance);

	int num_iterations = solver.solve(vg, A, u, f, M, true);

	if (my_global_rank == 0) {
		cout << "Iterations: " << num_iterations << endl;
	}
	timer->stop("Linear Solve");

	// calculate residual
	shared_ptr<ValVector<3>> au              = ValVector<3>::GetNewVector(domain, num_components);
	shared_ptr<ValVector<3>> residual_vector = ValVector<3>::GetNewVector(domain, num_components);

	A->apply(u, au);

	residual_vector->addScaled(-1, au, 1, f);

	double residual = residual_vector->twoNorm();
	double fnorm    = f->twoNorm();

	// calculate error
	shared_ptr<ValVector<3>> error = ValVector<3>::GetNewVector(domain, 1);
	error->addScaled(-1, exact, 1, u);
	if (neumann) {
		double u_average     = domain->integrate(u) / domain->volume();
		double exact_average = domain->integrate(exact) / domain->volume();

		if (my_global_rank == 0) {
			cout << "Average of computed solution: " << u_average << endl;
			cout << "Average of exact solution: " << exact_average << endl;
		}

		error->shift(exact_average - u_average);
	}
	double error_norm     = error->twoNorm();
	double error_norm_inf = error->infNorm();
	double exact_norm     = exact->twoNorm();

	double au_sum = domain->integrate(au);
	double f_sum  = domain->integrate(f);
	if (my_global_rank == 0) {
		std::cout << std::scientific;
		std::cout.precision(13);
		std::cout << "Error (2-norm):   " << error_norm / exact_norm << endl;
		std::cout << "Error (inf-norm): " << error_norm_inf << endl;
		std::cout << "Residual: " << residual / fnorm << endl;
		std::cout << u8"ΣAu-Σf: " << au_sum - f_sum << endl;
		cout.unsetf(std::ios_base::floatfield);
		int total_cells = domain->getNumGlobalCells();
		cout << "Total cells: " << total_cells << endl;
		cout << "Number of Processors: " << num_procs << endl;
	}

#ifdef HAVE_VTK
	// output
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

	// output timer
	cout << *timer;
	MPI_Finalize();
	return 0;
}
