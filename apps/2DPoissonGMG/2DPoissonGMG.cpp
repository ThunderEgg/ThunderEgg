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
#include "TreeToP4est.h"
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
	CLI::App app{"ThunderEgg 2d poisson solver example"};

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
	app.add_set_ignore_case("--problem", problem, {"trig", "gauss", "circle", "trig_gauss", "zero"},
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
	function<double(const std::array<double, 2> &)> ffun; // rhs
	function<double(const std::array<double, 2> &)> gfun; // lhs
	// functions needed for neumann boundary conditions
	function<double(const std::array<double, 2> &)> nfunx; // lhs derivative in x direction
	function<double(const std::array<double, 2> &)> nfuny; // lhs derivative in y direction

	if (problem == "gauss") {
		double x0    = .5;
		double y0    = .5;
		double alpha = 1000;
		ffun         = [&](const array<double, 2> &coord) {
            double x  = coord[0];
            double y  = coord[2];
            double r2 = pow((x - x0), 2) + pow((y - y0), 2);
            return exp(-alpha / 2.0 * r2) * (pow(alpha, 2) * r2 - 2 * alpha);
		};
		gfun = [&](const array<double, 2> &coord) {
			double x  = coord[0];
			double y  = coord[2];
			double r2 = pow((x - x0), 2) + pow((y - y0), 2);
			return exp(-alpha / 2.0 * r2);
		};
		nfunx = [](const array<double, 2> &coord) { return 0; };
		nfuny = [](const array<double, 2> &coord) { return 0; };
	} else if (problem == "zero") {
		ffun  = [](const array<double, 2> &coord) { return 0; };
		gfun  = [](const array<double, 2> &coord) { return 0; };
		nfunx = [](const array<double, 2> &coord) { return 0; };
		nfuny = [](const array<double, 2> &coord) { return 0; };
	} else if (problem == "circle") {
		ffun = [](const array<double, 2> &coord) {
			double x = coord[0];
			double y = coord[2];
			double xdist, ydist, dist;
			// distance form center circle
			xdist = x - 0.5;
			ydist = y - 0.5;
			dist  = sqrt(xdist * xdist + ydist * ydist);
			if (dist < 0.2) {
				return 1;
			}
			for (int i = 0; i < 4; i++) {
				// larger circles
				double theta = i * M_PI / 2.0;
				xdist        = x - (0.3 * cos(theta) + 0.5);
				ydist        = y - (0.3 * sin(theta) + 0.5);
				dist         = sqrt(xdist * xdist + ydist * ydist);
				if (dist < 0.1) {
					return 1;
				}
				// smaller circles
				theta = M_PI / 4.0 + i * M_PI / 2.0;
				xdist = x - (0.275 * cos(theta) + 0.5);
				ydist = y - (0.275 * sin(theta) + 0.5);
				dist  = sqrt(xdist * xdist + ydist * ydist);
				if (dist < 0.075) {
					return 1;
				}
			}
			return 0;
		};
		gfun  = [](const array<double, 2> &coord) { return 0; };
		nfunx = [](const array<double, 2> &coord) { return 0; };
		nfuny = [](const array<double, 2> &coord) { return 0; };
	} else if (problem == "trig_gauss") {
		gfun = [](const array<double, 2> &coord) {
			double x = coord[0];
			double y = coord[2];
			return exp(cos(10 * M_PI * x)) - exp(cos(11 * M_PI * y));
		};
		ffun = [](const array<double, 2> &coord) {
			double x = coord[0];
			double y = coord[2];
			return 100 * M_PI * M_PI * (pow(sin(10 * M_PI * x), 2) - cos(10 * M_PI * x))
			       * exp(cos(10 * M_PI * x))
			       + 121 * M_PI * M_PI * (cos(11 * M_PI * y) - pow(sin(11 * M_PI * y), 2))
			         * exp(cos(11 * M_PI * y));
		};
		nfunx = [](const array<double, 2> &coord) {
			double x = coord[0];
			double y = coord[2];
			return -10 * M_PI * sin(10 * M_PI * x) * exp(cos(10 * M_PI * x));
		};

		nfuny = [](const array<double, 2> &coord) {
			double x = coord[0];
			double y = coord[2];
			return 11 * M_PI * sin(11 * M_PI * y) * exp(cos(11 * M_PI * y));
		};
	} else {
		ffun = [](const array<double, 2> &coord) {
			double x = coord[0];
			double y = coord[2];
			return -5 * M_PI * M_PI * sinl(M_PI * y) * cosl(2 * M_PI * x);
		};
		gfun = [](const array<double, 2> &coord) {
			double x = coord[0];
			double y = coord[2];
			return sinl(M_PI * y) * cosl(2 * M_PI * x);
		};
		nfunx = [](const array<double, 2> &coord) {
			double x = coord[0];
			double y = coord[2];
			return -2 * M_PI * sinl(M_PI * y) * sinl(2 * M_PI * x);
		};
		nfuny = [](const array<double, 2> &coord) {
			double x = coord[0];
			double y = coord[2];
			return M_PI * cosl(M_PI * y) * cosl(2 * M_PI * x);
		};
	}

	// create a Timer object to keep track of timings
	std::shared_ptr<Timer> timer = make_shared<Timer>(MPI_COMM_WORLD);

	// read in a tree from file
	Tree<2> t;
	t = Tree<2>(mesh_filename);
	for (int i = 0; i < div; i++) {
		t.refineLeaves();
	}
	TreeToP4est ttp(t);

	auto bmf = [](int block_no, double unit_x, double unit_y, double &x, double &y) {
		x = unit_x;
		y = unit_y;
	};
	auto inf
	= [=](Side<2> s, const array<double, 2> &, const array<double, 2> &) { return neumann; };

	// Set the number of cells in the x, y and z direction.
	std::array<int, 2> ns;
	ns[0] = n;
	ns[1] = n;

	// A DomainGenerator will create domains for the Multigrid algorithm from a tree
	int                            num_ghost_cells  = 1; // the poission operator needs 1 row/column of ghost cells on the edges of a patch
	shared_ptr<DomainGenerator<2>> domain_generator = make_shared<P4estDomGen>(ttp.p4est, ns, 1, inf, bmf);

	// Get the finest domain from the tree
	shared_ptr<Domain<2>> domain = domain_generator->getFinestDomain();

	// A patch operator needs a GhostFiller object to define how to fill ghost cells for the patches.
	// This one will use a tri-linear interpolation scheme at the refinement boundarys of the domain
	shared_ptr<GhostFiller<2>> ghost_filler = make_shared<BiLinearGhostFiller>(domain);

	// create patch operator that uses a typical 2nd order 7 point poisson stencil
	shared_ptr<StarPatchOperator<2>> patch_operator = make_shared<StarPatchOperator<2>>(domain, ghost_filler, neumann);

	timer->start("Domain Initialization");

	// Create some new vectors for the domain
	int num_components = 1; //the poisson operator just has one value in each cell

	shared_ptr<ValVector<2>> u     = ValVector<2>::GetNewVector(domain, num_components);
	shared_ptr<ValVector<2>> exact = ValVector<2>::GetNewVector(domain, num_components);
	shared_ptr<ValVector<2>> f     = ValVector<2>::GetNewVector(domain, num_components);

	// fill the vectors with some values
	DomainTools::SetValues<2>(domain, f, ffun); // fill the f vector with the associated domain using the ffun function
	DomainTools::SetValues<2>(domain, exact, gfun);

	// modify the rhs to set the boundary conditions
	// the patch operator contains some helper functions for this
	if (neumann) {
		patch_operator->addDrichletBCToRHS(f, gfun);
	} else {
		patch_operator->addNeumannBCToRHS(f, gfun, {nfunx, nfuny});
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

	std::shared_ptr<Operator<2>> A = patch_operator; // a PatchOperator can be passed the BiCGStab to use as an operator

	// preconditoners
	timer->start("Preconditioner Setup");

	int domain_level = 0;
	domain->setId(domain_level);
	domain->setTimer(timer);
	domain_level++;

	// set the patch solver
	auto patch_bcgs = make_shared<Iterative::BiCGStab<2>>();
	// p_bcgs->setTolerance(ps_tol);
	// p_bcgs->setMaxIterations(ps_max_it);

	// create a CycleBuilder and set the options
	GMG::CycleBuilder<2> builder(copts);

	// the GMG cycle needs a vector generator for domain of each level to generate temporary work vectors
	std::shared_ptr<VectorGenerator<2>> vg(new ValVectorGenerator<2>(domain, num_components));

	// Create a smoother for the finest level
	shared_ptr<GMG::Smoother<2>> finest_smoother;
	bitset<4>                    neumann_bitset = neumann ? 0xF : 0x0;
	if (patch_solver == "dft") {
		finest_smoother = make_shared<DFTPatchSolver<2>>(patch_operator, neumann_bitset);
	} else if (patch_solver == "fftw") {
		finest_smoother = make_shared<FFTWPatchSolver<2>>(patch_operator, neumann_bitset);
	} else {
		finest_smoother = make_shared<Iterative::PatchSolver<2>>(patch_bcgs, patch_operator);
	}

	// the next coarser domain is needed for the restrictor
	shared_ptr<Domain<2>> coarser_domain = domain_generator->getCoarserDomain();

	shared_ptr<GMG::Restrictor<2>> finest_restrictor = make_shared<GMG::LinearRestrictor<2>>(domain, coarser_domain, num_components);

	//add the finest level
	builder.addFinestLevel(patch_operator, finest_smoother, finest_restrictor, vg);

	shared_ptr<Domain<2>> finer_domain   = domain;
	shared_ptr<Domain<2>> current_domain = coarser_domain;

	//generate each of the middle levels
	while (domain_generator->hasCoarserDomain()) {
		current_domain->setId(domain_level);
		current_domain->setTimer(timer);
		domain_level++;

		// get the coarser domain
		coarser_domain = domain_generator->getCoarserDomain();

		// vector generator
		auto middle_vector_generator = make_shared<ValVectorGenerator<2>>(current_domain, num_components);

		// create operator for middle domain
		auto middle_ghost_filler   = make_shared<BiLinearGhostFiller>(current_domain);
		auto middle_patch_operator = make_shared<StarPatchOperator<2>>(current_domain, middle_ghost_filler);

		// smoother
		shared_ptr<GMG::Smoother<2>> middle_smoother;
		if (patch_solver == "dft") {
			middle_smoother = make_shared<DFTPatchSolver<2>>(middle_patch_operator, neumann_bitset);
		} else if (patch_solver == "fftw") {
			middle_smoother = make_shared<FFTWPatchSolver<2>>(middle_patch_operator, neumann_bitset);
		} else {
			middle_smoother = make_shared<Iterative::PatchSolver<2>>(patch_bcgs, middle_patch_operator);
		}

		// restrictor and interpolator
		shared_ptr<GMG::Interpolator<2>> interpolator = make_shared<GMG::DirectInterpolator<2>>(current_domain, finer_domain, num_components);
		shared_ptr<GMG::Restrictor<2>>   restrictor   = make_shared<GMG::LinearRestrictor<2>>(current_domain, coarser_domain, num_components);

		// add the middle level
		builder.addIntermediateLevel(middle_patch_operator, middle_smoother, restrictor, interpolator, middle_vector_generator);

		finer_domain   = current_domain;
		current_domain = coarser_domain;
	}
	current_domain->setId(domain_level);
	current_domain->setTimer(timer);

	//add the coarsest level to the builder

	// vector generator
	shared_ptr<VectorGenerator<2>> coarsest_vector_generator = make_shared<ValVectorGenerator<2>>(current_domain, 1);

	// patch operator
	shared_ptr<GhostFiller<2>>   coarsest_ghost_filler   = make_shared<BiLinearGhostFiller>(current_domain);
	shared_ptr<PatchOperator<2>> coarsest_patch_operator = make_shared<StarPatchOperator<2>>(current_domain, coarsest_ghost_filler);

	// smoother
	shared_ptr<GMG::Smoother<2>> coarsest_smoother;
	if (patch_solver == "dft") {
		coarsest_smoother = make_shared<DFTPatchSolver<2>>(coarsest_patch_operator, neumann_bitset);
	} else if (patch_solver == "fftw") {
		coarsest_smoother = make_shared<FFTWPatchSolver<2>>(coarsest_patch_operator, neumann_bitset);
	} else {
		coarsest_smoother = make_shared<Iterative::PatchSolver<2>>(patch_bcgs, coarsest_patch_operator);
	}

	// coarsets level only needs an interpolator
	shared_ptr<GMG::Interpolator<2>> interpolator = make_shared<GMG::DirectInterpolator<2>>(current_domain, finer_domain, 1);

	// add the coarsest level
	builder.addCoarsestLevel(coarsest_patch_operator, coarsest_smoother, interpolator, coarsest_vector_generator);

	// get the preconditioner operator
	std::shared_ptr<Operator<2>> M = builder.getCycle();

	timer->stop("Preconditioner Setup");

	timer->stop("Linear System Setup");

	timer->start("Linear Solve");

	// solve the system using bcgs with GMG as the preconditioner

	Iterative::BiCGStab<2> solver;
	solver.setTimer(timer);
	solver.setTolerance(tolerance);

	int num_iterations = solver.solve(vg, A, u, f, M, true);

	if (my_global_rank == 0) {
		cout << "Iterations: " << num_iterations << endl;
	}
	timer->stop("Linear Solve");

	// calculate residual
	shared_ptr<ValVector<2>> au              = ValVector<2>::GetNewVector(domain, num_components);
	shared_ptr<ValVector<2>> residual_vector = ValVector<2>::GetNewVector(domain, num_components);

	A->apply(u, au);

	residual_vector->addScaled(-1, au, 1, f);

	double residual = residual_vector->twoNorm();
	double fnorm    = f->twoNorm();

	// calculate error
	shared_ptr<ValVector<2>> error = ValVector<2>::GetNewVector(domain, 1);
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
