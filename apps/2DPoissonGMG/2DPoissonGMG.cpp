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
#include "GetP4estDomainGenerator.h"
#include <cmath>

// =========== //
// main driver //
// =========== //

using namespace std;
using namespace ThunderEgg;
using namespace ThunderEgg::Poisson;

int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);

	Communicator comm(MPI_COMM_WORLD);

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

	double tolerance = 1e-8;
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
            double y  = coord[1];
            double r2 = pow((x - x0), 2) + pow((y - y0), 2);
            return exp(-alpha / 2.0 * r2) * (pow(alpha, 2) * r2 - 2 * alpha);
		};
		gfun = [&](const array<double, 2> &coord) {
			double x  = coord[0];
			double y  = coord[1];
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
			double y = coord[1];
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
			double y = coord[1];
			return exp(cos(10 * M_PI * x)) - exp(cos(11 * M_PI * y));
		};
		ffun = [](const array<double, 2> &coord) {
			double x = coord[0];
			double y = coord[1];
			return 100 * M_PI * M_PI * (pow(sin(10 * M_PI * x), 2) - cos(10 * M_PI * x))
			       * exp(cos(10 * M_PI * x))
			       + 121 * M_PI * M_PI * (cos(11 * M_PI * y) - pow(sin(11 * M_PI * y), 2))
			         * exp(cos(11 * M_PI * y));
		};
		nfunx = [](const array<double, 2> &coord) {
			double x = coord[0];
			return -10 * M_PI * sin(10 * M_PI * x) * exp(cos(10 * M_PI * x));
		};

		nfuny = [](const array<double, 2> &coord) {
			double y = coord[2];
			return 11 * M_PI * sin(11 * M_PI * y) * exp(cos(11 * M_PI * y));
		};
	} else {
		ffun = [](const array<double, 2> &coord) {
			double x = coord[0];
			double y = coord[1];
			return -5 * M_PI * M_PI * sin(M_PI * y) * cos(2 * M_PI * x);
		};
		gfun = [](const array<double, 2> &coord) {
			double x = coord[0];
			double y = coord[1];
			return sin(M_PI * y) * cos(2 * M_PI * x);
		};
		nfunx = [](const array<double, 2> &coord) {
			double x = coord[0];
			double y = coord[1];
			return -2 * M_PI * sin(M_PI * y) * sin(2 * M_PI * x);
		};
		nfuny = [](const array<double, 2> &coord) {
			double x = coord[0];
			double y = coord[1];
			return M_PI * cos(M_PI * y) * cos(2 * M_PI * x);
		};
	}

	// create a Timer object to keep track of timings
	std::shared_ptr<Timer> timer = make_shared<Timer>(comm);

	auto bmf = [](int block_no, double unit_x, double unit_y, double &x, double &y) {
		x = unit_x;
		y = unit_y;
	};

	// Set the number of cells in the x, y and z direction.
	std::array<int, 2> ns;
	ns[0] = n;
	ns[1] = n;

	// A DomainGenerator will create domains for the Multigrid algorithm from a tree
	int num_ghost_cells = 1; // the poission operator needs 1 row/column of ghost cells on the edges of a patch

	p4est_connectivity * conn             = p4est_connectivity_new_brick(1, 1, false, false);
	P4estDomainGenerator domain_generator = getP4estDomainGenerator(conn, mesh_filename, ns, num_ghost_cells, bmf);

	// Get the finest domain from the tree
	Domain<2> domain = domain_generator.getFinestDomain();

	// A patch operator needs a GhostFiller object to define how to fill ghost cells for the patches.
	// This one will use a tri-linear interpolation scheme at the refinement boundarys of the domain
	BiLinearGhostFiller ghost_filler(domain, GhostFillingType::Faces);

	// create patch operator that uses a typical 2nd order 7 point poisson stencil
	StarPatchOperator<2> patch_operator(domain, ghost_filler, neumann);

	timer->start("Domain Initialization");

	// Create some new vectors for the domain
	int num_components = 1; //the poisson operator just has one value in each cell

	Vector<2> exact(domain, num_components);
	Vector<2> f(domain, num_components);

	// fill the vectors with some values
	DomainTools::SetValues<2>(domain, f, ffun); // fill the f vector with the associated domain using the ffun function
	DomainTools::SetValues<2>(domain, exact, gfun);

	// modify the rhs to set the boundary conditions
	// the patch operator contains some helper functions for this
	if (neumann) {
		patch_operator.addNeumannBCToRHS(f, gfun, {nfunx, nfuny});
	} else {
		patch_operator.addDrichletBCToRHS(f, gfun);
	}

	timer->stop("Domain Initialization");

	if (neumann && !no_zero_rhs_avg) {
		double fdiff = DomainTools::Integrate<2>(domain, f) / domain.volume();
		if (comm.getRank() == 0)
			cout << "Fdiff: " << fdiff << endl;
		f.shift(-fdiff);
	}

	///////////////////
	// setup start
	///////////////////
	timer->start("Linear System Setup");

	// preconditoners
	timer->start("Preconditioner Setup");

	int domain_level = 0;
	domain.setTimer(timer);
	domain_level++;

	// set the patch solver
	Iterative::BiCGStab<2> patch_bcgs;
	//p_bcgs->setTolerance(ps_tol);
	// p_bcgs->setMaxIterations(ps_max_it);

	// create a CycleBuilder and set the options
	GMG::CycleBuilder<2> builder(copts);

	// Create a smoother for the finest level
	shared_ptr<PatchSolver<2>> finest_smoother;
	bitset<4>                  neumann_bitset = neumann ? 0xF : 0x0;
	if (patch_solver == "dft") {
		finest_smoother.reset(new DFTPatchSolver<2>(patch_operator, neumann_bitset));
	} else if (patch_solver == "fftw") {
		finest_smoother.reset(new FFTWPatchSolver<2>(patch_operator, neumann_bitset));
	} else {
		finest_smoother.reset(new Iterative::PatchSolver<2>(patch_bcgs, patch_operator));
	}

	std::shared_ptr<Operator<2>> M;
	if (true) {
		// the next coarser domain is needed for the restrictor
		Domain<2> coarser_domain = domain_generator.getCoarserDomain();

		GMG::LinearRestrictor<2> finest_restrictor(domain, coarser_domain);

		//add the finest level
		builder.addFinestLevel(patch_operator, *finest_smoother, finest_restrictor);

		Domain<2> finer_domain   = domain;
		Domain<2> current_domain = coarser_domain;

		//generate each of the middle levels
		while (domain_generator.hasCoarserDomain()) {
			current_domain.setTimer(timer);
			domain_level++;

			// get the coarser domain
			coarser_domain = domain_generator.getCoarserDomain();

			// create operator for middle domain
			BiLinearGhostFiller  middle_ghost_filler(current_domain, GhostFillingType::Faces);
			StarPatchOperator<2> middle_patch_operator(current_domain, middle_ghost_filler);

			// smoother
			unique_ptr<GMG::Smoother<2>> middle_smoother;
			if (patch_solver == "dft") {
				middle_smoother.reset(new DFTPatchSolver<2>(middle_patch_operator, neumann_bitset));
			} else if (patch_solver == "fftw") {
				middle_smoother.reset(new FFTWPatchSolver<2>(middle_patch_operator, neumann_bitset));
			} else {
				middle_smoother.reset(new Iterative::PatchSolver<2>(patch_bcgs, middle_patch_operator));
			}

			// restrictor and interpolator
			GMG::DirectInterpolator<2> interpolator(current_domain, finer_domain);
			GMG::LinearRestrictor<2>   restrictor(current_domain, coarser_domain);

			// add the middle level
			builder.addIntermediateLevel(middle_patch_operator, *middle_smoother, restrictor, interpolator);

			finer_domain   = current_domain;
			current_domain = coarser_domain;
		}
		current_domain.setTimer(timer);

		//add the coarsest level to the builder

		// patch operator
		BiLinearGhostFiller  coarsest_ghost_filler(current_domain, GhostFillingType::Faces);
		StarPatchOperator<2> coarsest_patch_operator(current_domain, coarsest_ghost_filler);

		// smoother
		unique_ptr<GMG::Smoother<2>> coarsest_smoother;
		if (patch_solver == "dft") {
			coarsest_smoother.reset(new DFTPatchSolver<2>(coarsest_patch_operator, neumann_bitset));
		} else if (patch_solver == "fftw") {
			coarsest_smoother.reset(new FFTWPatchSolver<2>(coarsest_patch_operator, neumann_bitset));
		} else {
			coarsest_smoother.reset(new Iterative::PatchSolver<2>(patch_bcgs, coarsest_patch_operator));
		}

		// coarsets level only needs an interpolator
		GMG::DirectInterpolator<2> interpolator(current_domain, finer_domain);

		// add the coarsest level
		builder.addCoarsestLevel(coarsest_patch_operator, *coarsest_smoother, interpolator);

		// get the preconditioner operator
		M = builder.getCycle();
	} else {
		M = finest_smoother;
	}

	timer->stop("Preconditioner Setup");

	timer->stop("Linear System Setup");

	timer->start("Linear Solve");

	// solve the system using bcgs with GMG as the preconditioner

	Iterative::BiCGStab<2> solver;
	solver.setTimer(timer);
	solver.setTolerance(tolerance);

	Vector<2> u(domain, num_components);

	int num_iterations = solver.solve(patch_operator, u, f, M.get(), true);

	if (comm.getRank() == 0) {
		cout << "Iterations: " << num_iterations << endl;
	}
	timer->stop("Linear Solve");

	// calculate residual
	Vector<2> au(domain, num_components);
	Vector<2> residual_vector(domain, num_components);

	patch_operator.apply(u, au);

	residual_vector.addScaled(-1, au, 1, f);

	double residual = residual_vector.twoNorm();
	double fnorm    = f.twoNorm();

	// calculate error
	Vector<2> error(domain, 1);
	error.addScaled(-1, exact, 1, u);
	if (neumann) {
		double u_average     = DomainTools::Integrate<2>(domain, u) / domain.volume();
		double exact_average = DomainTools::Integrate<2>(domain, exact) / domain.volume();

		if (comm.getRank() == 0) {
			cout << "Average of computed solution: " << u_average << endl;
			cout << "Average of exact solution: " << exact_average << endl;
		}

		error.shift(exact_average - u_average);
	}
	double error_norm     = error.twoNorm();
	double error_norm_inf = error.infNorm();
	double exact_norm     = exact.twoNorm();

	double au_sum = DomainTools::Integrate<2>(domain, au);
	double f_sum  = DomainTools::Integrate<2>(domain, f);
	if (comm.getRank() == 0) {
		std::cout << std::scientific;
		std::cout.precision(13);
		std::cout << "Error (2-norm):   " << error_norm / exact_norm << endl;
		std::cout << "Error (inf-norm): " << error_norm_inf << endl;
		std::cout << "Residual: " << residual / fnorm << endl;
		std::cout << u8"ΣAu-Σf: " << au_sum - f_sum << endl;
		cout.unsetf(std::ios_base::floatfield);
		int total_cells = domain.getNumGlobalCells();
		cout << "Total cells: " << total_cells << endl;
		cout << "Number of Processors: " << comm.getSize() << endl;
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
