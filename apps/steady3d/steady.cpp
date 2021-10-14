/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured 
 *  Cartesian grids.
 *
 *  Copyright (c) 2017-2020 Scott Aiton
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
#include "Init.h"
#include "Writers/ClawWriter.h"
#include <ThunderEgg/BiCGStab.h>
#include <ThunderEgg/Domain.h>
#include <ThunderEgg/DomainWrapOp.h>
#include <ThunderEgg/Experimental/DomGen.h>
#include <ThunderEgg/Experimental/OctTree.h>
#include <ThunderEgg/GMG/CycleFactory3d.h>
#include <ThunderEgg/GMG/CycleOpts.h>
#include <ThunderEgg/MatWrapper.h>
#include <ThunderEgg/PetscShellCreator.h>
#include <ThunderEgg/Poisson/MatrixHelper.h>
#include <ThunderEgg/Poisson/Schur/DftPatchSolver.h>
#include <ThunderEgg/Poisson/Schur/FftwPatchSolver.h>
#include <ThunderEgg/Poisson/Schur/SevenPtPatchOperator.h>
#include <ThunderEgg/Poisson/Schur/StarPatchOperator.h>
#include <ThunderEgg/Schur/InterfaceDomain.h>
#include <ThunderEgg/Schur/PolyChebPrec.h>
#include <ThunderEgg/Schur/SchurMatrixHelper.h>
#include <ThunderEgg/Schur/SchurWrapOp.h>
#include <ThunderEgg/Schur/TriLinInterp.h>
#include <ThunderEgg/SchwarzPrec.h>
#include <ThunderEgg/Timer.h>
#ifdef ENABLE_AMGX
#include "AmgxWrapper.h"
#endif
#ifdef ENABLE_MUELU
#include <MueLuWrapper.h>
#endif
#ifdef ENABLE_MUELU_CUDA
#include <MueLuCudaWrapper.h>
#endif
#ifdef HAVE_VTK
#include "Writers/VtkWriter.h"
#endif
#include <cmath>
#include <iostream>
#include <memory>
#include <petscksp.h>
#include <petscsys.h>
#include <petscvec.h>
#include <petscviewer.h>
#include <string>
#include <unistd.h>

// =========== //
// main driver //
// =========== //

using namespace std;
using namespace ThunderEgg;
using namespace ThunderEgg::Schur;
using namespace ThunderEgg::Experimental;
using namespace ThunderEgg::Poisson;
using namespace ThunderEgg::Poisson::Schur;

int main(int argc, char *argv[])
{
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
	app.add_set_ignore_case("--matrix_type", matrix_type, {"wrap", "crs", "pbm"},
	                        "Which type of matrix operator to use");

	int div = 0;
	app.add_option("--divide", div, "Number of levels to add to octtree");

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
	app.add_set_ignore_case("--problem", problem, {"trig", "gauss", "zero"},
	                        "Which problem to solve");

	string solver_type = "thunderegg";
	app.add_set_ignore_case("--solver", solver_type, {"petsc", "thunderegg"},
	                        "Which Solver to use");
	string preconditioner = "";
	app.add_set_ignore_case("--prec", preconditioner, {"GMG", "Schwarz", "BlockJacobi", "cheb"},
	                        "Which Preconditoner to use");

	string patch_solver = "fftw";
	app.add_option("--patch_solver", patch_solver, "Which patch solver to use");

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
	std::array<int, 3> ns;
	ns.fill(n);

	///////////////
	// Create Mesh
	///////////////
	shared_ptr<Domain<3>> dc;
	Tree<3>               t;
	t = Tree<3>(mesh_filename);
	for (int i = 0; i < div; i++) {
		t.refineLeaves();
	}

	// the functions that we are using
	function<double(double, double, double)> ffun;
	function<double(double, double, double)> gfun;
	function<double(double, double, double)> nfunx;
	function<double(double, double, double)> nfuny;
	function<double(double, double, double)> nfunz;

	if (problem == "zero") {
		ffun  = [](double x, double y, double z) { return 0; };
		gfun  = [](double x, double y, double z) { return 0; };
		nfunx = [](double x, double y, double z) { return 0; };
		nfuny = [](double x, double y, double z) { return 0; };
		nfunz = [](double x, double y, double z) { return 0; };
	} else if (problem == "gauss") {
		gfun = [](double x, double y, double z) {
			return exp(cos(10 * M_PI * x)) - exp(cos(11 * M_PI * y)) + exp(cos(12 * M_PI * z));
		};
		ffun = [](double x, double y, double z) {
			return -M_PI * M_PI
			       * (100 * exp(cos(10 * M_PI * x)) * cos(10 * M_PI * x)
			          - 100 * exp(cos(10 * M_PI * x)) * pow(sin(10 * M_PI * x), 2)
			          - 121 * exp(cos(11 * M_PI * y)) * cos(11 * M_PI * y)
			          + 121 * exp(cos(11 * M_PI * y)) * pow(sin(11 * M_PI * y), 2)
			          + 144 * exp(cos(12 * M_PI * z)) * cos(12 * M_PI * z)
			          - 144 * exp(cos(12 * M_PI * z)) * pow(sin(12 * M_PI * z), 2));
		};
		nfunx = [](double x, double y, double z) {
			return -10 * M_PI * sin(10 * M_PI * x) * exp(cos(10 * M_PI * x));
		};
		nfuny = [](double x, double y, double z) {
			return 11 * M_PI * sin(11 * M_PI * y) * exp(cos(11 * M_PI * y));
		};
		nfunz = [](double x, double y, double z) {
			return -12 * M_PI * sin(12 * M_PI * z) * exp(cos(12 * M_PI * z));
		};
	} else {
		ffun = [](double x, double y, double z) {
			x += .3;
			y += .3;
			z += .3;
			return -77.0 / 36 * M_PI * M_PI * sin(M_PI * x) * cos(2.0 / 3 * M_PI * y)
			       * sin(5.0 / 6 * M_PI * z);
		};
		gfun = [](double x, double y, double z) {
			x += .3;
			y += .3;
			z += .3;
			return sin(M_PI * x) * cos(2.0 / 3 * M_PI * y) * sin(5.0 / 6 * M_PI * z);
		};
		nfunx = [](double x, double y, double z) {
			x += .3;
			y += .3;
			z += .3;
			return M_PI * cos(M_PI * x) * cos(2.0 / 3 * M_PI * y) * sin(5.0 / 6 * M_PI * z);
		};
		nfuny = [](double x, double y, double z) {
			x += .3;
			y += .3;
			z += .3;
			return -2.0 / 3 * M_PI * sin(M_PI * x) * sin(2.0 / 3 * M_PI * y)
			       * sin(5.0 / 6 * M_PI * z);
		};
		nfunz = [](double x, double y, double z) {
			x += .3;
			y += .3;
			z += .3;
			return 5.0 / 6 * M_PI * sin(M_PI * x) * cos(2.0 / 3 * M_PI * y)
			       * cos(5.0 / 6 * M_PI * z);
		};
	}

	shared_ptr<DomainGenerator<3>> dcg(new DomGen<3>(t, ns, neumann));

	dc = dcg->getFinestDomain();
	shared_ptr<InterfaceDomain<3>> sch(new InterfaceDomain<3>(dc));

	// set the patch solver
	shared_ptr<PatchSolver<3>> p_solver;

	// patch operator
	shared_ptr<PatchOperator<3>> p_operator(new StarPatchOperator<3>(sch));

	// interface interpolator
	shared_ptr<IfaceInterp<3>> p_interp(new TriLinInterp(sch));

	if (patch_solver == "dft") {
		p_solver.reset(new DftPatchSolver<3>(sch));
	} else {
		p_solver.reset(new FftwPatchSolver<3>(sch));
	}

	Timer timer;
	for (int loop = 0; loop < loop_count; loop++) {
		timer.start("Domain Initialization");

		// Initialize Vectors
		shared_ptr<PetscVector<3>> u     = PetscVector<3>::GetNewVector(dc);
		shared_ptr<PetscVector<3>> exact = PetscVector<3>::GetNewVector(dc);
		shared_ptr<PetscVector<3>> f     = PetscVector<3>::GetNewVector(dc);
		shared_ptr<PetscVector<3>> au    = PetscVector<3>::GetNewVector(dc);

		if (neumann) {
			Init::initNeumann(*dc, f->vec, exact->vec, ffun, gfun, nfunx, nfuny, nfunz);
		} else {
			Init::initDirichlet(*dc, f->vec, exact->vec, ffun, gfun);
		}

		timer.stop("Domain Initialization");

		shared_ptr<PetscVector<2>> gamma = sch->getNewGlobalInterfaceVector();
		shared_ptr<PetscVector<2>> diff  = sch->getNewGlobalInterfaceVector();
		shared_ptr<PetscVector<2>> b     = sch->getNewGlobalInterfaceVector();

		if (neumann && !no_zero_rhs_avg) {
			double fdiff = dc->integrate(f) / dc->volume();
			if (my_global_rank == 0)
				cout << "Fdiff: " << fdiff << endl;
			f->shift(-fdiff);
		}

		if (solve_schur) {
			// initialize vectors for Schur compliment system

			// Get the b vector
			gamma->set(0);
			p_solver->solve(f, u, gamma);
			p_interp->interpolateToInterface(u, b);
			b->addScaled(-1, gamma);
			b->scale(-1);

			if (rhs_filename != "") {
				PetscViewer viewer;
				PetscViewerBinaryOpen(PETSC_COMM_WORLD, rhs_filename.c_str(), FILE_MODE_WRITE,
				                      &viewer);
				VecView(b->vec, viewer);
				PetscViewerDestroy(&viewer);
			}
			std::shared_ptr<Operator<2>> A;
			PW<Mat>                      A_petsc;
			std::shared_ptr<Operator<2>> M;
			PW<PC>                       M_petsc;

			///////////////////
			// setup start
			///////////////////
			timer.start("Linear System Setup");

			timer.start("Matrix Formation");
			if (matrix_type == "wrap") {
				A.reset(new SchurWrapOp<3>(dc, sch, p_solver, p_interp));
			} else if (matrix_type == "crs") {
				SchurMatrixHelper smh(sch, p_solver, p_interp);
				A_petsc = smh.formCRSMatrix();
				A.reset(new MatWrapper<2>(A_petsc));
				if (setrow) {
					int row = 0;
					MatZeroRows(A_petsc, 1, &row, 1.0, nullptr, nullptr);
				}

				if (matrix_filename != "") {
					PetscViewer viewer;
					PetscViewerBinaryOpen(PETSC_COMM_WORLD, matrix_filename.c_str(),
					                      FILE_MODE_WRITE, &viewer);
					MatView(A_petsc, viewer);
					PetscViewerDestroy(&viewer);
				}
			} else if (matrix_type == "pbm") {
				SchurMatrixHelper smh(sch, p_solver, p_interp);
				A_petsc = smh.getPBMatrix();
				A.reset(new MatWrapper<2>(A_petsc));
			}
			timer.stop("Matrix Formation");
			// preconditoners
			timer.start("Preconditioner Setup");
			if (preconditioner == "Scwharz") {
				throw 3;
			} else if (preconditioner == "GMG") {
				throw 3;
			} else if (preconditioner == "cheb") {
				M.reset(new PolyChebPrec(dc, sch, p_interp, p_solver));
			}
			timer.stop("Preconditioner Setup");

			PW<KSP> solver;
			// setup petsc if needed
			if (solver_type == "petsc") {
				timer.start("Petsc Setup");

				KSPCreate(MPI_COMM_WORLD, &solver);
				KSPSetFromOptions(solver);
				KSPSetOperators(solver, A_petsc, A_petsc);
				if (M != nullptr) {
					PC M_petsc;
					KSPGetPC(solver, &M_petsc);
					PetscShellCreator::getPCShell(M_petsc, M, sch);
				}
				KSPSetUp(solver);
				KSPSetTolerances(solver, tolerance, PETSC_DEFAULT, PETSC_DEFAULT, 5000);

				timer.stop("Petsc Setup");
			}
			///////////////////
			// setup end
			///////////////////
			timer.stop("Linear System Setup");

			timer.start("Linear Solve");
			if (solver_type == "petsc") {
				KSPSolve(solver, b->vec, gamma->vec);
				int its;
				KSPGetIterationNumber(solver, &its);
				if (my_global_rank == 0) {
					cout << "Iterations: " << its << endl;
				}
			} else {
				std::shared_ptr<VectorGenerator<2>> vg(new InterfaceDomainVG<2>(sch));

				int its = BiCGStab<2>::solve(vg, A, gamma, b, M);
				if (my_global_rank == 0) {
					cout << "Iterations: " << its << endl;
				}
			}
			timer.stop("Linear Solve");

			// Do one last solve
			timer.start("Patch Solve");

			p_solver->solve(f, u, gamma);

			timer.stop("Patch Solve");

			p_operator->apply(u, gamma, au);
		} else {
			std::shared_ptr<Operator<3>> A;
			PW<Mat>                      A_petsc;
			std::shared_ptr<Operator<3>> M;
			///////////////////
			// setup start
			///////////////////
			timer.start("Linear System Setup");

			timer.start("Matrix Formation");
			if (matrix_type == "wrap") {
				A.reset(new DomainWrapOp<3>(sch, p_interp, p_operator));
				A_petsc = PetscShellCreator::getMatShell(A, dc);
			} else if (matrix_type == "crs") {
				MatrixHelper mh(dc);
				A_petsc = mh.formCRSMatrix();
				A.reset(new MatWrapper<3>(A_petsc));
				if (setrow) {
					int row = 0;
					MatZeroRows(A_petsc, 1, &row, 1.0, nullptr, nullptr);
				}

				if (matrix_filename != "") {
					PetscViewer viewer;
					PetscViewerBinaryOpen(PETSC_COMM_WORLD, matrix_filename.c_str(),
					                      FILE_MODE_WRITE, &viewer);
					MatView(A_petsc, viewer);
					PetscViewerDestroy(&viewer);
				}
			} else if (matrix_type == "pbm") {
				throw 3;
			}
			timer.stop("Matrix Formation");
			// preconditoners
			timer.start("Preconditioner Setup");
			if (preconditioner == "Scwharz") {
				M.reset(new SchwarzPrec<3>(sch, p_solver, p_interp));
			} else if (preconditioner == "GMG") {
				timer.start("GMG Setup");

				M = GMG::CycleFactory3d::getCycle(copts, dcg, p_solver, p_operator, p_interp);

				timer.stop("GMG Setup");
			} else if (preconditioner == "cheb") {
				throw 3;
			}
			timer.stop("Preconditioner Setup");

			PW<KSP> solver;
			// setup petsc if needed
			if (solver_type == "petsc") {
				timer.start("Petsc Setup");

				KSPCreate(MPI_COMM_WORLD, &solver);
				KSPSetFromOptions(solver);
				KSPSetOperators(solver, A_petsc, A_petsc);
				if (M != nullptr) {
					PC M_petsc;
					KSPGetPC(solver, &M_petsc);
					PetscShellCreator::getPCShell(M_petsc, M, dc);
				}
				KSPSetUp(solver);
				KSPSetTolerances(solver, tolerance, PETSC_DEFAULT, PETSC_DEFAULT, 5000);

				timer.stop("Petsc Setup");
			}
			///////////////////
			// setup end
			///////////////////
			timer.stop("Linear System Setup");

			timer.start("Linear Solve");
			if (solver_type == "petsc") {
				KSPSolve(solver, f->vec, u->vec);
				int its;
				KSPGetIterationNumber(solver, &its);
				if (my_global_rank == 0) {
					cout << "Iterations: " << its << endl;
				}
			} else {
				std::shared_ptr<VectorGenerator<3>> vg(new ValVectorGenerator<3>(dc));

				int its = BiCGStab<3>::solve(vg, A, u, f, M);
				if (my_global_rank == 0) {
					cout << "Iterations: " << its << endl;
				}
			}
			timer.stop("Linear Solve");

			A->apply(u, au);
		}

		shared_ptr<PetscVector<3>> resid = PetscVector<3>::GetNewVector(dc);
		resid->addScaled(-1, au, 1, f);

		double residual = resid->twoNorm();
		double fnorm    = f->twoNorm();

		// error
		shared_ptr<PetscVector<3>> error = PetscVector<3>::GetNewVector(dc);
		error->addScaled(-1, exact, 1, u);
		if (neumann) {
			double uavg = dc->integrate(u) / dc->volume();
			double eavg = dc->integrate(exact) / dc->volume();

			if (my_global_rank == 0) {
				cout << "Average of computed solution: " << uavg << endl;
				cout << "Average of exact solution: " << eavg << endl;
			}

			error->shift(eavg - uavg);
		}
		double error_norm     = error->twoNorm();
		double error_norm_inf = error->infNorm();
		double exact_norm     = exact->twoNorm();

		double ausum = dc->integrate(au);
		double fsum  = dc->integrate(f);
		if (my_global_rank == 0) {
			std::cout << std::scientific;
			std::cout.precision(13);
			std::cout << "Error (2-norm):   " << error_norm / exact_norm << endl;
			std::cout << "Error (inf-norm): " << error_norm_inf << endl;
			std::cout << "Residual: " << residual / fnorm << endl;
			std::cout << u8"ΣAu-Σf: " << ausum - fsum << endl;
			cout.unsetf(std::ios_base::floatfield);
			int total_cells = dc->getNumGlobalCells();
			cout << "Total cells: " << total_cells << endl;
			cout << "Cores: " << num_procs << endl;
		}

		// output
		if (gamma_filename != "") {
			PetscViewer viewer;
			PetscViewerASCIIOpen(PETSC_COMM_WORLD, gamma_filename.c_str(), &viewer);
			VecView(gamma->vec, viewer);
		}
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
	}

#ifdef ENABLE_AMGX
	if (amgxsolver != nullptr) {
		delete amgxsolver;
	}
#endif
	if (my_global_rank == 0) {
		cout << timer;
	}
	PetscFinalize();
	return 0;
}
