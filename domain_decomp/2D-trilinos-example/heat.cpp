#include "DomainCollection.h"
#include "Factory.h"
#include "SchurHelper.h"
//#include "FunctionWrapper.h"
#include "FivePtPatchOperator.h"
#include "Init.h"
#include "MyTypeDefs.h"
#include "PatchSolvers/FftwPatchSolver.h"
#include "PatchSolvers/FishpackPatchSolver.h"
#include "QuadInterpolator.h"
#include "Writers/ClawWriter.h"
#include "Writers/MMWriter.h"
#ifdef ENABLE_AMGX
#include "AmgxWrapper.h"
#endif
#ifdef ENABLE_HYPRE
#include "HypreWrapper.h"
#endif
#ifdef HAVE_VTK
#include "Writers/VtkWriter.h"
#endif
#include "Timer.h"
#include "args.h"
#include <Amesos2.hpp>
#include <Amesos2_Version.hpp>
#include <BelosBiCGStabSolMgr.hpp>
#include <BelosBlockCGSolMgr.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosConfigDefs.hpp>
#include <BelosGCRODRSolMgr.hpp>
#include <BelosLSQRSolMgr.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosTpetraAdapter.hpp>
#include <Ifpack2_BlockRelaxation_decl.hpp>
#include <Ifpack2_BlockRelaxation_def.hpp>
#include <Ifpack2_Factory.hpp>
#include <Ifpack2_ILUT_decl.hpp>
#include <Ifpack2_ILUT_def.hpp>
#include <Ifpack2_Relaxation_decl.hpp>
#include <Ifpack2_Relaxation_def.hpp>
#include <MatrixMarket_Tpetra.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Tpetra_Experimental_BlockCrsMatrix_Helpers.hpp>
#include <Xpetra_TpetraBlockCrsMatrix.hpp>
#include <chrono>
#include <cmath>
#include <iostream>
//#include <mpi.h>
#include <string>
#include <unistd.h>
#ifndef M_PIl
#define M_PIl 3.141592653589793238462643383279502884L /* pi */
#endif

using Teuchos::RCP;
using Teuchos::rcp;

// =========== //
// main driver //
// =========== //

using namespace std;

int main(int argc, char *argv[])
{
	//    using Belos::FuncWrap;
	using Teuchos::RCP;
	using Teuchos::rcp;
	using namespace std::chrono;

	Teuchos::GlobalMPISession     global(&argc, &argv, nullptr);
	RCP<const Teuchos::Comm<int>> comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();

	int num_procs = comm->getSize();

	int my_global_rank = comm->getRank();

	// parse input
	args::ArgumentParser parser("");

	// program options
	args::HelpFlag       help(parser, "help", "Display this help menu", {'h', "help"});
	args::ValueFlag<int> f_n(parser, "n", "number of cells in the x direction, in each domain",
	                         {'n'});
	args::ValueFlag<int> f_l(parser, "n", "run the program n times and print out the average",
	                         {'l'});
	args::Flag           f_wrapper(parser, "wrapper", "use a function wrapper", {"wrap"});
	args::ValueFlag<int> f_div(parser, "divide", "use iterative method", {"divide"});
	args::Flag f_nozerof(parser, "zerou", "don't  make the rhs match the boundary conditions",
	                     {"nozerof"});
	args::ValueFlag<double> f_t(
	parser, "tolerance", "set the tolerance of the iterative solver (default is 1e-10)", {'t'});

	// mesh options
	args::ValueFlag<string> f_mesh(parser, "file_name", "read in a mesh", {"mesh"});
	args::ValueFlag<int>    f_square(parser, "num_domains",
	                              "create a num_domains x num_domains square of grids", {"square"});
	args::ValueFlag<int> f_amr(parser, "num_domains", "create a num_domains x num_domains square "
	                                                  "of grids, and a num_domains*2 x "
	                                                  "num_domains*2 refined square next to it",
	                           {"amr"});

	// output options
	args::Flag f_outclaw(parser, "outclaw", "output amrclaw ascii file", {"outclaw"});
#ifdef HAVE_VTK
	args::ValueFlag<string> f_outvtk(parser, "", "output to vtk format", {"outvtk"});
#endif
	args::ValueFlag<string> f_m(parser, "matrix filename", "the file to write the matrix to",
	                            {'m'});
	args::ValueFlag<string> f_s(parser, "solution filename", "the file to write the solution to",
	                            {'s'});
	args::ValueFlag<string> f_resid(parser, "residual filename",
	                                "the file to write the residual to", {"residual"});
	args::ValueFlag<string> f_error(parser, "error filename", "the file to write the error to",
	                                {"error"});
	args::ValueFlag<string> f_r(parser, "rhs filename", "the file to write the rhs vector to",
	                            {'r'});
	args::ValueFlag<string> f_g(parser, "gamma filename", "the file to write the gamma vector to",
	                            {'g'});
	args::ValueFlag<string> f_read_gamma(parser, "gamma filename",
	                                     "the file to read gamma vector from", {"readgamma"});

	// problem options
	args::Flag f_gauss(parser, "gauss", "solve gaussian function", {"gauss"});
	args::Flag f_zero(parser, "zero", "solve zero function", {"zero"});

	// preconditioners
	args::Flag f_precj(parser, "prec", "use jacobi preconditioner", {"precj"});
	args::Flag f_prec(parser, "prec", "use block diagonal jacobi preconditioner", {"prec"});
	args::Flag f_precmuelu(parser, "prec", "use MueLu AMG preconditioner", {"muelu"});
	args::Flag f_neumann(parser, "neumann", "use neumann boundary conditions", {"neumann"});
	args::Flag f_cg(parser, "gmres", "use CG for iterative solver", {"cg"});
	args::Flag f_gmres(parser, "gmres", "use GMRES for iterative solver", {"gmres"});
	args::Flag f_lsqr(parser, "gmres", "use least squares for iterative solver", {"lsqr"});
	args::Flag f_rgmres(parser, "rgmres", "use GCRO-DR (Recycling GMRES) for iterative solver",
	                    {"rgmres"});
	args::Flag f_bicg(parser, "gmres", "use BiCGStab for iterative solver", {"bicg"});

	// direct solvers
	args::Flag f_iter(parser, "iterative", "use iterative method", {"iterative"});

	// patch solvers
	args::Flag f_fish(parser, "fishpack", "use fishpack as the patch solver", {"fishpack"});
#ifdef __NVCC__
	args::Flag f_cufft(parser, "cufft", "use CuFFT as the patch solver", {"cufft"});
#endif

	// third-party preconditioners
	args::Flag f_amgx(parser, "amgx", "solve schur compliment system with amgx", {"amgx"});
	args::Flag f_hypre(parser, "hypre", "solve schur compliment system with hypre", {"hypre"});

	if (argc < 5) {
		if (my_global_rank == 0) std::cout << parser;
		return 0;
	}
	try {
		parser.ParseCLI(argc, argv);
	} catch (args::Help) {
		if (my_global_rank == 0) std::cout << parser;
		return 0;
	} catch (args::ParseError e) {
		if (my_global_rank == 0) {
			std::cerr << e.what() << std::endl;
			std::cerr << parser;
		}
		return 1;
	} catch (args::ValidationError e) {
		if (my_global_rank == 0) {
			std::cerr << e.what() << std::endl;
			std::cerr << parser;
		}
		return 1;
	}


	DomainCollection dc;
	if (f_mesh) {
		string d = args::get(f_mesh);
		dc       = DomainCollection(comm, d, comm->getRank());
	} else if (f_amr) {
		int d = args::get(f_amr);
		dc    = DomainCollection(comm, d, d, comm->getRank(), true);
	} else {
		int d = args::get(f_square);
		dc    = DomainCollection(comm, d, d, comm->getRank());
	}
	if (f_div) {
		for (int i = 0; i < args::get(f_div); i++) {
			dc.divide();
		}
	}
	if (f_neumann) {
		dc.setNeumann();
	}
	// Set the number of discretization points in the x and y direction.
	int nx = args::get(f_n);
	int ny = args::get(f_n);
	dc.n   = nx;
	for (auto &p : dc.domains) {
		p.second.n = nx;
	}
	int total_cells = dc.num_global_domains * nx * ny;
	cerr << "Total cells: " << total_cells << endl;

	if (dc.num_global_domains < num_procs) {
		std::cerr << "number of domains must be greater than or equal to the number of processes\n";
		return 1;
	}
	// partition domains
	if (num_procs > 1) {
		dc.zoltanBalance();
	}

	scalar_type tol = 1e-12;
	if (f_t) {
		tol = args::get(f_t);
	}

	int loop_count = 1;
	if (f_l) {
		loop_count = args::get(f_l);
	}

	string save_matrix_file = "";
	if (f_m) {
		save_matrix_file = args::get(f_m);
	}

	string save_solution_file = "";
	if (f_s) {
		save_solution_file = args::get(f_s);
	}

	string save_residual_file = "";
	if (f_resid) {
		save_residual_file = args::get(f_resid);
	}
	string save_error_file = "";
	if (f_error) {
		save_error_file = args::get(f_error);
	}
	string save_rhs_file = "";
	if (f_r) {
		save_rhs_file = args::get(f_r);
	}
	string save_gamma_file = "";
	if (f_g) {
		save_gamma_file = args::get(f_g);
	}

	// the functions that we are using
	function<double(double, double)> ffun;
	function<double(double, double)> gfun;
	function<double(double, double)> nfunx;
	function<double(double, double)> nfuny;

	if (f_zero) {
		ffun  = [](double x, double y) { return 0; };
		gfun  = [](double x, double y) { return 0; };
		nfunx = [](double x, double y) { return 0; };
		nfuny = [](double x, double y) { return 0; };
	} else if (f_gauss) {
		gfun
		= [](double x, double y) { return exp(cos(10 * M_PIl * x)) - exp(cos(11 * M_PIl * y)); };
		ffun = [](double x, double y) {
			return 100 * M_PIl * M_PIl * (pow(sin(10 * M_PIl * x), 2) - cos(10 * M_PIl * x))
			       * exp(cos(10 * M_PIl * x))
			       + 121 * M_PIl * M_PIl * (cos(11 * M_PIl * y) - pow(sin(11 * M_PIl * y), 2))
			         * exp(cos(11 * M_PIl * y));
		};
		nfunx = [](double x, double y) {
			return -10 * M_PIl * sin(10 * M_PIl * x) * exp(cos(10 * M_PIl * x));
		};

		nfuny = [](double x, double y) {
			return 11 * M_PIl * sin(11 * M_PIl * y) * exp(cos(11 * M_PIl * y));
		};
	} else {
		ffun = [](double x, double y) {
			return -5 * M_PIl * M_PIl * sinl(M_PIl * y) * cosl(2 * M_PIl * x);
		};
		gfun = [](double x, double y) { return sinl(M_PIl * y) * cosl(2 * M_PIl * x); };
		nfunx
		= [](double x, double y) { return -2 * M_PIl * sinl(M_PIl * y) * sinl(2 * M_PIl * x); };
		nfuny = [](double x, double y) { return M_PIl * cosl(M_PIl * y) * cosl(2 * M_PIl * x); };
	}

	// set the patch solver
	RCP<PatchSolver> p_solver;
	if (f_fish) {
		p_solver = rcp(new FishpackPatchSolver());
	} else {
		p_solver = rcp(new FftwPatchSolver(dc));
	}

	// patch operator
	RCP<PatchOperator> p_operator = rcp(new FivePtPatchOperator());

	// interface interpolator
	RCP<Interpolator> p_interp = rcp(new QuadInterpolator());

	Tools::Timer timer;
	for (int loop = 0; loop < loop_count; loop++) {
		timer.start("Domain Initialization");

		SchurHelper sch(dc, comm, p_solver, p_operator, p_interp);

		RCP<map_type>    domain_map = dc.getDomainRowMap();
		RCP<vector_type> u          = rcp(new vector_type(domain_map, 1));
		RCP<vector_type> exact      = rcp(new vector_type(domain_map, 1));
		RCP<vector_type> f          = rcp(new vector_type(domain_map, 1));

		if (f_neumann) {
			double *f_ptr     = &f->get1dViewNonConst()[0];
			double *exact_ptr = &exact->get1dViewNonConst()[0];
			Init::initNeumann(dc, nx, f_ptr, exact_ptr, ffun, gfun, nfunx, nfuny);
		} else {
			double *f_ptr     = &f->get1dViewNonConst()[0];
			double *exact_ptr = &exact->get1dViewNonConst()[0];
			Init::initDirichlet(dc, nx, f_ptr, exact_ptr, ffun, gfun);
		}

		timer.stop("Domain Initialization");

		// Create a map that will be used in the iterative solver
		RCP<map_type>       matrix_map       = dc.getSchurRowMap();
		RCP<const map_type> matrix_map_const = matrix_map;

		// Create the gamma and diff vectors
		RCP<vector_type>                   gamma = rcp(new vector_type(matrix_map, 1));
		RCP<vector_type>                   diff  = rcp(new vector_type(matrix_map, 1));
		RCP<vector_type>                   b     = rcp(new vector_type(matrix_map, 1));
		RCP<matrix_type>                   A;
		RCP<Tpetra::Operator<scalar_type>> op;
		RCP<const Tpetra::RowMatrix<>>     rm;
		RCP<Amesos2::Solver<matrix_type, vector_type>> dsolver;

		// Create linear problem for the Belos solver
		RCP<Belos::LinearProblem<scalar_type, vector_type, Tpetra::Operator<scalar_type>>> problem;
		RCP<Belos::SolverManager<scalar_type, vector_type, Tpetra::Operator<scalar_type>>> solver;
		Teuchos::ParameterList belosList;

		typedef Ifpack2::Preconditioner<scalar_type> Preconditioner;
		RCP<Preconditioner>                          prec;
		if (f_neumann && !f_nozerof) {
			double fdiff = (dc.integrate(*f)) / dc.area();
			if (my_global_rank == 0) cout << "Fdiff: " << fdiff << endl;

			RCP<vector_type> diff = rcp(new vector_type(domain_map, 1));
			diff->putScalar(fdiff);
			f->update(1.0, *diff, 1.0);
		}

#ifdef ENABLE_AMGX
		Teuchos::RCP<AmgxWrapper> amgxsolver;
#endif
#ifdef ENABLE_HYPRE
		Teuchos::RCP<HypreWrapper> hypresolver;
#endif

		steady_clock::time_point tsolve_start;
		if (dc.num_global_domains != 1) {
			// do iterative solve

			// Get the b vector
			sch.solveWithInterface(*f, *u, *gamma, *b);
			b->scale(-1.0);

			if (save_rhs_file != "") {
				Tpetra::MatrixMarket::Writer<matrix_type>::writeDenseFile(save_rhs_file, b, "", "");
			}

			///////////////////
			// setup start
			///////////////////
			timer.start("Linear System Setup");

			if (f_wrapper) {
				// Create a function wrapper
				// op = rcp(new FuncWrap(b, &dc));
			} else {
				timer.start("Matrix Formation");

				A = sch.formCRSMatrix();

				timer.stop("Matrix Formation");

				op = A;
				rm = A;

				if (save_matrix_file != "")
					Tpetra::MatrixMarket::Writer<matrix_type>::writeSparseFile(save_matrix_file, A,
					                                                           "", "");
			}

			problem
			= rcp(new Belos::LinearProblem<scalar_type, vector_type, Tpetra::Operator<scalar_type>>(
			op, gamma, b));
			if (f_amgx) {
#ifdef ENABLE_AMGX
				timer.start("AMGX Setup");
				amgxsolver = rcp(new AmgxWrapper(A, dc, nx));
				timer.stop("AMGX Setup");
#endif
			} else if (f_hypre) {
#ifdef ENABLE_HYPRE
				timer.start("Hypre Setup");
				hypresolver = rcp(new HypreWrapper(A, dc, nx, tol, true));
				timer.stop("Hypre Setup");
#endif
			} else {
				timer.start("MueLu Preconditioner Formation");

				Teuchos::RCP<op_type> P = Factory::getAmgPreconditioner(A);

				problem->setLeftPrec(P);

				timer.stop("MueLu Preconditioner Formation");
			}
			///////////////////
			// setup end
			///////////////////
			timer.stop("Linear System Setup");
		}
		///////////////////
		// solve start
		///////////////////
		timer.start("Complete Solve");

		if (dc.num_global_domains != 1) {
			timer.start("Gamma Solve");
			if (f_amgx) {
// solve
#ifdef ENABLE_AMGX
				amgxsolver->solve(gamma, b);
#endif
			} else if (f_hypre) {
#ifdef ENABLE_HYPRE
				hypresolver->solve(gamma, b);
#endif
			} else if (f_read_gamma) {
				gamma = Tpetra::MatrixMarket::Reader<matrix_type>::readDenseFile(
				args::get(f_read_gamma), comm, matrix_map_const);
			} else {
				problem->setProblem();

				// Set the parameters
				belosList.set("Block Size", 1);
				belosList.set("Maximum Iterations", 5000);
				belosList.set("Convergence Tolerance", tol);
				belosList.set("Output Frequency", 1);
				int verbosity = Belos::Errors + /* Belos::StatusTestDetails +*/ Belos::Warnings
				                + Belos::TimingDetails + Belos::Debug + Belos::IterationDetails;
				belosList.set("Verbosity", verbosity);
				// belosList.set("Orthogonalization", "ICGS");
				//
				belosList.set("Rel RHS Err", 0.0);
				belosList.set("Rel Mat Err", 0.0);

				// Create solver and solve
				if (f_rgmres) {
					solver = rcp(new Belos::GCRODRSolMgr<scalar_type, vector_type,
					                                     Tpetra::Operator<scalar_type>>(
					problem, rcp(&belosList, false)));
				} else if (f_cg) {
					solver = rcp(new Belos::BlockCGSolMgr<scalar_type, vector_type,
					                                      Tpetra::Operator<scalar_type>>(
					problem, rcp(&belosList, false)));
				} else if (f_bicg) {
					solver = rcp(new Belos::BiCGStabSolMgr<scalar_type, vector_type,
					                                       Tpetra::Operator<scalar_type>>(
					problem, rcp(&belosList, false)));
				} else {
					solver = rcp(new Belos::BlockGmresSolMgr<scalar_type, vector_type,
					                                         Tpetra::Operator<scalar_type>>(
					problem, rcp(&belosList, false)));
				}
				solver->solve();
			}
			timer.stop("Gamma Solve");
		}

		// Do one last solve
		timer.start("Patch Solve");

		sch.solveWithInterface(*f, *u, *gamma, *diff);

		timer.stop("Patch Solve");

		/*
		double ausum2 = sch.integrateAU();
		double fsum2  = sch.integrateF();
		double bflux  = sch.integrateBoundaryFlux();
		if (my_global_rank == 0) {
		    std::cout << u8"Σf-Au: " << fsum2 - ausum2 << endl;
		    std::cout << u8"Σf: " << fsum2 << endl;
		    std::cout << u8"ΣAu: " << ausum2 << endl;
		    if (f_neumann) {
		        std::cout << u8"∮ du/dn: " << bflux << endl;
		        std::cout << u8"∮ du/dn - Σf: " << bflux - fsum2 << endl;
		        std::cout << u8"∮ du/dn - ΣAu: " << bflux - ausum2 << endl;
		    }
		}
		*/
		/*
		if (f_iter && !direct_solve) {
		    timer.start("Iterative Refinement Step");
		    sch.residual();
		    sch.swapResidSol();

		    if (dc.num_global_domains != 1) {
		        x->putScalar(0);
		        sch.solveWithInterface(*x, *r);
		        // op->apply(*gamma, *r);
		        // r->update(1.0, *b, -1.0);

		        solver->reset(Belos::ResetType::Problem);
		        if (f_wrapper) {
		            ((FuncWrap *) op.getRawPtr())->setB(r);
		        }
		        problem->setProblem(x, r);
		        solver->setProblem(problem);
		        solver->solve();
		    }
		    sch.solveWithInterface(*x, *d);
		    sch.sumResidIntoSol();
		    timer.stop("Iterative Refinement Step");
		}
		*/

		///////////////////
		// solve end
		///////////////////
		timer.stop("Complete Solve");

		// residual
		RCP<vector_type> resid = rcp(new vector_type(domain_map, 1));
		RCP<vector_type> au    = rcp(new vector_type(domain_map, 1));
		sch.applyWithInterface(*u, *gamma, *au);
		resid->update(-1.0, *f, 1.0, *au, 0.0);
		double residual = resid->getVector(0)->norm2();
		double fnorm    = f->getVector(0)->norm2();

		// error
		RCP<vector_type> error = rcp(new vector_type(domain_map, 1));
		error->update(-1.0, *exact, 1.0, *u, 0.0);
		if (f_neumann) {
			double uavg = dc.integrate(*u) / dc.area();
			double eavg = dc.integrate(*exact) / dc.area();

			if (my_global_rank == 0) {
				cout << "Average of computed solution: " << uavg << endl;
				cout << "Average of exact solution: " << eavg << endl;
			}

			vector_type ones(domain_map, 1);
			ones.putScalar(1);
			error->update(eavg - uavg, ones, 1.0);
		}
		double error_norm = error->getVector(0)->norm2();
		double exact_norm = exact->getVector(0)->norm2();

		double ausum = dc.integrate(*au);
		double fsum  = dc.integrate(*f);
		if (my_global_rank == 0) {
			std::cout << std::scientific;
			std::cout.precision(13);
			std::cout << "Error: " << error_norm / exact_norm << endl;
			std::cout << "Residual: " << residual / fnorm << endl;
			std::cout << u8"ΣAu-Σf: " << ausum - fsum << endl;
			cout.unsetf(std::ios_base::floatfield);
		}
		MMWriter mmwriter(dc, f_amr);
		if (save_solution_file != "") {
			mmwriter.write(*u, save_solution_file);
		}
		if (save_residual_file != "") {
			mmwriter.write(*resid, save_residual_file);
		}
		if (save_error_file != "") {
			mmwriter.write(*error, save_error_file);
		}
		if (save_gamma_file != "") {
			Tpetra::MatrixMarket::Writer<matrix_type>::writeDenseFile(save_gamma_file, gamma, "",
			                                                          "");
		}
		if (f_outclaw) {
			ClawWriter writer(dc);
			writer.write(*u, *resid);
		}
#ifdef HAVE_VTK
		if (f_outvtk) {
			VtkWriter writer(dc);
			writer.write(args::get(f_outvtk),*u, *error, *resid);
		}
#endif
		cout.unsetf(std::ios_base::floatfield);
	}

	if (my_global_rank == 0) {
		cout << timer;
	}
	return 0;
}
