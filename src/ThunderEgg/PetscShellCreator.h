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

#ifndef THUNDEREGG_PETSCSHELLCREATOR_H
#define THUNDEREGG_PETSCSHELLCREATOR_H
#include <ThunderEgg/Operator.h>
#include <petscpc.h>
namespace ThunderEgg
{
/**
 * @brief Class that contains methods for wrapping ThunderEgg operators for use in petsc
 */
class PetscShellCreator
{
	private:
	template <size_t D> class PetscPCShellOpDomain
	{
		private:
		std::shared_ptr<Operator<D>> op;
		std::shared_ptr<Domain<D>>   domain;

		public:
		PetscPCShellOpDomain(std::shared_ptr<Operator<D>> op, std::shared_ptr<Domain<D>> domain)
		{
			this->op = op;

			this->domain = domain;
		}
		static int applyMat(Mat A, Vec x, Vec b)
		{
			PetscPCShellOpDomain<D> *wrap = nullptr;
			MatShellGetContext(A, &wrap);
			std::shared_ptr<const Vector<D>> x_vec(
			new PetscVector<D>(x, -1, wrap->domain->getNs(), false));
			std::shared_ptr<Vector<D>> b_vec(
			new PetscVector<D>(b, -1, wrap->domain->getNs(), false));
			wrap->op->apply(x_vec, b_vec);
			return 0;
		}
		static int destroyMat(Mat A)
		{
			PetscPCShellOpDomain<D> *wrap = nullptr;
			MatShellGetContext(A, &wrap);
			delete wrap;
			return 0;
		}
		static int apply(PC P, Vec x, Vec b)
		{
			PetscPCShellOpDomain<D> *wrap = nullptr;
			PCShellGetContext(P, (void **) &wrap);
			std::shared_ptr<Vector<D>> x_vec(
			new PetscVector<D>(x, -1, wrap->domain->getNs(), false));
			std::shared_ptr<Vector<D>> b_vec(
			new PetscVector<D>(b, -1, wrap->domain->getNs(), false));
			wrap->op->apply(x_vec, b_vec);
			return 0;
		}
		static int destroy(PC P)
		{
			PetscPCShellOpDomain<D> *wrap = nullptr;
			PCShellGetContext(P, (void **) &wrap);
			delete wrap;
			return 0;
		}
	};
	template <size_t D> class PetscPCShellOpSchur
	{
		private:
		std::shared_ptr<Operator<D>>               op;
		std::shared_ptr<Schur::SchurHelper<D + 1>> sh;

		public:
		PetscPCShellOpSchur(std::shared_ptr<Operator<D>>               op,
		                    std::shared_ptr<Schur::SchurHelper<D + 1>> sh)
		{
			this->op = op;
			this->sh = sh;
		}
		static int applyMat(Mat A, Vec x, Vec b)
		{
			PetscPCShellOpSchur<D> *wrap = nullptr;
			MatShellGetContext(A, &wrap);
			std::shared_ptr<Vector<D>> x_vec(
			new PetscVector<D>(x, -1, wrap->sh->getLengths(), false));
			std::shared_ptr<Vector<D>> b_vec(
			new PetscVector<D>(b, -1, wrap->sh->getLengths(), false));
			wrap->op->apply(x_vec, b_vec);
			return 0;
		}
		static int destroyMat(Mat A)
		{
			PetscPCShellOpSchur<D> *wrap = nullptr;
			MatShellGetContext(A, &wrap);
			delete wrap;
			return 0;
		}
		static int apply(PC P, Vec x, Vec b)
		{
			PetscPCShellOpSchur<D> *wrap = nullptr;
			PCShellGetContext(P, (void **) &wrap);
			std::shared_ptr<Vector<D>> x_vec(
			new PetscVector<D>(x, -1, wrap->sh->getLengths(), false));
			std::shared_ptr<Vector<D>> b_vec(
			new PetscVector<D>(b, -1, wrap->sh->getLengths(), false));
			wrap->op->apply(x_vec, b_vec);
			return 0;
		}
		static int destroy(PC P)
		{
			PetscPCShellOpSchur<D> *wrap = nullptr;
			PCShellGetContext(P, (void **) &wrap);
			delete wrap;
			return 0;
		}
	};

	public:
	/**
	 * @brief Wrap a ThunderEgg operator for use as a preconditioner in petsc
	 *
	 * @tparam D the number of cartesian dimensions on a patch
	 * @param pc the PC object from petsc
	 * @param op the operator that we are wrapping
	 * @param domain the domain of the problem
	 */
	template <size_t D>
	static void getPCShell(PC pc, std::shared_ptr<Operator<D>> op,
	                       std::shared_ptr<Domain<D>> domain)
	{
		PetscPCShellOpDomain<D> *wrap = new PetscPCShellOpDomain<D>(op, domain);
		PCSetType(pc, PCSHELL);
		PCShellSetContext(pc, wrap);
		PCShellSetApply(pc, PetscPCShellOpDomain<D>::apply);
		PCShellSetDestroy(pc, PetscPCShellOpDomain<D>::destroy);
	}
	/**
	 * @brief Wrap a ThunderEgg operator for use as a preconditioner in petsc
	 *
	 * This is for wrapping preconditioners of the Schur matrix.
	 *
	 * @tparam D the number of cartesian dimensions on a patch
	 * @param pc the PC object from petsc
	 * @param op the operator that we are wrapping
	 * @param sh the SchurHelper of the problem
	 */
	template <size_t D>
	static void getPCShell(PC pc, std::shared_ptr<Operator<D>> op,
	                       std::shared_ptr<Schur::SchurHelper<D + 1>> sh)
	{
		PetscPCShellOpSchur<D> *wrap = new PetscPCShellOpSchur<D>(op, sh);
		PCSetType(pc, PCSHELL);
		PCShellSetContext(pc, wrap);
		PCShellSetApply(pc, PetscPCShellOpSchur<D>::apply);
		PCShellSetDestroy(pc, PetscPCShellOpSchur<D>::destroy);
	}
	/**
	 * @brief Wrap a ThunderEgg operator for use in petsc solver
	 *
	 * @tparam D the number of cartesian dimensions in a patch
	 * @param op the operator that we are wrapping
	 * @param domain the domain of the problem
	 * @return PW_explicit<Mat> Return a new Mat object
	 */
	template <size_t D>
	static PW_explicit<Mat> getMatShell(std::shared_ptr<Operator<D>> op,
	                                    std::shared_ptr<Domain<D>>   domain)
	{
		PetscPCShellOpDomain<D> *wrap = new PetscPCShellOpDomain<D>(op, domain);
		int                      M    = domain->getNumGlobalCells();
		int                      m    = domain->getNumLocalCells();
		PW<Mat>                  A;
		MatCreateShell(MPI_COMM_WORLD, m, m, M, M, wrap, &A);
		MatShellSetOperation(A, MATOP_MULT, (void (*)(void)) PetscPCShellOpDomain<D>::applyMat);
		MatShellSetOperation(A, MATOP_DESTROY,
		                     (void (*)(void)) PetscPCShellOpDomain<D>::destroyMat);
		return A;
	}
	/**
	 * @brief Wrap a ThunderEgg operator for use in petsc solver
	 *
	 * This is for use in schur matrix operators
	 *
	 * @tparam D the number of cartesian dimensions in a patch
	 * @param op the operator that we are wrapping
	 * @param sh the SchurHelper of the problem
	 * @return PW_explicit<Mat> Return a new Mat object
	 */
	template <size_t D>
	static PW_explicit<Mat> getMatShell(std::shared_ptr<Operator<D>>           op,
	                                    std::shared_ptr<Schur::SchurHelper<D>> sh)
	{
		PetscPCShellOpSchur<D> *wrap = new PetscPCShellOpSchur<D>(op, sh);
		int                     M    = sh->getSchurVecGlobalSize();
		int                     m    = sh->getSchurVecLocalSize();
		PW<Mat>                 A;
		MatCreateShell(MPI_COMM_WORLD, m, m, M, M, wrap, &A);
		MatShellSetOperation(A, MATOP_MULT, (void (*)(void)) PetscPCShellOpSchur<D>::applyMat);
		MatShellSetOperation(A, MATOP_DESTROY, (void (*)(void)) PetscPCShellOpSchur<D>::destroyMat);
		return A;
	}
};
} // namespace ThunderEgg
#endif
