/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured
 *  Cartesian grids.
 *
 *  Copyright (c) 2021      Scott Aiton
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

#ifndef THUNDEREGG_PETSC_KSPSOLVER_H
#define THUNDEREGG_PETSC_KSPSOLVER_H
/**
 * @file
 *
 * @brief KSPSolver class
 */

#include <ThunderEgg/Iterative/Solver.h>
#include <ThunderEgg/PETSc/MatShellCreator.h>
#include <ThunderEgg/PETSc/PCShellCreator.h>
#include <petscksp.h>
#include <petscmat.h>

namespace ThunderEgg::PETSc {
/**
 * @brief Wraps a PETSc KSP solver for use use as a Solver
 */
template<int D>
class KSPSolver : public Iterative::Solver<D>
{
private:
  /**
   * @brief The PETSc matrix
   */
  KSPType type = KSPGMRES;
  /**
   * @brief Allocate a PETSc Vec without the additional padding for ghost cells
   *
   * @param vec the ThunderEgg Vector
   * @return Vec the PETSc Vec
   */
  Vec getPetscVecWithoutGhost(const Vector<D>& vec) const
  {
    Vec petsc_vec;
    VecCreateMPI(vec.getCommunicator().getMPIComm(),
                 vec.getNumLocalCells() * vec.getNumComponents(),
                 PETSC_DETERMINE,
                 &petsc_vec);
    return petsc_vec;
  }
  /**
   * @brief Copy to a PETSc Vec from a ThunderEgg Vector
   *
   * @param vec the ThunderEgg Vector
   * @param petsc_vec the PETSc Vec to copy to
   */
  void copyToPetscVec(const Vector<D>& vec, Vec petsc_vec) const
  {
    double* petsc_vec_view;
    size_t curr_index = 0;
    VecGetArray(petsc_vec, &petsc_vec_view);
    for (int i = 0; i < vec.getNumLocalPatches(); i++) {
      for (int c = 0; c < vec.getNumComponents(); c++) {
        const ComponentView<const double, D> ld = vec.getComponentView(c, i);
        Loop::Nested<D>(ld.getStart(), ld.getEnd(), [&](const std::array<int, D>& coord) {
          petsc_vec_view[curr_index] = ld[coord];
          curr_index++;
        });
      }
    }
    VecRestoreArray(petsc_vec, &petsc_vec_view);
  }
  /**
   * @brief Copy from a PETSc Vec to a ThunderEgg Vector
   *
   * @param petsc_vec the PETSc Vec
   * @param vec the ThunderEgg Vector
   */
  void copyToVec(Vec petsc_vec, Vector<D>& vec) const
  {
    const double* petsc_vec_view;
    size_t curr_index = 0;
    VecGetArrayRead(petsc_vec, &petsc_vec_view);
    for (int i = 0; i < vec.getNumLocalPatches(); i++) {
      for (int c = 0; c < vec.getNumComponents(); c++) {
        ComponentView<double, D> ld = vec.getComponentView(c, i);
        Loop::Nested<D>(ld.getStart(), ld.getEnd(), [&](const std::array<int, D>& coord) {
          ld[coord] = petsc_vec_view[curr_index];
          curr_index++;
        });
      }
    }
    VecRestoreArrayRead(petsc_vec, &petsc_vec_view);
  }
  /**
   * @brief Deallocate the PETSc vector
   *
   * @param vec the corresponding ThunderEgg Vector
   * @param petsc_vec the PETSc Vec to deallocate
   */
  void destroyPetscVec(const Vector<D>& vec, Vec petsc_vec) const { VecDestroy(&petsc_vec); }
  static int MonitorResidual(KSP ksp, PetscInt n, PetscReal rnorm, std::ostream* os)
  {
    char buf[100];
    sprintf(buf, "%5d %16.8e\n", n, rnorm);
    (*os) << std::string(buf);
    return 0;
  }

public:
  /**
   * @brief Construct a new MatWrapper object
   *
   * This object will not deallocate the PETSc Mat, you are responsible for deallocating it.
   *
   * @param A_in the PETSc Mat to wrap
   */
  explicit KSPSolver() {}

  KSPSolver<D>* clone() const override { return new KSPSolver<D>(*this); }

  void setType(KSPType new_type) { type = new_type; }
  int solve(const Operator<D>& A,
            Vector<D>& x,
            const Vector<D>& b,
            const Operator<D>* Mr = nullptr,
            bool output = false,
            std::ostream& os = std::cout) const override
  {
    Mat A_PETSC = MatShellCreator<D>::GetNewMatShell(A, [&]() { return x.getZeroClone(); });

    PC Mr_PETSC = nullptr;
    if (Mr != nullptr) {
      Mr_PETSC = PCShellCreator<D>::GetNewPCShell(*Mr, A, [&]() { return x.getZeroClone(); });
    }
    Vec x_PETSC = getPetscVecWithoutGhost(x);
    copyToPetscVec(x, x_PETSC);

    Vec b_PETSC = getPetscVecWithoutGhost(x);
    copyToPetscVec(b, b_PETSC);

    KSP ksp;
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetType(ksp, type);
    KSPSetOperators(ksp, A_PETSC, A_PETSC);

    if (Mr != nullptr) {
      KSPSetPC(ksp, Mr_PETSC);
      // KSPSetPCSide(ksp, PC_RIGHT);
    }

    if (output) {
      KSPMonitorSet(
        ksp, (PetscErrorCode(*)(KSP, PetscInt, PetscReal, void*)) & MonitorResidual, &os, nullptr);
    }

    KSPSolve(ksp, b_PETSC, x_PETSC);

    int iterations;
    KSPGetIterationNumber(ksp, &iterations);

    const char* converged_reason;
    KSPGetConvergedReasonString(ksp, &converged_reason);
    os << converged_reason << std::endl;

    KSPDestroy(&ksp);
    VecDestroy(&x_PETSC);
    VecDestroy(&b_PETSC);
    PCDestroy(&Mr_PETSC);
    MatDestroy(&A_PETSC);
    return iterations;
  }
};
} // namespace ThunderEgg::PETSc
#endif