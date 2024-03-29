/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured
 *  Cartesian grids.
 *
 *  Copyright (c) 2019-2021 Scott Aiton
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

#ifndef THUNDEREGG_POISSON_FASTSCHURMATRIXASSEMBLE3D_H
#define THUNDEREGG_POISSON_FASTSCHURMATRIXASSEMBLE3D_H
/**
 * @file
 *
 * @brief FastSchurMatrixAssemble3D class
 */
#include <ThunderEgg/Poisson/FFTWPatchSolver.h>
#include <ThunderEgg/Schur/InterfaceDomain.h>
#include <petscmat.h>
namespace ThunderEgg::Poisson {
/**
 * @brief A fast algorithm for forming the Schur compliment matrix
 *
 * Currently this algorithm only supports the FFTWPatchSolver and it has to use
 * TriLinearGhostFiller
 *
 * @param iface_domain the interface domain to form the schur compliment matrix for
 * @param solver the patch solver to use for the formation
 * @return Mat the PETSc matrix, user is responsible for destroying
 */
Mat
FastSchurMatrixAssemble3D(const Schur::InterfaceDomain<3>& iface_domain,
                          Poisson::FFTWPatchSolver<3>& solver);
} // namespace ThunderEgg::Poisson
#endif
