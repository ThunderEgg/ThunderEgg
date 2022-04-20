/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured
 *  Cartesian grids.
 *
 *  Copyright (c) 2020-2021 Scott Aiton
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

#ifndef THUNDEREGG_GMG_MPIRESTRICTOR_H
#define THUNDEREGG_GMG_MPIRESTRICTOR_H
/**
 * @file
 *
 * @brief MPIRestrictor class
 */
#include <ThunderEgg/Domain.h>
#include <ThunderEgg/GMG/InterLevelComm.h>
#include <ThunderEgg/GMG/Level.h>
#include <ThunderEgg/GMG/Restrictor.h>

namespace ThunderEgg::GMG {
/**
 * @brief Base class that makes the necessary mpi calls, derived classes only have to
 * implement restrictPatches() method
 */
template<int D>
  requires is_supported_dimension<D>
class MPIRestrictor : public Restrictor<D>
{
private:
  /**
   * @brief Implimentation class
   */
  class Implimentation;

  /**
   * @brief pointer to the implimentation
   */
  std::unique_ptr<const Implimentation> implimentation;

public:
  /**
   * @brief Create new MPIRestrictor object.
   *
   * @param ilc the communcation package for the two levels.
   */
  MPIRestrictor(const Domain<D>& coarser_domain, const Domain<D>& finer_domain);

  /**
   * @brief Destroy the MPIRestrictor object
   */
  ~MPIRestrictor();

  /**
   * @brief Copy construct a new Inter Level Comm object
   *
   * @param other the MPIRestrictor to copy
   */
  MPIRestrictor(const MPIRestrictor& other);

  /**
   * @brief Move construct a new Inter Level Comm object
   *
   * @param other the MPIRestrictor to copy
   */
  MPIRestrictor(MPIRestrictor&& other);

  /**
   * @brief Copy assign a new Inter Level Comm object
   *
   * @param other the MPIRestrictor to copy
   */
  MPIRestrictor&
  operator=(const MPIRestrictor& other);

  /**
   * @brief Move assign a new MPIRestrictor object
   *
   * @param other the InterLevelComm to copy
   */
  MPIRestrictor&
  operator=(MPIRestrictor&& other);

  Vector<D> restrict(const Vector<D>& fine) const override;

  /**
   * @brief Restrict values into coarse vector
   *
   * The idea behind this is that this function will be called twice. Once to fill in the ghost
   * values, and once to fill in the local values. The ghost values will be filled first and the
   * local values will be fill while MPI communication is happening.
   *
   * @param patches pairs where the first value is the index in the coarse vector and the second
   * value is a reference to the PatchInfo object
   * @param finer_vector the finer vector
   * @param coarser_vector the coarser vector
   */
  virtual void
  restrictPatches(
    const std::vector<std::pair<int, std::reference_wrapper<const PatchInfo<D>>>>& patches,
    const Vector<D>& finer_vector,
    Vector<D>& coarser_vector) const = 0;
};
} // namespace ThunderEgg::GMG
// explicit instantiation
extern template class ThunderEgg::GMG::MPIRestrictor<2>;
extern template class ThunderEgg::GMG::MPIRestrictor<3>;
#endif