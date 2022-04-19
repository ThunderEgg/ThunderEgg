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

#ifndef THUNDEREGG_GMG_MPIINTERPOLATOR_H
#define THUNDEREGG_GMG_MPIINTERPOLATOR_H
/**
 * @file
 *
 * @brief MPIInterpolator class
 */
#include <ThunderEgg/Domain.h>
#include <ThunderEgg/GMG/InterLevelComm.h>
#include <ThunderEgg/GMG/Interpolator.h>
#include <ThunderEgg/GMG/Level.h>

namespace ThunderEgg::GMG {
/**
 * @brief Base class that makes the necessary mpi calls, derived classes only have to
 * implement interpolatePatches() method
 */
template<int D>
  requires is_supported_dimension<D>
class MPIInterpolator : public Interpolator<D>
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
   * @brief Create new MPIInterpolator object.
   *
   * @param ilc the communcation package for the two levels.
   */
  MPIInterpolator(const Domain<D>& coarser_domain, const Domain<D>& finer_domain);

  /**
   * @brief Destroy the MPIInterpolator object
   */
  ~MPIInterpolator();

  /**
   * @brief Copy constructor
   *
   * @param other other MPIInterpolator
   */
  MPIInterpolator(const MPIInterpolator& other);

  /**
   * @brief Move constructor
   *
   * @param other MPIInterpolator
   */
  MPIInterpolator(MPIInterpolator&& other);

  /**
   * @brief copy assignment
   *
   * @param other MPIInterpolator
   * @return MPIInterpolator& this
   */
  MPIInterpolator&
  operator=(const MPIInterpolator& other);

  /**
   * @brief move assignment
   *
   * @param other MPIInterpolator
   * @return MPIInterpolator& this
   */
  MPIInterpolator&
  operator=(MPIInterpolator&& other);

  /**
   * @brief Interpolate values from coarse vector to the finer vector
   *
   * The idea behind this is that this function will be called twice. Once to interpolate from the
   * local values, and once to interpolate from the ghost values. The local values will be
   * interpolated from first, while MPI communication is happening, and the ghost values will be
   * interpolated from last.
   *
   * @param patches pairs where the first value is the index in the coarse vector and the second
   * value is a reference to the PatchInfo object
   * @param finer_vector the finer vector
   * @param coarser_vector the coarser vector
   */
  virtual void
  interpolatePatches(
    const std::vector<std::pair<int, std::reference_wrapper<const PatchInfo<D>>>>& patches,
    const Vector<D>& coarser_vector,
    Vector<D>& finer_vector) const = 0;

  /**
   * @brief interpolation function
   *
   * @param coarse the input vector that is interpolated from
   * @param fine the output vector that is interpolated to.
   */
  void
  interpolate(const Vector<D>& coarse, Vector<D>& fine) const;
};

// explicit instantiation

extern template class MPIInterpolator<2>;
extern template class MPIInterpolator<3>;

} // namespace ThunderEgg::GMG
#endif