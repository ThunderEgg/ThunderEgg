/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured
 *  Cartesian grids.
 *
 *  Copyright (c) 2018-2021 Scott Aiton
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

#ifndef THUNDEREGG_GMG_INTERLEVELCOMM_H
#define THUNDEREGG_GMG_INTERLEVELCOMM_H
/**
 * @file
 *
 * @brief InterLevelComm class
 */

#include <ThunderEgg/RuntimeError.h>
#include <ThunderEgg/Vector.h>

namespace ThunderEgg::GMG {
/**
 * @brief Facilitates communication between a finer domain and a coarser domain.
 *
 * This class will determine the following things:
 *
 * Which patches in the finer domain have a parent patch in the coarser domain on the same rank?
 * 	- getPatchesWithLocalParent() will return a vector of these patches and the local indexes of
 *    their parent patches.
 *
 * Which patches in the finer domain have a parent patch in the coarser domain on a different rank?
 * 	- getPatchesWithGhostParent() will return a vector of these patches and the local indexes in
 * the ghost vector.
 * 	- getNewGhostVector() will allocate a new vector for these coarse ghost values.
 *
 *
 * Scatter functions are provided for scattering to and from the coarse vector to the coarse ghost
 * vector.
 */
template<int D>
  requires is_supported_dimension<D>
class InterLevelComm
{
private:
  /**
   * @brief Implimentation class
   */
  class Implimentation;

  /**
   * @brief pointer to the implimentation
   */
  std::unique_ptr<Implimentation> implimentation;

public:
  /**
   * @brief Create a new InterLevelComm object.
   *
   * @param coarse_domain the coarser DomainCollection.
   * @param fine_domain the finer DomainCollection.
   */
  InterLevelComm(const Domain<D>& coarser_domain, const Domain<D>& finer_domain);

  /**
   * @brief Destroy the Inter Level Comm object
   */
  ~InterLevelComm();

  /**
   * @brief Copy construct a new Inter Level Comm object
   *
   * @param other the InterLevelComm to copy
   */
  InterLevelComm(const InterLevelComm& other);

  /**
   * @brief Move construct a new Inter Level Comm object
   *
   * @param other the InterLevelComm to copy
   */
  InterLevelComm(InterLevelComm&& other);

  /**
   * @brief Copy assign a new Inter Level Comm object
   *
   * @param other the InterLevelComm to copy
   */
  InterLevelComm&
  operator=(const InterLevelComm& other);

  /**
   * @brief Move assign a new Inter Level Comm object
   *
   * @param other the InterLevelComm to copy
   */
  InterLevelComm&
  operator=(InterLevelComm&& other);

  /**
   * @brief Allocate a new vector for ghost patch values
   * @param num_components the number of components
   * @return the newly allocated vector.
   */
  Vector<D>
  getNewGhostVector(int num_components) const;

  /**
   * @brief Get the vector of finer patches that have a local parent
   *
   * The vector will consist of pair values:
   * 		- First value: the local index of the parent patch
   * 		- Second value: a reference to the child patch
   *
   * @return const std::vector<std::pair<int, std::reference_wrapper<const PatchInfo<D>>>>& the
   * vector
   */
  const std::vector<std::pair<int, std::reference_wrapper<const PatchInfo<D>>>>&
  getPatchesWithLocalParent() const;

  /**
   * @brief Get the vector of finer patches that have a ghost parent
   *
   * The vector will consist of pair values:
   * 		- First value: the local index in the ghost vector of the parent patch
   * 		- Second value: a reference to the child patch
   *
   * @return const std::vector<std::pair<int, std::reference_wrapper<const PatchInfo<D>>>>& the
   * vector
   */
  const std::vector<std::pair<int, std::reference_wrapper<const PatchInfo<D>>>>&
  getPatchesWithGhostParent() const;

  /**
   * @brief Start the communication for sending ghost values.
   *
   * This will send the values in the ghost vector and will add them to the values in the vector.
   * This is essentially a reverse scatter.
   *
   * This function is seperated into a Start and Finish function, this allows for other
   * computations to happen while the communicating.
   *
   * @param vector the vector
   * @param ghost_vector the associated ghost vector
   */
  void
  sendGhostPatchesStart(Vector<D>& vector, const Vector<D>& ghost_vector);

  /**
   * @brief Finish the communication for sending ghost values.
   *
   * This will send the values in the ghost vector and will add them to the values in the
   * vector. This is essentially a reverse scatter.
   *
   * This function is seperated into a Start and Finish function, this allows for other
   * computations to happen while the communicating.
   *
   * @param vector the vector
   * @param ghost_vector the associated ghost vector
   */
  void
  sendGhostPatchesFinish(Vector<D>& vector, const Vector<D>& ghost_vector);

  /**
   * @brief Start the communication for getting ghost values.
   *
   * This will send the values in the vector to the ghost vector, and will overwrite the
   * values in the ghost vector. This is essentially a forward scatter.
   *
   * This function is seperated into a Start and Finish function, this allows for other
   * computations to happen while the communicating.
   *
   * @param vector the vector
   * @param ghost_vector the associated ghost vector
   */
  void
  getGhostPatchesStart(const Vector<D>& vector, Vector<D>& ghost_vector);

  /**
   * @brief Finish the communication for getting ghost values.
   *
   * This will send the values in the vector to the ghost vector, and will overwrite the
   * values in the ghost vector. This is essentially a forward scatter.
   *
   * This function is seperated into a Start and Finish function, this allows for other
   * computations to happen while communicating.
   *
   * @param vector the vector
   * @param ghost_vector the associated ghost vector
   */
  void
  getGhostPatchesFinish(const Vector<D>& vector, Vector<D>& ghost_vector);

  /**
   * @brief Get the coarser Domain
   *
   * @return const Domain<D>& the coarser Domain
   */
  const Domain<D>&
  getCoarserDomain() const;

  /**
   * @brief Get the finer Domain
   *
   * @return const Domain<D>& the finer Domain
   */
  const Domain<D>&
  getFinerDomain() const;
};
extern template class InterLevelComm<2>;
extern template class InterLevelComm<3>;
} // namespace ThunderEgg::GMG
#endif
