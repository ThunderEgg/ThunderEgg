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

#ifndef THUNDEREGG_DOMAIN_H
#define THUNDEREGG_DOMAIN_H
/**
 * @file
 *
 * @brief Domain class
 */
#include <ThunderEgg/PatchInfo.h>
#include <ThunderEgg/Timer.h>
#include <map>
#include <set>
#include <vector>

namespace ThunderEgg {
/**
 * @brief Implimentation of Domain class
 *
 * @tparam D the number of Cartesian dimensions
 */
template<int D>
  requires is_supported_dimension<D>
class DomainImpl;

/**
 * @brief Uses a collection of PatchInfo objects to represent the domain of the problem.
 *
 * Each patch passed to the constructor needs to have the following complete information:
 * 	* Each patch within a domain have to have a unique id set.
 *  * Each patch needs to have the neighboring patch information filled out, with the id of neighbor
 * patches.
 *
 *  When passed to the constructor, the constructor will create an indexing for the
 * This class mainly manages a set of patches that makes up the domain. It is responsible for
 * setting up the indexing of the domains, which is used in the rest of the ThunderEgg library.
 *
 * @tparam D the number of Cartesian dimensions
 */
template<int D>
  requires is_supported_dimension<D>
class Domain
{
private:
  /**
   * @brief implimenation of Domain
   */
  std::shared_ptr<const DomainImpl<D>> implimentation;
  /**
   * @brief The id of the domain
   */
  int id = -1;
  /**
   * @brief The timer
   */
  mutable std::shared_ptr<Timer> timer;

  /**
   * @brief
   *
   * @param comm
   * @param id
   * @param ns
   * @param num_ghost_cells
   * @param pinfos
   * @return DomainImpl<D>*
   */
  static std::shared_ptr<const DomainImpl<D>>
  NewDomainImpl(const Communicator& comm,
                int id,
                const std::array<int, D>& ns,
                int num_ghost_cells,
                std::vector<PatchInfo<D>>&& pinfos);

public:
  /**
   * @brief Construct a new Domain object
   *
   * @tparam InputIterator the iterator for PatchInfo objects
   * @param id the id of the domain should be unique within a multigrid cycle
   * @param ns the number of cells in each direction
   * @param num_ghost_cells the number of ghost cells on each side of the patch
   * @param first_pinfo start iterator for PatchInfo objects
   * @param last_pinfo end iterator for PatchInfo objects
   */
  template<class InputIterator>
  Domain(Communicator comm,
         int id,
         std::array<int, D> ns,
         int num_ghost_cells,
         InputIterator first_pinfo,
         InputIterator last_pinfo)
    : implimentation(NewDomainImpl(comm,
                                   id,
                                   ns,
                                   num_ghost_cells,
                                   std::vector<PatchInfo<D>>(first_pinfo, last_pinfo)))
    , id(id)
  {
  }

  /**
   * @brief Get the Communicator object associated with this domain
   *
   * @return const Communicator& the Communicator
   */
  const Communicator&
  getCommunicator() const;

  /**
   * @brief Get a vector of PatchInfo pointers where index in the vector corresponds to the
   * patch's local index
   */
  const std::vector<PatchInfo<D>>&
  getPatchInfoVector() const;

  /**
   * @brief Get the number of cells in each direction
   *
   */
  const std::array<int, D>&
  getNs() const;

  /**
   * @brief Get the number of global patches
   */
  int
  getNumGlobalPatches() const;

  /**
   * @brief Get the number of local patches
   */
  int
  getNumLocalPatches() const;

  /**
   * @brief get the number of global cells
   */
  int
  getNumGlobalCells() const;

  /**
   * @brief Get get the number of local cells
   */
  int
  getNumLocalCells() const;

  /**
   * @brief Get get the number of local cells (including ghost cells)
   */
  int
  getNumLocalCellsWithGhost() const;

  /**
   * @brief Get the number of cells in a patch
   */
  int
  getNumCellsInPatch() const;

  /**
   * @brief get the number of ghost cell on each side of a patch
   */
  int
  getNumGhostCells() const;

  /**
   * @brief Get the volume of the domain.
   *
   * For 2D, this will be the area.
   */
  double
  volume() const;

  /**
   * @brief Set the Timer object
   *
   * @param timer the timer
   */
  void
  setTimer(std::shared_ptr<Timer> timer) const;

  /**
   * @brief Get the Timer object
   *
   * @return std::shared_ptr<Timer> the timer
   */
  std::shared_ptr<Timer>
  getTimer() const;

  /**
   * @brief Check if the Domain has a timer associated with it
   * @return true if the Domain has a timer associated with it
   */
  bool
  hasTimer() const;

  /**
   * @brief Get the domain's id
   *
   * @return int the id
   */
  int
  getId() const;
};

template<int D>
void
to_json(tpl::nlohmann::json& j, const Domain<D>& domain);

extern template class Domain<2>;
extern template class Domain<3>;

extern template void
to_json<2>(tpl::nlohmann::json& j, const Domain<2>& domain);
extern template void
to_json<3>(tpl::nlohmann::json& j, const Domain<3>& domain);

} // namespace ThunderEgg
#endif
