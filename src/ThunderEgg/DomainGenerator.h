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

#ifndef THUNDEREGG_DOMAINGENERATOR_H
#define THUNDEREGG_DOMAINGENERATOR_H
/**
 * @file
 *
 * @brief DomainGenerator class
 */
#include <ThunderEgg/Domain.h>
namespace ThunderEgg {
/**
 * @brief Generates Domain objects.
 *
 * This class is intended to wrap around octree/quadtree libraries and provide ThunderEgg with the
 * necessary patch information. See P4estDG for the implimentation with the p4est library.
 *
 * @tparam D the number of Cartesian dimensions
 */
template<int D>
class DomainGenerator
{
public:
  /**
   * @brief Destroy the DomainGenerator object
   */
  virtual ~DomainGenerator(){};
  /**
   * @brief Return the finest domain
   */
  virtual Domain<D> getFinestDomain() = 0;
  /**
   * @brief return true if there is a coarser domain to be generated.
   */
  virtual bool hasCoarserDomain() = 0;
  /**
   * @brief Return a new coarser domain
   */
  virtual Domain<D> getCoarserDomain() = 0;
};
} // namespace ThunderEgg
#endif
