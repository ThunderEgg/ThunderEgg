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

#include <ThunderEgg/Domain.h>

class PatchVector
{
public:
  int max_level;
  std::vector<ThunderEgg::PatchInfo<3>> pinfos;
  std::vector<std::vector<std::vector<const ThunderEgg::PatchInfo<3>*>>> array;
  PatchVector(const ThunderEgg::Domain<3>& domain, int max_level);
  const ThunderEgg::PatchInfo<3>* operator[](const std::string& str) const;
};

std::vector<ThunderEgg::PatchInfo<3>>
GetAllPatchesOnRank0(const ThunderEgg::Domain<3>& domain);
void
Ident(int block_no, double unit_x, double unit_y, double unit_z, double& x, double& y, double& z);

void
CheckParentAndChildIdsAndRanks(const ThunderEgg::Domain<3>& coarser_domain,
                               int coarser_max_level,
                               const ThunderEgg::Domain<3>& finer_domain,
                               int finer_max_level);
void
CheckParentIdsAndRanksNull(const ThunderEgg::Domain<3>& domain);
void
CheckChildIdsAndRanksNull(const ThunderEgg::Domain<3>& domain);

void
CheckParentAndChildIdsAndRanksRefined(const ThunderEgg::Domain<3>& coarser_domain,
                                      int coarser_max_level,
                                      const ThunderEgg::Domain<3>& finer_domain,
                                      int finer_max_level);

void
CheckRootDomainNeighbors(const ThunderEgg::Domain<3>& domain);
void
Check2x2x2DomainNeighbors(const ThunderEgg::Domain<3>& domain);
void
Check2x2x2RefinedBSWDomainNeighbors(const ThunderEgg::Domain<3>& domain);
void
Check4x4x4DomainNeighbors(const ThunderEgg::Domain<3>& domain);

void
Check4x4x4RefinedBSWDomainNeighbors(const ThunderEgg::Domain<3>& domain);
