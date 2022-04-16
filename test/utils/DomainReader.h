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
#include <TestConfig.h>
#include <ThunderEgg/Domain.h>

template<int D>
class DomainReader
{
private:
  std::array<int, D> ns;
  int num_ghost;
  std::shared_ptr<ThunderEgg::Domain<D>> coarser_domain;
  std::shared_ptr<ThunderEgg::Domain<D>> finer_domain;
  ThunderEgg::PatchInfo<D>
  parsePatch(ThunderEgg::tpl::nlohmann::json& patch_j);

public:
  DomainReader(std::string file_name, std::array<int, D> ns_in, int num_ghost_in);
  ThunderEgg::Domain<D>
  getCoarserDomain();
  ThunderEgg::Domain<D>
  getFinerDomain();
};
extern template class DomainReader<2>;
extern template class DomainReader<3>;