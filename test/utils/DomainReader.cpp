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
#include "DomainReader.h"

#include <ThunderEgg/tpl/json.hpp>
#include <fstream>

template<int D>
ThunderEgg::PatchInfo<D>
DomainReader<D>::parsePatch(ThunderEgg::tpl::nlohmann::json& patch_j)
{
  ThunderEgg::PatchInfo<D> pinfo = patch_j.get<ThunderEgg::PatchInfo<D>>();
  pinfo.num_ghost_cells = num_ghost;
  pinfo.ns = ns;
  for (size_t d = 0; d < D; d++) {
    pinfo.spacings[d] /= ns[d];
  }
  return pinfo;
}

template<int D>
DomainReader<D>::DomainReader(std::string file_name, std::array<int, D> ns_in, int num_ghost_in)
  : ns(ns_in)
  , num_ghost(num_ghost_in)
{
  ThunderEgg::Communicator comm(MPI_COMM_WORLD);

  ThunderEgg::tpl::nlohmann::json j;
  std::ifstream input_stream(ThunderEgg::TEST_SRC_DIR + file_name);
  if (!input_stream.good()) {
    throw "could not open file";
  }
  input_stream >> j;
  input_stream.close();
  std::deque<ThunderEgg::PatchInfo<D>> finer_patches;
  for (ThunderEgg::tpl::nlohmann::json& patch_j : j.at("levels")[0]) {
    auto patch = parsePatch(patch_j);
    if (patch.rank == comm.getRank())
      finer_patches.push_back(patch);
  }
  finer_domain = std::make_shared<ThunderEgg::Domain<D>>(
    comm, 0, ns, num_ghost, finer_patches.begin(), finer_patches.end());
  std::deque<ThunderEgg::PatchInfo<D>> coarser_patches;
  for (ThunderEgg::tpl::nlohmann::json& patch_j : j.at("levels")[1]) {
    auto patch = parsePatch(patch_j);
    if (patch.rank == comm.getRank())
      coarser_patches.push_back(patch);
  }
  coarser_domain = std::make_shared<ThunderEgg::Domain<D>>(
    comm, 1, ns, num_ghost, coarser_patches.begin(), coarser_patches.end());
}

template<int D>
ThunderEgg::Domain<D>
DomainReader<D>::getCoarserDomain()
{
  return *coarser_domain;
}

template<int D>
ThunderEgg::Domain<D>
DomainReader<D>::getFinerDomain()
{
  return *finer_domain;
}

template class DomainReader<2>;
template class DomainReader<3>;