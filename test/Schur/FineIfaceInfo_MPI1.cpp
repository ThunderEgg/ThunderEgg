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

#include <ThunderEgg/Schur/FineIfaceInfo.h>
#include <algorithm>

#include <doctest.h>

using namespace std;
using namespace ThunderEgg;

template<typename Container, typename Value>
bool
contains(Container& deque, Value a)
{
  return find(deque.begin(), deque.end(), a) != deque.end();
}
TEST_CASE("Schur::FineIfaceInfo constructor")
{
  for (Side<2> s : Side<2>::getValues()) {
    int id = 1;
    array<int, 2> nbr_ids = { 2, 3 };
    PatchInfo<2> pinfo;
    pinfo.rank = 0;
    pinfo.id = id;
    pinfo.setNbrInfo(s, new FineNbrInfo<1>(nbr_ids));
    pinfo.getFineNbrInfo(s).ranks[0] = 1;
    pinfo.getFineNbrInfo(s).ranks[1] = 2;
    Schur::FineIfaceInfo<2> iface_info(pinfo, s);
    CHECK_EQ(iface_info.rank, 0);
    CHECK_EQ(iface_info.fine_ranks[0], 1);
    CHECK_EQ(iface_info.fine_ranks[1], 2);
    // check that the id is encoded as expected
    CHECK_EQ(iface_info.id / (int)Side<2>::number_of, id);
    CHECK_EQ(iface_info.id % Side<2>::number_of, s.getIndex());
    // check that iface belongs to nbr
    CHECK_EQ(iface_info.fine_ids[0] / (int)Side<2>::number_of, nbr_ids[0]);
    CHECK_EQ(iface_info.fine_ids[0] % Side<2>::number_of, s.opposite().getIndex());
    CHECK_EQ(iface_info.fine_ids[1] / (int)Side<2>::number_of, nbr_ids[1]);
    CHECK_EQ(iface_info.fine_ids[1] % Side<2>::number_of, s.opposite().getIndex());
    // local and global index should be set to -1
    CHECK_EQ(iface_info.patch_local_index, -1);
    CHECK_EQ(iface_info.row_local_index, -1);
    CHECK_EQ(iface_info.col_local_index, -1);
    CHECK_EQ(iface_info.fine_col_local_indexes[0], -1);
    CHECK_EQ(iface_info.fine_col_local_indexes[1], -1);
    CHECK_EQ(iface_info.global_index, -1);
    CHECK_EQ(iface_info.fine_global_indexes[0], -1);
    CHECK_EQ(iface_info.fine_global_indexes[1], -1);
  }
}
