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

#include <ThunderEgg/Schur/NormalIfaceInfo.h>

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
TEST_CASE("Schur::NormalIfaceInfo constructor")
{
  for (Side<2> s : Side<2>::getValues()) {
    int id = 1;
    int nbr_id = 2;
    PatchInfo<2> pinfo;
    pinfo.rank = 0;
    pinfo.id = id;
    pinfo.setNbrInfo(s, new NormalNbrInfo<1>(nbr_id));
    pinfo.getNormalNbrInfo(s).rank = 1;
    Schur::NormalIfaceInfo<2> iface_info(pinfo, s);
    if (s.isHigherOnAxis()) {
      // check that iface belongs to this
      CHECK_EQ(iface_info.rank, 0);
    } else {
      // check that iface belongs to nbr
      CHECK_EQ(iface_info.rank, 1);
    }
    // check that the id is encoded as expected
    if (s.isHigherOnAxis()) {
      // check that iface belongs to this
      CHECK_EQ(iface_info.id / (int)Side<2>::number_of, id);
      CHECK_EQ(iface_info.id % Side<2>::number_of, s.getIndex());
    } else {
      // check that iface belongs to nbr
      CHECK_EQ(iface_info.id / (int)Side<2>::number_of, nbr_id);
      CHECK_EQ(iface_info.id % Side<2>::number_of, s.opposite().getIndex());
    }
    // local and global index should be set to -1
    CHECK_EQ(iface_info.patch_local_index, -1);
    CHECK_EQ(iface_info.row_local_index, -1);
    CHECK_EQ(iface_info.col_local_index, -1);
    CHECK_EQ(iface_info.global_index, -1);
  }
}
