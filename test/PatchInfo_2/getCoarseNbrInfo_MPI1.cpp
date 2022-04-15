/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured
 *  Cartesian grids.
 *
 *  Copyright (c) 2020-2022 Scott Aiton
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
#include <ThunderEgg/PatchInfo.h>

#include <catch2/catch_test_macros.hpp>

using namespace std;
using namespace ThunderEgg;
using namespace ThunderEgg::tpl;

TEST_CASE("PatchInfo<2> getCoarseNbrInfo side", "[PatchInfo]")
{
  for (Side<2> s : Side<2>::getValues()) {
    PatchInfo<2> pinfo;

    CoarseNbrInfo<1>* nbr_info = new CoarseNbrInfo<1>(2, Orthant<1>::lower());
    pinfo.setNbrInfo(s, nbr_info);
    CHECK(&pinfo.getCoarseNbrInfo(s) == nbr_info);
  }
}
TEST_CASE("PatchInfo<2> getCoarseNbrInfo corner", "[PatchInfo]")
{
  for (Corner<2> c : Corner<2>::getValues()) {
    PatchInfo<2> pinfo;

    CoarseNbrInfo<0>* nbr_info = new CoarseNbrInfo<0>(2, Orthant<0>::null());
    pinfo.setNbrInfo(c, nbr_info);
    CHECK(&pinfo.getCoarseNbrInfo(c) == nbr_info);
  }
}
TEST_CASE("PatchInfo<2> getCoarseNbrInfo side throws on null", "[PatchInfo]")
{
  for (Side<2> s : Side<2>::getValues()) {
    PatchInfo<2> pinfo;

    CHECK_THROWS(pinfo.getCoarseNbrInfo(s));
  }
}
TEST_CASE("PatchInfo<2> getCoarseNbrInfo corner throws on null", "[PatchInfo]")
{
  for (Corner<2> c : Corner<2>::getValues()) {
    PatchInfo<2> pinfo;

    CHECK_THROWS(pinfo.getCoarseNbrInfo(c));
  }
}
TEST_CASE("PatchInfo<2> getCoarseNbrInfo side throws on wrong type", "[PatchInfo]")
{
  for (Side<2> s : Side<2>::getValues()) {
    PatchInfo<2> pinfo;

    FineNbrInfo<1>* nbr_info = new FineNbrInfo<1>({ 1, 2 });
    pinfo.setNbrInfo(s, nbr_info);
    CHECK_THROWS(pinfo.getCoarseNbrInfo(s));
  }
}
TEST_CASE("PatchInfo<2> getCoarseNbrInfo corner throws on wrong type", "[PatchInfo]")
{
  for (Corner<2> c : Corner<2>::getValues()) {
    PatchInfo<2> pinfo;

    FineNbrInfo<0>* nbr_info = new FineNbrInfo<0>({ 1 });
    pinfo.setNbrInfo(c, nbr_info);
    CHECK_THROWS(pinfo.getCoarseNbrInfo(c));
  }
}