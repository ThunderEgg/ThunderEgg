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

#include <doctest.h>

using namespace std;
using namespace ThunderEgg;
using namespace ThunderEgg::tpl;

TEST_CASE("PatchInfo<2> getNormalNbrInfo side")
{
  for (Side<2> s : Side<2>::getValues()) {
    PatchInfo<2> pinfo;

    NormalNbrInfo<1>* nbr_info = new NormalNbrInfo<1>(2);
    pinfo.setNbrInfo(s, nbr_info);
    CHECK_EQ(&pinfo.getNormalNbrInfo(s), nbr_info);
  }
}

TEST_CASE("PatchInfo<2> getNormalNbrInfo corner")
{
  for (Corner<2> c : Corner<2>::getValues()) {
    PatchInfo<2> pinfo;

    NormalNbrInfo<0>* nbr_info = new NormalNbrInfo<0>(2);
    pinfo.setNbrInfo(c, nbr_info);
    CHECK_EQ(&pinfo.getNormalNbrInfo(c), nbr_info);
  }
}

TEST_CASE("PatchInfo<2> getNormalNbrInfo side throws on null")
{
  for (Side<2> s : Side<2>::getValues()) {
    PatchInfo<2> pinfo;

    CHECK_THROWS(pinfo.getNormalNbrInfo(s));
  }
}

TEST_CASE("PatchInfo<2> getNormalNbrInfo corner throws on null")
{
  for (Corner<2> c : Corner<2>::getValues()) {
    PatchInfo<2> pinfo;

    CHECK_THROWS(pinfo.getNormalNbrInfo(c));
  }
}

TEST_CASE("PatchInfo<2> getNormalNbrInfo side throws on wrong type")
{
  for (Side<2> s : Side<2>::getValues()) {
    for (NbrType nbr_type : { NbrType::Normal, NbrType::Coarse, NbrType::Fine }) {
      PatchInfo<2> pinfo;

      FineNbrInfo<1>* nbr_info = new FineNbrInfo<1>({ 1, 2 });
      pinfo.setNbrInfo(s, nbr_info);
      CHECK_THROWS(pinfo.getNormalNbrInfo(s));
    }
  }
}

TEST_CASE("PatchInfo<2> getNormalNbrInfo corner throws on wrong type")
{
  for (Corner<2> c : Corner<2>::getValues()) {
    PatchInfo<2> pinfo;

    FineNbrInfo<0>* nbr_info = new FineNbrInfo<0>({ 1 });
    pinfo.setNbrInfo(c, nbr_info);
    CHECK_THROWS(pinfo.getNormalNbrInfo(c));
  }
}