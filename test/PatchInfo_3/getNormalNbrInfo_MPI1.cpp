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

TEST_CASE("PatchInfo<3> getNormalNbrInfo side")
{
  for (Side<3> e : Side<3>::getValues()) {
    PatchInfo<3> pinfo;

    NormalNbrInfo<2>* nbr_info = new NormalNbrInfo<2>(2);
    pinfo.setNbrInfo(e, nbr_info);
    CHECK_EQ(&pinfo.getNormalNbrInfo(e), nbr_info);
  }
}

TEST_CASE("PatchInfo<3> getNormalNbrInfo edge")
{
  for (Edge e : Edge::getValues()) {
    PatchInfo<3> pinfo;

    NormalNbrInfo<1>* nbr_info = new NormalNbrInfo<1>(2);
    pinfo.setNbrInfo(e, nbr_info);
    CHECK_EQ(&pinfo.getNormalNbrInfo(e), nbr_info);
  }
}

TEST_CASE("PatchInfo<3> getNormalNbrInfo corner")
{
  for (Corner<3> c : Corner<3>::getValues()) {
    PatchInfo<3> pinfo;

    NormalNbrInfo<0>* nbr_info = new NormalNbrInfo<0>(2);
    pinfo.setNbrInfo(c, nbr_info);
    CHECK_EQ(&pinfo.getNormalNbrInfo(c), nbr_info);
  }
}

TEST_CASE("PatchInfo<3> getNormalNbrInfo side throws on null")
{
  for (Side<3> s : Side<3>::getValues()) {
    PatchInfo<3> pinfo;

    CHECK_THROWS(pinfo.getNormalNbrInfo(s));
  }
}

TEST_CASE("PatchInfo<3> getNormalNbrInfo edge throws on null")
{
  for (Edge e : Edge::getValues()) {
    PatchInfo<3> pinfo;

    CHECK_THROWS(pinfo.getNormalNbrInfo(e));
  }
}

TEST_CASE("PatchInfo<3> getNormalNbrInfo corner throws on null")
{
  for (Corner<3> c : Corner<3>::getValues()) {
    PatchInfo<3> pinfo;

    CHECK_THROWS(pinfo.getNormalNbrInfo(c));
  }
}

TEST_CASE("PatchInfo<3> getNormalNbrInfo side throws on wrong type")
{
  for (Side<3> s : Side<3>::getValues()) {
    PatchInfo<3> pinfo;

    FineNbrInfo<2>* nbr_info = new FineNbrInfo<2>({ 1, 2, 3, 4 });
    pinfo.setNbrInfo(s, nbr_info);
    CHECK_THROWS(pinfo.getNormalNbrInfo(s));
  }
}

TEST_CASE("PatchInfo<3> getNormalNbrInfo edge throws on wrong type")
{
  for (Edge e : Edge::getValues()) {
    PatchInfo<3> pinfo;

    FineNbrInfo<1>* nbr_info = new FineNbrInfo<1>({ 1, 2 });
    pinfo.setNbrInfo(e, nbr_info);
    CHECK_THROWS(pinfo.getNormalNbrInfo(e));
  }
}

TEST_CASE("PatchInfo<3> getNormalNbrInfo corner throws on wrong type")
{
  for (Corner<3> c : Corner<3>::getValues()) {
    PatchInfo<3> pinfo;

    FineNbrInfo<0>* nbr_info = new FineNbrInfo<0>({ 1 });
    pinfo.setNbrInfo(c, nbr_info);
    CHECK_THROWS(pinfo.getNormalNbrInfo(c));
  }
}