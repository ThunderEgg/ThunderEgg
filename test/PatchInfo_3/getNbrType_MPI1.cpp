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

TEST_CASE("PatchInfo<3> getNbrType side", "[PatchInfo]")
{
  for (Side<3> s : Side<3>::getValues()) {
    for (NbrType nbr_type : { NbrType::Normal, NbrType::Coarse, NbrType::Fine }) {
      PatchInfo<3> pinfo;

      NbrInfo<2>* nbr_info = nullptr;
      switch (nbr_type) {
        case NbrType::Normal:
          nbr_info = new NormalNbrInfo<2>(2);
          break;
        case NbrType::Coarse:
          nbr_info = new CoarseNbrInfo<2>(2, Orthant<2>::sw());
          break;
        case NbrType::Fine:
          nbr_info = new FineNbrInfo<2>({ 1 });
          break;
      }
      pinfo.setNbrInfo(s, nbr_info);
      CHECK(pinfo.getNbrType(s) == nbr_type);
    }
  }
}
TEST_CASE("PatchInfo<3> getNbrType edge", "[PatchInfo]")
{
  for (Edge e : Edge::getValues()) {
    for (NbrType nbr_type : { NbrType::Normal, NbrType::Coarse, NbrType::Fine }) {
      PatchInfo<3> pinfo;

      NbrInfo<1>* nbr_info = nullptr;
      switch (nbr_type) {
        case NbrType::Normal:
          nbr_info = new NormalNbrInfo<1>(2);
          break;
        case NbrType::Coarse:
          nbr_info = new CoarseNbrInfo<1>(2, Orthant<1>::lower());
          break;
        case NbrType::Fine:
          nbr_info = new FineNbrInfo<1>({ 1 });
          break;
      }
      pinfo.setNbrInfo(e, nbr_info);
      CHECK(pinfo.getNbrType(e) == nbr_type);
    }
  }
}
TEST_CASE("PatchInfo<3> getNbrType corner", "[PatchInfo]")
{
  for (Corner<3> c : Corner<3>::getValues()) {
    for (NbrType nbr_type : { NbrType::Normal, NbrType::Coarse, NbrType::Fine }) {
      PatchInfo<3> pinfo;

      NbrInfo<0>* nbr_info = nullptr;
      switch (nbr_type) {
        case NbrType::Normal:
          nbr_info = new NormalNbrInfo<0>(2);
          break;
        case NbrType::Coarse:
          nbr_info = new CoarseNbrInfo<0>(2, Orthant<0>::null());
          break;
        case NbrType::Fine:
          nbr_info = new FineNbrInfo<0>({ 1 });
          break;
      }
      pinfo.setNbrInfo(c, nbr_info);
      CHECK(pinfo.getNbrType(c) == nbr_type);
    }
  }
}
TEST_CASE("PatchInfo<3> getNbrType side throws on null", "[PatchInfo]")
{
  for (Side<3> s : Side<3>::getValues()) {
    for (NbrType nbr_type : { NbrType::Normal, NbrType::Coarse, NbrType::Fine }) {
      PatchInfo<3> pinfo;

      CHECK_THROWS(pinfo.getNbrType(s));
    }
  }
}
TEST_CASE("PatchInfo<3> getNbrType edge throws on null", "[PatchInfo]")
{
  for (Edge e : Edge::getValues()) {
    for (NbrType nbr_type : { NbrType::Normal, NbrType::Coarse, NbrType::Fine }) {
      PatchInfo<3> pinfo;

      CHECK_THROWS(pinfo.getNbrType(e));
    }
  }
}
TEST_CASE("PatchInfo<3> getNbrType corner throws on null", "[PatchInfo]")
{
  for (Corner<3> c : Corner<3>::getValues()) {
    for (NbrType nbr_type : { NbrType::Normal, NbrType::Coarse, NbrType::Fine }) {
      PatchInfo<3> pinfo;

      CHECK_THROWS(pinfo.getNbrType(c));
    }
  }
}
TEST_CASE("PatchInfo<3> getNbrType side throws on null after set", "[PatchInfo]")
{
  for (Side<3> s : Side<3>::getValues()) {
    for (NbrType nbr_type : { NbrType::Normal, NbrType::Coarse, NbrType::Fine }) {
      PatchInfo<3> pinfo;

      NbrInfo<2>* nbr_info = nullptr;
      switch (nbr_type) {
        case NbrType::Normal:
          nbr_info = new NormalNbrInfo<2>(2);
          break;
        case NbrType::Coarse:
          nbr_info = new CoarseNbrInfo<2>(2, Orthant<2>::sw());
          break;
        case NbrType::Fine:
          nbr_info = new FineNbrInfo<2>({ 1, 2 });
          break;
      }
      pinfo.setNbrInfo(s, nbr_info);
      pinfo.setNbrInfo(s, nullptr);
      CHECK_THROWS(pinfo.getNbrType(s));
    }
  }
}
TEST_CASE("PatchInfo<3> getNbrType edge throws on null after set", "[PatchInfo]")
{
  for (Edge e : Edge::getValues()) {
    for (NbrType nbr_type : { NbrType::Normal, NbrType::Coarse, NbrType::Fine }) {
      PatchInfo<3> pinfo;

      NbrInfo<1>* nbr_info = nullptr;
      switch (nbr_type) {
        case NbrType::Normal:
          nbr_info = new NormalNbrInfo<1>(2);
          break;
        case NbrType::Coarse:
          nbr_info = new CoarseNbrInfo<1>(2, Orthant<2>::lower());
          break;
        case NbrType::Fine:
          nbr_info = new FineNbrInfo<1>({ 1 });
          break;
      }
      pinfo.setNbrInfo(e, nbr_info);
      pinfo.setNbrInfo(e, nullptr);
      CHECK_THROWS(pinfo.getNbrType(e));
    }
  }
}
TEST_CASE("PatchInfo<3> getNbrType corner throws on null after set", "[PatchInfo]")
{
  for (Corner<3> c : Corner<3>::getValues()) {
    for (NbrType nbr_type : { NbrType::Normal, NbrType::Coarse, NbrType::Fine }) {
      PatchInfo<3> pinfo;

      NbrInfo<0>* nbr_info = nullptr;
      switch (nbr_type) {
        case NbrType::Normal:
          nbr_info = new NormalNbrInfo<0>(2);
          break;
        case NbrType::Coarse:
          nbr_info = new CoarseNbrInfo<0>(2, Orthant<0>::null());
          break;
        case NbrType::Fine:
          nbr_info = new FineNbrInfo<0>({ 1 });
          break;
      }
      pinfo.setNbrInfo(c, nbr_info);
      pinfo.setNbrInfo(c, nullptr);
      CHECK_THROWS(pinfo.getNbrType(c));
    }
  }
}