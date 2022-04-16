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

#include <ThunderEgg/Schur/IfaceType.h>

#include <set>

#include <catch2/catch_test_macros.hpp>

using namespace std;
using namespace ThunderEgg;

TEST_CASE("Schur::IfaceType default constructor", "[Schur::IfaceType]")
{
  Schur::IfaceType<2> type;
  CHECK_FALSE(type.isNormal());
  CHECK_FALSE(type.isCoarseToCoarse());
  CHECK_FALSE(type.isCoarseToFine());
  CHECK_FALSE(type.isFineToFine());
  CHECK_FALSE(type.isFineToCoarse());
  CHECK(type.getOrthant() == Orthant<1>::null());
}
TEST_CASE("Schur::IfaceType Normal constructor", "[Schur::IfaceType]")
{
  Schur::IfaceType<2> type = Schur::IfaceType<2>::Normal();
  CHECK(type.isNormal());
  CHECK_FALSE(type.isCoarseToCoarse());
  CHECK_FALSE(type.isCoarseToFine());
  CHECK_FALSE(type.isFineToFine());
  CHECK_FALSE(type.isFineToCoarse());
  CHECK(type.getOrthant() == Orthant<1>::null());
}
TEST_CASE("Schur::IfaceType CoarseToCoarse constructor", "[Schur::IfaceType]")
{
  Schur::IfaceType<2> type = Schur::IfaceType<2>::CoarseToCoarse();
  CHECK_FALSE(type.isNormal());
  CHECK(type.isCoarseToCoarse());
  CHECK_FALSE(type.isCoarseToFine());
  CHECK_FALSE(type.isFineToFine());
  CHECK_FALSE(type.isFineToCoarse());
  CHECK(type.getOrthant() == Orthant<1>::null());
}
TEST_CASE("Schur::IfaceType CoarseToFine constructor", "[Schur::IfaceType]")
{
  for (Orthant<1> orth : Orthant<1>::getValues()) {
    Schur::IfaceType<2> type = Schur::IfaceType<2>::CoarseToFine(orth);
    CHECK_FALSE(type.isNormal());
    CHECK_FALSE(type.isCoarseToCoarse());
    CHECK(type.isCoarseToFine());
    CHECK_FALSE(type.isFineToFine());
    CHECK_FALSE(type.isFineToCoarse());
    CHECK(type.getOrthant() == orth);
  }
}
TEST_CASE("Schur::IfaceType FineToFine constructor", "[Schur::IfaceType]")
{
  for (Orthant<1> orth : Orthant<1>::getValues()) {
    Schur::IfaceType<2> type = Schur::IfaceType<2>::FineToFine(orth);
    CHECK_FALSE(type.isNormal());
    CHECK_FALSE(type.isCoarseToCoarse());
    CHECK_FALSE(type.isCoarseToFine());
    CHECK(type.isFineToFine());
    CHECK_FALSE(type.isFineToCoarse());
    CHECK(type.getOrthant() == orth);
  }
}
TEST_CASE("Schur::IfaceType FineToCoarse constructor", "[Schur::IfaceType]")
{
  for (Orthant<1> orth : Orthant<1>::getValues()) {
    Schur::IfaceType<2> type = Schur::IfaceType<2>::FineToCoarse(orth);
    CHECK_FALSE(type.isNormal());
    CHECK_FALSE(type.isCoarseToCoarse());
    CHECK_FALSE(type.isCoarseToFine());
    CHECK_FALSE(type.isFineToFine());
    CHECK(type.isFineToCoarse());
    CHECK(type.getOrthant() == orth);
  }
}
TEST_CASE("Schur::IfaceType <", "[Schur::IfaceType]")
{
  // this is an indirect test
  set<Schur::IfaceType<2>> ifaces;
  ifaces.insert(Schur::IfaceType<2>());
  ifaces.insert(Schur::IfaceType<2>::Normal());
  ifaces.insert(Schur::IfaceType<2>::CoarseToCoarse());
  for (Orthant<1> orth : Orthant<1>::getValues()) {
    ifaces.insert(Schur::IfaceType<2>::CoarseToFine(orth));
    ifaces.insert(Schur::IfaceType<2>::FineToFine(orth));
    ifaces.insert(Schur::IfaceType<2>::FineToCoarse(orth));
  }
  CHECK(ifaces.size() == 9);
}
TEST_CASE("Schur::IfaceType setOrthant", "[Schur::IfaceType]")
{
  Schur::IfaceType<2> type = Schur::IfaceType<2>::CoarseToFine(Orthant<1>::null());
  for (Orthant<1> orth : Orthant<1>::getValues()) {
    type.setOrthant(orth);
    CHECK(type.getOrthant() == orth);
  }
}
