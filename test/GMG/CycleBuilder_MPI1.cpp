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

#include <ThunderEgg/GMG/CycleBuilder.h>

#include <catch2/catch_test_macros.hpp>

using namespace std;
using namespace ThunderEgg;

class MockOperator : public Operator<2>
{
public:
  int magic_number = 42;
  MockOperator* clone() const override { return new MockOperator(*this); }
  void apply(const Vector<2>& x, Vector<2>& b) const override {}
};
class MockSmoother : public GMG::Smoother<2>
{
public:
  int magic_number = 42;
  MockSmoother* clone() const override { return new MockSmoother(*this); }
  void smooth(const Vector<2>& x, Vector<2>& b) const override {}
};
class MockInterpolator : public GMG::Interpolator<2>
{
public:
  int magic_number = 42;
  MockInterpolator* clone() const override { return new MockInterpolator(*this); }
  void interpolate(const Vector<2>& x, Vector<2>& b) const override {}
};
class MockRestrictor : public GMG::Restrictor<2>
{
public:
  int magic_number = 42;
  MockRestrictor* clone() const override { return new MockRestrictor(*this); }
  Vector<2> restrict(const Vector<2>& x) const override { return x.getZeroClone(); }
};
TEST_CASE("CycleBuilder with two levels", "[GMG::CycleBuilder]")
{
  GMG::CycleOpts opts;
  GMG::CycleBuilder<2> builder(opts);
  std::array<MockOperator, 2> ops;
  std::array<MockSmoother, 2> smoothers;
  for (int i = 0; i < 2; i++) {
    ops[i].magic_number += i;
    smoothers[i].magic_number += i;
  }
  MockRestrictor restrictor;
  MockInterpolator interpolator;

  builder.addFinestLevel(ops[0], smoothers[0], restrictor);
  builder.addCoarsestLevel(ops[1], smoothers[1], interpolator);

  auto cycle = builder.getCycle();

  const GMG::Level<2>& finest_level = cycle->getFinestLevel();
  CHECK(dynamic_cast<const MockOperator&>(finest_level.getOperator()).magic_number == ops[0].magic_number);
  CHECK(dynamic_cast<const MockSmoother&>(finest_level.getSmoother()).magic_number == smoothers[0].magic_number);
  CHECK(dynamic_cast<const MockRestrictor&>(finest_level.getRestrictor()).magic_number == restrictor.magic_number);
  CHECK_THROWS_AS(finest_level.getInterpolator(), RuntimeError);
  CHECK(finest_level.finest());
  CHECK_FALSE(finest_level.coarsest());

  const GMG::Level<2>& coarsest_level = finest_level.getCoarser();
  CHECK(dynamic_cast<const MockOperator&>(coarsest_level.getOperator()).magic_number == ops[1].magic_number);
  CHECK(dynamic_cast<const MockSmoother&>(coarsest_level.getSmoother()).magic_number == smoothers[1].magic_number);
  CHECK_THROWS_AS(coarsest_level.getRestrictor(), RuntimeError);
  CHECK(dynamic_cast<const MockInterpolator&>(coarsest_level.getInterpolator()).magic_number == interpolator.magic_number);
  CHECK_FALSE(coarsest_level.finest());
  CHECK(coarsest_level.coarsest());
  CHECK_THROWS_AS(coarsest_level.getCoarser(), RuntimeError);
}
TEST_CASE("CycleBuilder with three levels", "[GMG::CycleBuilder]")
{
  std::array<MockOperator, 3> ops;
  std::array<MockSmoother, 3> smoothers;
  for (int i = 0; i < 3; i++) {
    ops[i].magic_number += i;
    smoothers[i].magic_number += i;
  }
  std::array<MockRestrictor, 2> restrictors;
  std::array<MockInterpolator, 2> interpolators;
  for (int i = 0; i < 2; i++) {
    restrictors[i].magic_number += i;
    interpolators[i].magic_number += i;
  }

  GMG::CycleOpts opts;
  GMG::CycleBuilder<2> builder(opts);
  builder.addFinestLevel(ops[0], smoothers[0], restrictors[0]);
  builder.addIntermediateLevel(ops[1], smoothers[1], restrictors[1], interpolators[0]);
  builder.addCoarsestLevel(ops[2], smoothers[2], interpolators[1]);

  auto cycle = builder.getCycle();

  const GMG::Level<2>& finest_level = cycle->getFinestLevel();
  CHECK(dynamic_cast<const MockOperator&>(finest_level.getOperator()).magic_number == ops[0].magic_number);
  CHECK(dynamic_cast<const MockSmoother&>(finest_level.getSmoother()).magic_number == smoothers[0].magic_number);
  CHECK(dynamic_cast<const MockRestrictor&>(finest_level.getRestrictor()).magic_number == restrictors[0].magic_number);
  CHECK_THROWS_AS(finest_level.getInterpolator(), RuntimeError);
  CHECK(finest_level.finest());
  CHECK_FALSE(finest_level.coarsest());

  const GMG::Level<2>& second_level = finest_level.getCoarser();
  CHECK(dynamic_cast<const MockOperator&>(second_level.getOperator()).magic_number == ops[1].magic_number);
  CHECK(dynamic_cast<const MockSmoother&>(second_level.getSmoother()).magic_number == smoothers[1].magic_number);
  CHECK(dynamic_cast<const MockRestrictor&>(second_level.getRestrictor()).magic_number == restrictors[1].magic_number);
  CHECK(dynamic_cast<const MockInterpolator&>(second_level.getInterpolator()).magic_number == interpolators[0].magic_number);
  CHECK_FALSE(second_level.finest());
  CHECK_FALSE(second_level.coarsest());

  const GMG::Level<2>& coarsest_level = second_level.getCoarser();
  CHECK(dynamic_cast<const MockOperator&>(coarsest_level.getOperator()).magic_number == ops[2].magic_number);
  CHECK(dynamic_cast<const MockSmoother&>(coarsest_level.getSmoother()).magic_number == smoothers[2].magic_number);
  CHECK_THROWS_AS(coarsest_level.getRestrictor(), RuntimeError);
  CHECK(dynamic_cast<const MockInterpolator&>(coarsest_level.getInterpolator()).magic_number == interpolators[1].magic_number);
  CHECK_FALSE(coarsest_level.finest());
  CHECK(coarsest_level.coarsest());
  CHECK_THROWS_AS(coarsest_level.getCoarser(), RuntimeError);
}
TEST_CASE("CycleBuilder with four levels", "[GMG::CycleBuilder]")
{
  std::array<MockOperator, 4> ops;
  std::array<MockSmoother, 4> smoothers;
  for (int i = 0; i < 4; i++) {
    ops[i].magic_number += i;
    smoothers[i].magic_number += i;
  }
  std::array<MockRestrictor, 3> restrictors;
  std::array<MockInterpolator, 3> interpolators;
  for (int i = 0; i < 3; i++) {
    restrictors[i].magic_number += i;
    interpolators[i].magic_number += i;
  }

  GMG::CycleOpts opts;
  GMG::CycleBuilder<2> builder(opts);
  builder.addFinestLevel(ops[0], smoothers[0], restrictors[0]);
  builder.addIntermediateLevel(ops[1], smoothers[1], restrictors[1], interpolators[0]);
  builder.addIntermediateLevel(ops[2], smoothers[2], restrictors[2], interpolators[1]);
  builder.addCoarsestLevel(ops[3], smoothers[3], interpolators[2]);

  auto cycle = builder.getCycle();

  const GMG::Level<2>& finest_level = cycle->getFinestLevel();
  CHECK(dynamic_cast<const MockOperator&>(finest_level.getOperator()).magic_number == ops[0].magic_number);
  CHECK(dynamic_cast<const MockSmoother&>(finest_level.getSmoother()).magic_number == smoothers[0].magic_number);
  CHECK(dynamic_cast<const MockRestrictor&>(finest_level.getRestrictor()).magic_number == restrictors[0].magic_number);
  CHECK_THROWS_AS(finest_level.getInterpolator(), RuntimeError);
  CHECK(finest_level.finest());
  CHECK_FALSE(finest_level.coarsest());

  const GMG::Level<2>& second_level = finest_level.getCoarser();
  CHECK(dynamic_cast<const MockOperator&>(second_level.getOperator()).magic_number == ops[1].magic_number);
  CHECK(dynamic_cast<const MockSmoother&>(second_level.getSmoother()).magic_number == smoothers[1].magic_number);
  CHECK(dynamic_cast<const MockRestrictor&>(second_level.getRestrictor()).magic_number == restrictors[1].magic_number);
  CHECK(dynamic_cast<const MockInterpolator&>(second_level.getInterpolator()).magic_number == interpolators[0].magic_number);
  CHECK_FALSE(second_level.finest());
  CHECK_FALSE(second_level.coarsest());

  const GMG::Level<2>& third_level = second_level.getCoarser();
  CHECK(dynamic_cast<const MockOperator&>(third_level.getOperator()).magic_number == ops[2].magic_number);
  CHECK(dynamic_cast<const MockSmoother&>(third_level.getSmoother()).magic_number == smoothers[2].magic_number);
  CHECK(dynamic_cast<const MockRestrictor&>(third_level.getRestrictor()).magic_number == restrictors[2].magic_number);
  CHECK(dynamic_cast<const MockInterpolator&>(third_level.getInterpolator()).magic_number == interpolators[1].magic_number);
  CHECK_FALSE(third_level.finest());
  CHECK_FALSE(third_level.coarsest());

  const GMG::Level<2>& coarsest_level = third_level.getCoarser();
  CHECK(dynamic_cast<const MockOperator&>(coarsest_level.getOperator()).magic_number == ops[3].magic_number);
  CHECK(dynamic_cast<const MockSmoother&>(coarsest_level.getSmoother()).magic_number == smoothers[3].magic_number);
  CHECK_THROWS_AS(coarsest_level.getRestrictor(), RuntimeError);
  CHECK(dynamic_cast<const MockInterpolator&>(coarsest_level.getInterpolator()).magic_number == interpolators[2].magic_number);
  CHECK_FALSE(coarsest_level.finest());
  CHECK(coarsest_level.coarsest());
  CHECK_THROWS_AS(coarsest_level.getCoarser(), RuntimeError);
}
TEST_CASE("CycleBuilder addFinestLevel throws exception if called twice", "[GMG::CycleBuilder]")
{
  std::array<MockOperator, 4> ops;
  std::array<MockSmoother, 4> smoothers;
  std::array<MockRestrictor, 3> restrictors;
  std::array<MockInterpolator, 3> interpolators;

  GMG::CycleOpts opts;
  GMG::CycleBuilder<2> builder(opts);
  builder.addFinestLevel(ops[0], smoothers[0], restrictors[0]);
  CHECK_THROWS_AS(builder.addFinestLevel(ops[0], smoothers[0], restrictors[0]), RuntimeError);
}
TEST_CASE("CycleBuilder addFinestLevel throws exception if called after addIntermediateLevel", "[GMG::CycleBuilder]")
{
  std::array<MockOperator, 4> ops;
  std::array<MockSmoother, 4> smoothers;
  std::array<MockRestrictor, 3> restrictors;
  std::array<MockInterpolator, 3> interpolators;

  GMG::CycleOpts opts;
  GMG::CycleBuilder<2> builder(opts);
  builder.addFinestLevel(ops[0], smoothers[0], restrictors[0]);
  builder.addIntermediateLevel(ops[1], smoothers[1], restrictors[1], interpolators[0]);
  CHECK_THROWS_AS(builder.addFinestLevel(ops[0], smoothers[0], restrictors[0]), RuntimeError);
}
TEST_CASE("CycleBuilder addFinestLevel throws exception if called after addCoarsestLevel", "[GMG::CycleBuilder]")
{
  std::array<MockOperator, 4> ops;
  std::array<MockSmoother, 4> smoothers;
  std::array<MockRestrictor, 3> restrictors;
  std::array<MockInterpolator, 3> interpolators;

  GMG::CycleOpts opts;
  GMG::CycleBuilder<2> builder(opts);
  builder.addFinestLevel(ops[0], smoothers[0], restrictors[0]);
  builder.addIntermediateLevel(ops[1], smoothers[1], restrictors[1], interpolators[0]);
  builder.addIntermediateLevel(ops[2], smoothers[2], restrictors[2], interpolators[1]);
  builder.addCoarsestLevel(ops[3], smoothers[3], interpolators[2]);
  CHECK_THROWS_AS(builder.addFinestLevel(ops[0], smoothers[0], restrictors[0]), RuntimeError);
}
TEST_CASE("CycleBuilder addIntermediateLevel throws exception if addFinestLevel isn't called", "[GMG::CycleBuilder]")
{
  std::array<MockOperator, 4> ops;
  std::array<MockSmoother, 4> smoothers;
  std::array<MockRestrictor, 3> restrictors;
  std::array<MockInterpolator, 3> interpolators;

  GMG::CycleOpts opts;
  GMG::CycleBuilder<2> builder(opts);
  CHECK_THROWS_AS(builder.addIntermediateLevel(ops[1], smoothers[1], restrictors[1], interpolators[0]), RuntimeError);
}
TEST_CASE("CycleBuilder addIntermediateLevel throws exception if called after addCoarsestLevel", "[GMG::CycleBuilder]")
{
  std::array<MockOperator, 4> ops;
  std::array<MockSmoother, 4> smoothers;
  std::array<MockRestrictor, 3> restrictors;
  std::array<MockInterpolator, 3> interpolators;

  GMG::CycleOpts opts;
  GMG::CycleBuilder<2> builder(opts);
  builder.addFinestLevel(ops[0], smoothers[0], restrictors[0]);
  builder.addIntermediateLevel(ops[1], smoothers[1], restrictors[1], interpolators[0]);
  builder.addIntermediateLevel(ops[2], smoothers[2], restrictors[2], interpolators[1]);
  builder.addCoarsestLevel(ops[3], smoothers[3], interpolators[2]);
  CHECK_THROWS_AS(builder.addIntermediateLevel(ops[2], smoothers[2], restrictors[2], interpolators[1]), RuntimeError);
}
TEST_CASE("CycleBuilder addCoarsestLevel throws exception if addFinestLevel isn't called", "[GMG::CycleBuilder]")
{
  std::array<MockOperator, 4> ops;
  std::array<MockSmoother, 4> smoothers;
  std::array<MockRestrictor, 3> restrictors;
  std::array<MockInterpolator, 3> interpolators;

  GMG::CycleOpts opts;
  GMG::CycleBuilder<2> builder(opts);
  CHECK_THROWS_AS(builder.addCoarsestLevel(ops[3], smoothers[3], interpolators[2]), RuntimeError);
}
TEST_CASE("CycleBuilder addCoarsestLevel throws exception if called twice", "[GMG::CycleBuilder]")
{
  std::array<MockOperator, 4> ops;
  std::array<MockSmoother, 4> smoothers;
  std::array<MockRestrictor, 3> restrictors;
  std::array<MockInterpolator, 3> interpolators;

  GMG::CycleOpts opts;
  GMG::CycleBuilder<2> builder(opts);
  builder.addFinestLevel(ops[0], smoothers[0], restrictors[0]);
  builder.addIntermediateLevel(ops[1], smoothers[1], restrictors[1], interpolators[0]);
  builder.addIntermediateLevel(ops[2], smoothers[2], restrictors[2], interpolators[1]);
  builder.addCoarsestLevel(ops[3], smoothers[3], interpolators[2]);
  CHECK_THROWS_AS(builder.addCoarsestLevel(ops[3], smoothers[3], interpolators[2]), RuntimeError);
}
TEST_CASE("CycleBuilder getCycle throws exception if addCoarsestLevel isn't called", "[GMG::CycleBuilder]")
{
  std::array<MockOperator, 4> ops;
  std::array<MockSmoother, 4> smoothers;
  std::array<MockRestrictor, 3> restrictors;
  std::array<MockInterpolator, 3> interpolators;

  GMG::CycleOpts opts;
  GMG::CycleBuilder<2> builder(opts);
  builder.addFinestLevel(ops[0], smoothers[0], restrictors[0]);
  builder.addIntermediateLevel(ops[1], smoothers[1], restrictors[1], interpolators[0]);
  builder.addIntermediateLevel(ops[2], smoothers[2], restrictors[2], interpolators[1]);

  CHECK_THROWS_AS(builder.getCycle(), RuntimeError);
}
TEST_CASE("CycleBuilder getCycle throws exception if no calls are made", "[GMG::CycleBuilder]")
{
  GMG::CycleOpts opts;
  GMG::CycleBuilder<2> builder(opts);

  CHECK_THROWS_AS(builder.getCycle(), RuntimeError);
}
