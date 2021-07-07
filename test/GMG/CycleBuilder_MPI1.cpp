/***************************************************************************
 *  ThunderEgg, a library for solving Poisson's equation on adaptively
 *  refined block-structured Cartesian grids
 *
 *  Copyright (C) 2019  ThunderEgg Developers. See AUTHORS.md file at the
 *  top-level directory.
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
	void apply(const Vector<2> &x, Vector<2> &b) const override {}
};
class MockSmoother : public GMG::Smoother<2>
{
	public:
	int           magic_number = 42;
	MockSmoother *clone() const override
	{
		return new MockSmoother(*this);
	}
	void smooth(const Vector<2> &x, Vector<2> &b) const override {}
};
class MockInterpolator : public GMG::Interpolator<2>
{
	public:
	int               magic_number = 42;
	MockInterpolator *clone() const override
	{
		return new MockInterpolator(*this);
	}
	void interpolate(const Vector<2> &x, Vector<2> &b) const override
	{
	}
};
class MockRestrictor : public GMG::Restrictor<2>
{
	public:
	int             magic_number = 42;
	MockRestrictor *clone() const override
	{
		return new MockRestrictor(*this);
	}
	Vector<2> restrict(const Vector<2> &x) const override
	{
		return x.getZeroClone();
	}
	Vector<2> getNewCoarserVector(int num_components) const override
	{
		return Vector<2>();
	}
};
TEST_CASE("CycleBuilder with two levels", "[GMG::CycleBuilder]")
{
	GMG::CycleOpts                               opts;
	GMG::CycleBuilder<2>                         builder(opts);
	std::array<std::shared_ptr<MockOperator>, 2> ops;
	std::array<MockSmoother, 2>                  smoothers;
	for (int i = 0; i < 2; i++) {
		ops[i] = make_shared<MockOperator>();
		smoothers[i].magic_number += i;
	}
	MockRestrictor   restrictor;
	MockInterpolator interpolator;

	builder.addFinestLevel(ops[0], smoothers[0], restrictor);
	builder.addCoarsestLevel(ops[1], smoothers[1], interpolator);

	auto cycle = builder.getCycle();

	auto finest_level = cycle->getFinestLevel();
	CHECK(finest_level->getOperator() == ops[0]);
	CHECK(dynamic_cast<const MockSmoother &>(finest_level->getSmoother()).magic_number == smoothers[0].magic_number);
	CHECK(dynamic_cast<const MockRestrictor &>(finest_level->getRestrictor()).magic_number == restrictor.magic_number);
	CHECK_THROWS_AS(finest_level->getInterpolator(), RuntimeError);
	CHECK(finest_level->finest());
	CHECK_FALSE(finest_level->coarsest());
	CHECK_THROWS_AS(finest_level->getFiner(), std::bad_weak_ptr);
	CHECK(finest_level->getCoarser() != nullptr);

	auto coarsest_level = finest_level->getCoarser();
	CHECK(coarsest_level->getOperator() == ops[1]);
	CHECK(dynamic_cast<const MockSmoother &>(coarsest_level->getSmoother()).magic_number == smoothers[1].magic_number);
	CHECK_THROWS_AS(coarsest_level->getRestrictor(), RuntimeError);
	CHECK(dynamic_cast<const MockInterpolator &>(coarsest_level->getInterpolator()).magic_number == interpolator.magic_number);
	CHECK_FALSE(coarsest_level->finest());
	CHECK(coarsest_level->coarsest());
	CHECK(coarsest_level->getFiner() == finest_level);
	CHECK(coarsest_level->getCoarser() == nullptr);
}
TEST_CASE("CycleBuilder with three levels", "[GMG::CycleBuilder]")
{
	std::array<std::shared_ptr<MockOperator>, 3> ops;
	std::array<MockSmoother, 3>                  smoothers;
	for (int i = 0; i < 3; i++) {
		ops[i] = make_shared<MockOperator>();
		smoothers[i].magic_number += i;
	}
	std::array<MockRestrictor, 2>   restrictors;
	std::array<MockInterpolator, 2> interpolators;
	for (int i = 0; i < 2; i++) {
		restrictors[i].magic_number += i;
		interpolators[i].magic_number += i;
	}

	GMG::CycleOpts       opts;
	GMG::CycleBuilder<2> builder(opts);
	builder.addFinestLevel(ops[0], smoothers[0], restrictors[0]);
	builder.addIntermediateLevel(ops[1], smoothers[1], restrictors[1], interpolators[0]);
	builder.addCoarsestLevel(ops[2], smoothers[2], interpolators[1]);

	auto cycle = builder.getCycle();

	auto finest_level = cycle->getFinestLevel();
	CHECK(finest_level->getOperator() == ops[0]);
	CHECK(dynamic_cast<const MockSmoother &>(finest_level->getSmoother()).magic_number == smoothers[0].magic_number);
	CHECK(dynamic_cast<const MockRestrictor &>(finest_level->getRestrictor()).magic_number == restrictors[0].magic_number);
	CHECK_THROWS_AS(finest_level->getInterpolator(), RuntimeError);
	CHECK(finest_level->finest());
	CHECK_FALSE(finest_level->coarsest());
	CHECK_THROWS_AS(finest_level->getFiner(), std::bad_weak_ptr);
	CHECK(finest_level->getCoarser() != nullptr);

	auto second_level = finest_level->getCoarser();
	CHECK(second_level->getOperator() == ops[1]);
	CHECK(dynamic_cast<const MockSmoother &>(second_level->getSmoother()).magic_number == smoothers[1].magic_number);
	CHECK(dynamic_cast<const MockRestrictor &>(second_level->getRestrictor()).magic_number == restrictors[1].magic_number);
	CHECK(dynamic_cast<const MockInterpolator &>(second_level->getInterpolator()).magic_number == interpolators[0].magic_number);
	CHECK_FALSE(second_level->finest());
	CHECK_FALSE(second_level->coarsest());
	CHECK(second_level->getFiner() == finest_level);
	CHECK(second_level->getCoarser() != nullptr);

	auto coarsest_level = second_level->getCoarser();
	CHECK(coarsest_level->getOperator() == ops[2]);
	CHECK(dynamic_cast<const MockSmoother &>(coarsest_level->getSmoother()).magic_number == smoothers[2].magic_number);
	CHECK_THROWS_AS(coarsest_level->getRestrictor(), RuntimeError);
	CHECK(dynamic_cast<const MockInterpolator &>(coarsest_level->getInterpolator()).magic_number == interpolators[1].magic_number);
	CHECK_FALSE(coarsest_level->finest());
	CHECK(coarsest_level->coarsest());
	CHECK(coarsest_level->getFiner() == second_level);
	CHECK(coarsest_level->getCoarser() == nullptr);
}
TEST_CASE("CycleBuilder with four levels", "[GMG::CycleBuilder]")
{
	std::array<std::shared_ptr<MockOperator>, 4> ops;
	std::array<MockSmoother, 4>                  smoothers;
	for (int i = 0; i < 4; i++) {
		ops[i] = make_shared<MockOperator>();
		smoothers[i].magic_number += i;
	}
	std::array<MockRestrictor, 3>   restrictors;
	std::array<MockInterpolator, 3> interpolators;
	for (int i = 0; i < 3; i++) {
		restrictors[i].magic_number += i;
		interpolators[i].magic_number += i;
	}

	GMG::CycleOpts       opts;
	GMG::CycleBuilder<2> builder(opts);
	builder.addFinestLevel(ops[0], smoothers[0], restrictors[0]);
	builder.addIntermediateLevel(ops[1], smoothers[1], restrictors[1], interpolators[0]);
	builder.addIntermediateLevel(ops[2], smoothers[2], restrictors[2], interpolators[1]);
	builder.addCoarsestLevel(ops[3], smoothers[3], interpolators[2]);

	auto cycle = builder.getCycle();

	auto finest_level = cycle->getFinestLevel();
	CHECK(finest_level->getOperator() == ops[0]);
	CHECK(dynamic_cast<const MockSmoother &>(finest_level->getSmoother()).magic_number == smoothers[0].magic_number);
	CHECK(dynamic_cast<const MockRestrictor &>(finest_level->getRestrictor()).magic_number == restrictors[0].magic_number);
	CHECK_THROWS_AS(finest_level->getInterpolator(), RuntimeError);
	CHECK(finest_level->finest());
	CHECK_FALSE(finest_level->coarsest());
	CHECK_THROWS_AS(finest_level->getFiner(), std::bad_weak_ptr);
	CHECK(finest_level->getCoarser() != nullptr);

	auto second_level = finest_level->getCoarser();
	CHECK(second_level->getOperator() == ops[1]);
	CHECK(dynamic_cast<const MockSmoother &>(second_level->getSmoother()).magic_number == smoothers[1].magic_number);
	CHECK(dynamic_cast<const MockRestrictor &>(second_level->getRestrictor()).magic_number == restrictors[1].magic_number);
	CHECK(dynamic_cast<const MockInterpolator &>(second_level->getInterpolator()).magic_number == interpolators[0].magic_number);
	CHECK_FALSE(second_level->finest());
	CHECK_FALSE(second_level->coarsest());
	CHECK(second_level->getFiner() == finest_level);
	CHECK(second_level->getCoarser() != nullptr);

	auto third_level = second_level->getCoarser();
	CHECK(third_level->getOperator() == ops[2]);
	CHECK(dynamic_cast<const MockSmoother &>(third_level->getSmoother()).magic_number == smoothers[2].magic_number);
	CHECK(dynamic_cast<const MockRestrictor &>(third_level->getRestrictor()).magic_number == restrictors[2].magic_number);
	CHECK(dynamic_cast<const MockInterpolator &>(third_level->getInterpolator()).magic_number == interpolators[1].magic_number);
	CHECK_FALSE(third_level->finest());
	CHECK_FALSE(third_level->coarsest());
	CHECK(third_level->getFiner() == second_level);
	CHECK(third_level->getCoarser() != nullptr);

	auto coarsest_level = third_level->getCoarser();
	CHECK(coarsest_level->getOperator() == ops[3]);
	CHECK(dynamic_cast<const MockSmoother &>(coarsest_level->getSmoother()).magic_number == smoothers[3].magic_number);
	CHECK_THROWS_AS(coarsest_level->getRestrictor(), RuntimeError);
	CHECK(dynamic_cast<const MockInterpolator &>(coarsest_level->getInterpolator()).magic_number == interpolators[2].magic_number);
	CHECK_FALSE(coarsest_level->finest());
	CHECK(coarsest_level->coarsest());
	CHECK(coarsest_level->getFiner() == third_level);
	CHECK(coarsest_level->getCoarser() == nullptr);
}
TEST_CASE("CycleBuilder addFinestLevel throws exception with nullptr for operator",
          "[GMG::CycleBuilder]")
{
	std::array<std::shared_ptr<MockOperator>, 4> ops;
	std::array<MockSmoother, 4>                  smoothers;
	for (int i = 0; i < 4; i++) {
		ops[i] = make_shared<MockOperator>();
		smoothers[i].magic_number += i;
	}
	std::array<MockRestrictor, 3>   restrictors;
	std::array<MockInterpolator, 3> interpolators;

	GMG::CycleOpts       opts;
	GMG::CycleBuilder<2> builder(opts);
	CHECK_THROWS_AS(builder.addFinestLevel(nullptr, smoothers[0], restrictors[0]),
	                RuntimeError);
}
TEST_CASE("CycleBuilder addIntermediateLevel throws exception with nullptr for operator",
          "[GMG::CycleBuilder]")
{
	std::array<std::shared_ptr<MockOperator>, 4> ops;
	std::array<MockSmoother, 4>                  smoothers;
	for (int i = 0; i < 4; i++) {
		ops[i] = make_shared<MockOperator>();
		smoothers[i].magic_number += i;
	}
	std::array<MockRestrictor, 3>   restrictors;
	std::array<MockInterpolator, 3> interpolators;

	GMG::CycleOpts       opts;
	GMG::CycleBuilder<2> builder(opts);
	builder.addFinestLevel(ops[0], smoothers[0], restrictors[0]);
	CHECK_THROWS_AS(
	builder.addIntermediateLevel(nullptr, smoothers[1], restrictors[1], interpolators[0]),
	RuntimeError);
}
TEST_CASE("CycleBuilder addCoarsestLevel throws exception with nullptr for operator",
          "[GMG::CycleBuilder]")
{
	std::array<std::shared_ptr<MockOperator>, 4> ops;
	std::array<MockSmoother, 4>                  smoothers;
	for (int i = 0; i < 4; i++) {
		ops[i] = make_shared<MockOperator>();
		smoothers[i].magic_number += i;
	}
	std::array<MockRestrictor, 3>   restrictors;
	std::array<MockInterpolator, 3> interpolators;

	GMG::CycleOpts       opts;
	GMG::CycleBuilder<2> builder(opts);
	builder.addFinestLevel(ops[0], smoothers[0], restrictors[0]);
	builder.addIntermediateLevel(ops[1], smoothers[1], restrictors[1], interpolators[0]);
	builder.addIntermediateLevel(ops[2], smoothers[2], restrictors[2], interpolators[1]);
	CHECK_THROWS_AS(builder.addCoarsestLevel(nullptr, smoothers[3], interpolators[2]),
	                RuntimeError);
}
TEST_CASE("CycleBuilder addFinestLevel throws exception if called twice", "[GMG::CycleBuilder]")
{
	std::array<std::shared_ptr<MockOperator>, 4> ops;
	std::array<MockSmoother, 4>                  smoothers;
	for (int i = 0; i < 4; i++) {
		ops[i] = make_shared<MockOperator>();
		smoothers[i].magic_number += i;
	}
	std::array<MockRestrictor, 3>   restrictors;
	std::array<MockInterpolator, 3> interpolators;

	GMG::CycleOpts       opts;
	GMG::CycleBuilder<2> builder(opts);
	builder.addFinestLevel(ops[0], smoothers[0], restrictors[0]);
	CHECK_THROWS_AS(builder.addFinestLevel(ops[0], smoothers[0], restrictors[0]),
	                RuntimeError);
}
TEST_CASE("CycleBuilder addFinestLevel throws exception if called after addIntermediateLevel",
          "[GMG::CycleBuilder]")
{
	std::array<std::shared_ptr<MockOperator>, 4> ops;
	std::array<MockSmoother, 4>                  smoothers;
	for (int i = 0; i < 4; i++) {
		ops[i] = make_shared<MockOperator>();
		smoothers[i].magic_number += i;
	}
	std::array<MockRestrictor, 3>   restrictors;
	std::array<MockInterpolator, 3> interpolators;

	GMG::CycleOpts       opts;
	GMG::CycleBuilder<2> builder(opts);
	builder.addFinestLevel(ops[0], smoothers[0], restrictors[0]);
	builder.addIntermediateLevel(ops[1], smoothers[1], restrictors[1], interpolators[0]);
	CHECK_THROWS_AS(builder.addFinestLevel(ops[0], smoothers[0], restrictors[0]),
	                RuntimeError);
}
TEST_CASE("CycleBuilder addFinestLevel throws exception if called after addCoarsestLevel",
          "[GMG::CycleBuilder]")
{
	std::array<std::shared_ptr<MockOperator>, 4> ops;
	std::array<MockSmoother, 4>                  smoothers;
	for (int i = 0; i < 4; i++) {
		ops[i] = make_shared<MockOperator>();
		smoothers[i].magic_number += i;
	}
	std::array<MockRestrictor, 3>   restrictors;
	std::array<MockInterpolator, 3> interpolators;

	GMG::CycleOpts       opts;
	GMG::CycleBuilder<2> builder(opts);
	builder.addFinestLevel(ops[0], smoothers[0], restrictors[0]);
	builder.addIntermediateLevel(ops[1], smoothers[1], restrictors[1], interpolators[0]);
	builder.addIntermediateLevel(ops[2], smoothers[2], restrictors[2], interpolators[1]);
	builder.addCoarsestLevel(ops[3], smoothers[3], interpolators[2]);
	CHECK_THROWS_AS(builder.addFinestLevel(ops[0], smoothers[0], restrictors[0]),
	                RuntimeError);
}
TEST_CASE("CycleBuilder addIntermediateLevel throws exception if addFinestLevel isn't called",
          "[GMG::CycleBuilder]")
{
	std::array<std::shared_ptr<MockOperator>, 4> ops;
	std::array<MockSmoother, 4>                  smoothers;
	for (int i = 0; i < 4; i++) {
		ops[i] = make_shared<MockOperator>();
		smoothers[i].magic_number += i;
	}
	std::array<MockRestrictor, 3>   restrictors;
	std::array<MockInterpolator, 3> interpolators;

	GMG::CycleOpts       opts;
	GMG::CycleBuilder<2> builder(opts);
	CHECK_THROWS_AS(
	builder.addIntermediateLevel(ops[1], smoothers[1], restrictors[1], interpolators[0]),
	RuntimeError);
}
TEST_CASE("CycleBuilder addIntermediateLevel throws exception if called after addCoarsestLevel",
          "[GMG::CycleBuilder]")
{
	std::array<std::shared_ptr<MockOperator>, 4> ops;
	std::array<MockSmoother, 4>                  smoothers;
	for (int i = 0; i < 4; i++) {
		ops[i] = make_shared<MockOperator>();
		smoothers[i].magic_number += i;
	}
	std::array<MockRestrictor, 3>   restrictors;
	std::array<MockInterpolator, 3> interpolators;

	GMG::CycleOpts       opts;
	GMG::CycleBuilder<2> builder(opts);
	builder.addFinestLevel(ops[0], smoothers[0], restrictors[0]);
	builder.addIntermediateLevel(ops[1], smoothers[1], restrictors[1], interpolators[0]);
	builder.addIntermediateLevel(ops[2], smoothers[2], restrictors[2], interpolators[1]);
	builder.addCoarsestLevel(ops[3], smoothers[3], interpolators[2]);
	CHECK_THROWS_AS(
	builder.addIntermediateLevel(ops[2], smoothers[2], restrictors[2], interpolators[1]),
	RuntimeError);
}
TEST_CASE("CycleBuilder addCoarsestLevel throws exception if addFinestLevel isn't called",
          "[GMG::CycleBuilder]")
{
	std::array<std::shared_ptr<MockOperator>, 4> ops;
	std::array<MockSmoother, 4>                  smoothers;
	for (int i = 0; i < 4; i++) {
		ops[i] = make_shared<MockOperator>();
		smoothers[i].magic_number += i;
	}
	std::array<MockRestrictor, 3>   restrictors;
	std::array<MockInterpolator, 3> interpolators;

	GMG::CycleOpts       opts;
	GMG::CycleBuilder<2> builder(opts);
	CHECK_THROWS_AS(builder.addCoarsestLevel(ops[3], smoothers[3], interpolators[2]),
	                RuntimeError);
}
TEST_CASE("CycleBuilder addCoarsestLevel throws exception if called twice", "[GMG::CycleBuilder]")
{
	std::array<std::shared_ptr<MockOperator>, 4> ops;
	std::array<MockSmoother, 4>                  smoothers;
	for (int i = 0; i < 4; i++) {
		ops[i] = make_shared<MockOperator>();
		smoothers[i].magic_number += i;
	}
	std::array<MockRestrictor, 3>   restrictors;
	std::array<MockInterpolator, 3> interpolators;

	GMG::CycleOpts       opts;
	GMG::CycleBuilder<2> builder(opts);
	builder.addFinestLevel(ops[0], smoothers[0], restrictors[0]);
	builder.addIntermediateLevel(ops[1], smoothers[1], restrictors[1], interpolators[0]);
	builder.addIntermediateLevel(ops[2], smoothers[2], restrictors[2], interpolators[1]);
	builder.addCoarsestLevel(ops[3], smoothers[3], interpolators[2]);
	CHECK_THROWS_AS(builder.addCoarsestLevel(ops[3], smoothers[3], interpolators[2]),
	                RuntimeError);
}
TEST_CASE("CycleBuilder getCycle throws exception if addCoarsestLevel isn't called",
          "[GMG::CycleBuilder]")
{
	std::array<std::shared_ptr<MockOperator>, 4> ops;
	std::array<MockSmoother, 4>                  smoothers;
	for (int i = 0; i < 4; i++) {
		ops[i] = make_shared<MockOperator>();
		smoothers[i].magic_number += i;
	}
	std::array<MockRestrictor, 3>   restrictors;
	std::array<MockInterpolator, 3> interpolators;

	GMG::CycleOpts       opts;
	GMG::CycleBuilder<2> builder(opts);
	builder.addFinestLevel(ops[0], smoothers[0], restrictors[0]);
	builder.addIntermediateLevel(ops[1], smoothers[1], restrictors[1], interpolators[0]);
	builder.addIntermediateLevel(ops[2], smoothers[2], restrictors[2], interpolators[1]);

	CHECK_THROWS_AS(builder.getCycle(), RuntimeError);
}
TEST_CASE("CycleBuilder getCycle throws exception if no calls are made", "[GMG::CycleBuilder]")
{
	std::array<std::shared_ptr<MockOperator>, 4> ops;
	std::array<MockSmoother, 4>                  smoothers;
	for (int i = 0; i < 4; i++) {
		ops[i] = make_shared<MockOperator>();
		smoothers[i].magic_number += i;
	}
	std::array<std::shared_ptr<MockRestrictor>, 3>   restrictors;
	std::array<std::shared_ptr<MockInterpolator>, 3> interpolators;
	for (int i = 0; i < 3; i++) {
		restrictors[i]   = make_shared<MockRestrictor>();
		interpolators[i] = make_shared<MockInterpolator>();
	}

	GMG::CycleOpts       opts;
	GMG::CycleBuilder<2> builder(opts);

	CHECK_THROWS_AS(builder.getCycle(), RuntimeError);
}
