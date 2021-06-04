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

#include <ThunderEgg/GhostFiller.h>
#include <ThunderEgg/PatchSolver.h>
#include <set>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

namespace ThunderEgg
{
namespace
{
template <int D>
class MockGhostFiller : public GhostFiller<D>
{
	private:
	mutable bool called = false;

	public:
	void fillGhost(std::shared_ptr<const Vector<D>> u) const override
	{
		called = true;
	}
	bool wasCalled()
	{
		return called;
	}
};
template <int D>
class PatchFillingGhostFiller : public GhostFiller<D>
{
	private:
	mutable bool called = false;
	double       fill_value;

	public:
	PatchFillingGhostFiller(double fill_value)
	: fill_value(fill_value) {}
	void fillGhost(std::shared_ptr<const Vector<D>> u) const override
	{
		std::const_pointer_cast<Vector<D>>(u)->setWithGhost(fill_value);
		called = true;
	}
	bool wasCalled()
	{
		return called;
	}
};
template <int D>
class MockPatchSolver : public PatchSolver<D>
{
	private:
	mutable std::set<int> patch_ids_to_be_called;

	public:
	MockPatchSolver(std::shared_ptr<const Domain<D>>      domain_in,
	                std::shared_ptr<const GhostFiller<D>> ghost_filler_in)
	: PatchSolver<D>(domain_in, ghost_filler_in)
	{
		{
			for (const PatchInfo<D> &pinfo : this->domain->getPatchInfoVector()) {
				patch_ids_to_be_called.insert(pinfo.id);
			}
		}
	}
	void solveSinglePatch(const PatchInfo<D> &                 pinfo,
	                      const std::vector<ComponentView<D>> &fs,
	                      std::vector<ComponentView<D>> &      us) const override
	{
		CHECK(patch_ids_to_be_called.count(pinfo.id) == 1);
		patch_ids_to_be_called.erase(pinfo.id);
	}
	bool allPatchesCalled()
	{
		return patch_ids_to_be_called.empty();
	}
};
template <int D>
class RHSGhostCheckingPatchSolver : public PatchSolver<D>
{
	private:
	double       schur_fill_value;
	mutable bool was_called = false;

	public:
	RHSGhostCheckingPatchSolver(std::shared_ptr<const Domain<D>>      domain_in,
	                            std::shared_ptr<const GhostFiller<D>> ghost_filler_in,
	                            double                                schur_fill_value)
	: PatchSolver<D>(domain_in, ghost_filler_in), schur_fill_value(schur_fill_value)
	{
	}
	void solveSinglePatch(const PatchInfo<D> &                 pinfo,
	                      const std::vector<ComponentView<D>> &fs,
	                      std::vector<ComponentView<D>> &      us) const override
	{
		was_called = true;
		for (Side<D> s : Side<D>::getValues()) {
			if (pinfo.hasNbr(s)) {
				auto ghosts = us[0].getSliceOn(s, {-1});
				auto inner  = us[0].getSliceOn(s, {0});
				nested_loop<D - 1>(
				ghosts.getStart(), ghosts.getEnd(), [&](const std::array<int, D - 1> &coord) {
					CHECK((ghosts[coord] + inner[coord]) / 2 == Catch::Approx(schur_fill_value));
				});
			}
		}
	}
	bool wasCalled()
	{
		return was_called;
	}
};
} // namespace
} // namespace ThunderEgg