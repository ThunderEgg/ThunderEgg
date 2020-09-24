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

#include "catch.hpp"

namespace ThunderEgg
{
namespace
{
template <int D> class MockGhostFiller : public GhostFiller<D>
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
template <int D> class PatchFillingGhostFiller : public GhostFiller<D>
{
	private:
	mutable bool called = false;
	double       fill_value;

	public:
	PatchFillingGhostFiller(double fill_value) : fill_value(fill_value) {}
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
template <int D> class MockPatchSolver : public PatchSolver<D>
{
	private:
	mutable std::set<std::shared_ptr<const PatchInfo<D>>> patches_to_be_called;

	public:
	MockPatchSolver(std::shared_ptr<const Domain<D>>      domain_in,
	                std::shared_ptr<const GhostFiller<D>> ghost_filler_in)
	: PatchSolver<D>(domain_in, ghost_filler_in)
	{
		{
			for (auto pinfo : this->domain->getPatchInfoVector()) {
				patches_to_be_called.insert(pinfo);
			}
		}
	}
	void solveSinglePatch(std::shared_ptr<const PatchInfo<D>> pinfo, LocalData<D> u,
	                      const LocalData<D> f) const override
	{
		CHECK(patches_to_be_called.count(pinfo) == 1);
		patches_to_be_called.erase(pinfo);
	}
	bool allPatchesCalled()
	{
		return patches_to_be_called.empty();
	}
};
template <int D> class RHSGhostCheckingPatchSolver : public PatchSolver<D>
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
	void solveSinglePatch(std::shared_ptr<const PatchInfo<D>> pinfo, LocalData<D> u,
	                      const LocalData<D> f) const override
	{
		was_called = true;
		for (Side<D> s : Side<D>::getValues()) {
			if (pinfo->hasNbr(s)) {
				auto ghosts = u.getGhostSliceOnSide(s, 1);
				auto inner  = u.getSliceOnSide(s);
				nested_loop<D - 1>(
				ghosts.getStart(), ghosts.getEnd(), [&](const std::array<int, D - 1> &coord) {
					CHECK((ghosts[coord] + inner[coord]) / 2 == Approx(schur_fill_value));
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