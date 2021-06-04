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
#include <ThunderEgg/Iterative/Solver.h>
#include <ThunderEgg/PatchOperator.h>
#include <set>

#include <catch2/catch_test_macros.hpp>

namespace ThunderEgg
{
namespace
{
template <int D>
class MockGhostFiller : public GhostFiller<D>
{
	public:
	void fillGhost(std::shared_ptr<const Vector<D>> u) const override {}
};
template <int D>
class MockPatchOperator : public PatchOperator<D>
{
	private:
	mutable bool rhs_was_modified   = false;
	mutable bool bc_enforced        = false;
	mutable bool interior_dirichlet = false;
	mutable int  num_apply_calls    = 0;

	public:
	MockPatchOperator(std::shared_ptr<const Domain<D>>      domain,
	                  std::shared_ptr<const GhostFiller<D>> ghost_filler)
	: PatchOperator<D>(domain, ghost_filler)
	{
	}
	void applySinglePatch(const PatchInfo<D> &        pinfo,
	                      const std::vector<View<D>> &us, std::vector<View<D>> &fs) const override
	{
		num_apply_calls++;
	}
	void enforceBoundaryConditions(const PatchInfo<D> &pinfo, const std::vector<View<D>> &us) const override
	{
		bc_enforced = true;
	}
	void enforceZeroDirichletAtInternalBoundaries(const PatchInfo<D> &pinfo, const std::vector<View<D>> &us) const override
	{
		interior_dirichlet = true;
	}
	void modifyRHSForZeroDirichletAtInternalBoundaries(const PatchInfo<D> &        pinfo,
	                                                   const std::vector<View<D>> &us,
	                                                   std::vector<View<D>> &      fs) const override
	{
		rhs_was_modified = true;
	}
	bool rhsWasModified()
	{
		return rhs_was_modified;
	}
	bool boundaryConditionsEnforced()
	{
		return bc_enforced;
	}
	bool internalBoundaryConditionsEnforced()
	{
		return interior_dirichlet;
	}
	int getNumApplyCalls()
	{
		return num_apply_calls;
	}
};
template <int D>
class NonLinMockPatchOperator : public PatchOperator<D>
{
	private:
	mutable bool rhs_was_modified   = false;
	mutable bool bc_enforced        = false;
	mutable bool interior_dirichlet = false;

	public:
	NonLinMockPatchOperator(std::shared_ptr<const Domain<D>>      domain,
	                        std::shared_ptr<const GhostFiller<D>> ghost_filler)
	: PatchOperator<D>(domain, ghost_filler)
	{
	}
	void applySinglePatch(const PatchInfo<D> &        pinfo,
	                      const std::vector<View<D>> &us, std::vector<View<D>> &fs) const override
	{
		for (size_t c = 0; c < fs.size(); c++) {
			nested_loop<D>(fs[c].getStart(), fs[c].getEnd(),
			               [&](const std::array<int, D> &coord) { fs[c][coord] += 1; });
		}
	}
	void enforceBoundaryConditions(const PatchInfo<D> &pinfo, const std::vector<View<D>> &us) const override
	{
		bc_enforced = true;
	}
	void enforceZeroDirichletAtInternalBoundaries(const PatchInfo<D> &pinfo, const std::vector<View<D>> &us) const override
	{
		interior_dirichlet = true;
	}
	void modifyRHSForZeroDirichletAtInternalBoundaries(const PatchInfo<D> &        pinfo,
	                                                   const std::vector<View<D>> &us,
	                                                   std::vector<View<D>> &      fs) const override
	{
		rhs_was_modified = true;
	}
	bool rhsWasModified()
	{
		return rhs_was_modified;
	}
	bool boundaryConditionsEnforced()
	{
		return bc_enforced;
	}
	bool internalBoundaryConditionsEnforced()
	{
		return interior_dirichlet;
	}
};
template <int D>
class MockSolver : public Iterative::Solver<D>
{
	private:
	std::function<int(std::shared_ptr<VectorGenerator<D>>, std::shared_ptr<const Operator<D>>,
	                  std::shared_ptr<Vector<D>>, std::shared_ptr<const Vector<D>>,
	                  std::shared_ptr<const Operator<D>>)>
	callback;

	public:
	MockSolver(
	std::function<int(std::shared_ptr<VectorGenerator<D>>, std::shared_ptr<const Operator<D>>,
	                  std::shared_ptr<Vector<D>>, std::shared_ptr<const Vector<D>>,
	                  std::shared_ptr<const Operator<D>>)>
	callback)
	: callback(callback)
	{
	}
	int solve(std::shared_ptr<VectorGenerator<D>> vg, std::shared_ptr<const Operator<D>> A,
	          std::shared_ptr<Vector<D>> x, std::shared_ptr<const Vector<D>> b,
	          std::shared_ptr<const Operator<D>> Mr = nullptr, bool output = false,
	          std::ostream &os = std::cout) const override
	{
		return callback(vg, A, x, b, Mr);
	}
};

} // namespace
} // namespace ThunderEgg