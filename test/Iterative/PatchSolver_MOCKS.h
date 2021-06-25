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
	void fillGhost(const Vector<D> &u) const override {}
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
	void applySinglePatch(const PatchInfo<D> &              pinfo,
	                      const PatchView<const double, D> &us, const PatchView<double, D> &fs) const override
	{
		num_apply_calls++;
	}
	void enforceBoundaryConditions(const PatchInfo<D> &pinfo, const PatchView<const double, D> &us) const override
	{
		bc_enforced = true;
	}
	void enforceZeroDirichletAtInternalBoundaries(const PatchInfo<D> &pinfo, const PatchView<const double, D> &us) const override
	{
		interior_dirichlet = true;
	}
	void modifyRHSForZeroDirichletAtInternalBoundaries(const PatchInfo<D> &              pinfo,
	                                                   const PatchView<const double, D> &us,
	                                                   const PatchView<double, D> &      fs) const override
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
	void applySinglePatch(const PatchInfo<D> &              pinfo,
	                      const PatchView<const double, D> &us, const PatchView<double, D> &fs) const override
	{
		loop_over_interior_indexes<D>(fs,
		                              [&](const std::array<int, D + 1> &coord) { fs[coord] += 1; });
	}
	void enforceBoundaryConditions(const PatchInfo<D> &pinfo, const PatchView<const double, D> &us) const override
	{
		bc_enforced = true;
	}
	void enforceZeroDirichletAtInternalBoundaries(const PatchInfo<D> &pinfo, const PatchView<const double, D> &us) const override
	{
		interior_dirichlet = true;
	}
	void modifyRHSForZeroDirichletAtInternalBoundaries(const PatchInfo<D> &              pinfo,
	                                                   const PatchView<const double, D> &us,
	                                                   const PatchView<double, D> &      fs) const override
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
	std::function<int(const Operator<D> &,
	                  Vector<D> &, const Vector<D> &,
	                  const Operator<D> *)>
	callback;

	public:
	MockSolver(
	std::function<int(const Operator<D> &,
	                  Vector<D> &, const Vector<D> &,
	                  const Operator<D> *)>
	callback)
	: callback(callback)
	{
	}
	int solve(const Operator<D> &A,
	          Vector<D> &x, const Vector<D> &b,
	          const Operator<D> *Mr = nullptr, bool output = false,
	          std::ostream &os = std::cout) const override
	{
		return callback(A, x, b, Mr);
	}
};

} // namespace
} // namespace ThunderEgg