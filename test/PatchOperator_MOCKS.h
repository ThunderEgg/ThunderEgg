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

#ifndef MPIPATCHOPERATOR_MOCKS_H
#define MPIPATCHOPERATOR_MOCKS_H

#include <ThunderEgg/GhostFiller.h>
#include <ThunderEgg/PatchOperator.h>
#include <set>

#include "catch.hpp"

namespace ThunderEgg
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
template <int D> class MockPatchOperator : public PatchOperator<D>
{
	private:
	std::shared_ptr<Vector<D>>                            u_vec;
	std::shared_ptr<Vector<D>>                            f_vec;
	mutable std::set<std::shared_ptr<const PatchInfo<D>>> patches_to_be_called;

	public:
	MockPatchOperator(std::shared_ptr<const Domain<D>>      domain_in,
	                  std::shared_ptr<const GhostFiller<D>> ghost_filler_in,
	                  std::shared_ptr<Vector<D>> u_in, std::shared_ptr<Vector<D>> f_in)
	: PatchOperator<D>(domain_in, ghost_filler_in), u_vec(u_in), f_vec(f_in)
	{
		for (auto pinfo : this->domain->getPatchInfoVector()) {
			patches_to_be_called.insert(pinfo);
		}
	}
	void applySinglePatch(std::shared_ptr<const PatchInfo<D>> pinfo, LocalData<D> u,
	                      LocalData<D> f) const override
	{
		CHECK(patches_to_be_called.count(pinfo) == 1);
		patches_to_be_called.erase(pinfo);
		CHECK(u_vec->getLocalData(pinfo->local_index).getPtr() == u.getPtr());
		CHECK(f_vec->getLocalData(pinfo->local_index).getPtr() == f.getPtr());
	}
	void addGhostToRHS(std::shared_ptr<const PatchInfo<D>> pinfo, LocalData<D> u,
	                   LocalData<D> f) const override
	{
	}
	bool allPatchesCalled()
	{
		return patches_to_be_called.empty();
	}
};
} // namespace ThunderEgg
#endif