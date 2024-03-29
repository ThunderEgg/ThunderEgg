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

#ifndef MPIPATCHOPERATOR_MOCKS_H
#define MPIPATCHOPERATOR_MOCKS_H

#include <ThunderEgg/GhostFiller.h>
#include <ThunderEgg/PatchOperator.h>

#include <set>

#include <doctest.h>

namespace ThunderEgg {
namespace {
template<int D>
class MockGhostFiller : public GhostFiller<D>
{
private:
  std::shared_ptr<bool> called = std::make_shared<bool>(false);

public:
  MockGhostFiller<D>* clone() const override { return new MockGhostFiller<D>(*this); }
  void fillGhost(const Vector<D>& u) const override { *called = true; }
  bool wasCalled() { return *called; }
};
template<int D>
class MockPatchOperator : public PatchOperator<D>
{
private:
  std::shared_ptr<std::set<int>> patches_to_be_called;

public:
  MockPatchOperator(const Domain<D>& domain, const GhostFiller<D>& ghost_filler)
    : PatchOperator<D>(domain, ghost_filler)
  {
    patches_to_be_called.reset(new std::set<int>());
    for (const PatchInfo<D>& pinfo : domain.getPatchInfoVector()) {
      patches_to_be_called->insert(pinfo.id);
    }
  }
  MockPatchOperator<D>* clone() const override { return new MockPatchOperator<D>(*this); }
  void applySinglePatch(const PatchInfo<D>& pinfo, const PatchView<const double, D>& us, const PatchView<double, D>& fs) const override
  {
    CHECK(patches_to_be_called->count(pinfo.id) == 1);
    patches_to_be_called->erase(pinfo.id);
    INFO("LOCAL_INDEX: " << pinfo.local_index);
    std::array<int, D + 1> zero;
    zero.fill(0);
  }
  void applySinglePatchWithInternalBoundaryConditions(const PatchInfo<D>& pinfo, const PatchView<const double, D>& us, const PatchView<double, D>& fs) const override {}
  void modifyRHSForInternalBoundaryConditions(const PatchInfo<D>& pinfo, const PatchView<const double, D>& us, const PatchView<double, D>& fs) const override {}
  bool allPatchesCalled() { return patches_to_be_called->empty(); }
};
} // namespace
} // namespace ThunderEgg
#endif