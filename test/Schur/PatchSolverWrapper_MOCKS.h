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

#include <ThunderEgg/GhostFiller.h>
#include <ThunderEgg/PatchSolver.h>
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
class PatchFillingGhostFiller : public GhostFiller<D>
{
private:
  std::shared_ptr<bool> called = std::make_shared<bool>(false);
  double fill_value;

public:
  PatchFillingGhostFiller(double fill_value)
    : fill_value(fill_value)
  {
  }
  PatchFillingGhostFiller<D>* clone() const override { return new PatchFillingGhostFiller<D>(*this); }
  void fillGhost(const Vector<D>& u) const override
  {
    const_cast<Vector<D>&>(u).setWithGhost(fill_value);
    *called = true;
  }
  bool wasCalled() { return *called; }
};
template<int D>
class MockPatchSolver : public PatchSolver<D>
{
private:
  std::shared_ptr<std::set<int>> patch_ids_to_be_called;

public:
  MockPatchSolver(const Domain<D>& domain_in, const GhostFiller<D>& ghost_filler_in)
    : PatchSolver<D>(domain_in, ghost_filler_in)
  {
    patch_ids_to_be_called = std::make_shared<std::set<int>>();
    for (const PatchInfo<D>& pinfo : this->getDomain().getPatchInfoVector()) {
      patch_ids_to_be_called->insert(pinfo.id);
    }
  }
  MockPatchSolver<D>* clone() const override { return new MockPatchSolver<D>(*this); }
  void solveSinglePatch(const PatchInfo<D>& pinfo, const PatchView<const double, D>& f_view, const PatchView<double, D>& u_view) const override
  {
    CHECK(patch_ids_to_be_called->count(pinfo.id) == 1);
    patch_ids_to_be_called->erase(pinfo.id);
  }
  bool allPatchesCalled() { return patch_ids_to_be_called->empty(); }
};
template<int D>
class RHSGhostCheckingPatchSolver : public PatchSolver<D>
{
private:
  double schur_fill_value;
  std::shared_ptr<bool> was_called = std::make_shared<bool>(false);

public:
  RHSGhostCheckingPatchSolver(const Domain<D>& domain_in, const GhostFiller<D>& ghost_filler_in, double schur_fill_value)
    : PatchSolver<D>(domain_in, ghost_filler_in)
    , schur_fill_value(schur_fill_value)
  {
  }
  RHSGhostCheckingPatchSolver<D>* clone() const override { return new RHSGhostCheckingPatchSolver<D>(*this); }
  void solveSinglePatch(const PatchInfo<D>& pinfo, const PatchView<const double, D>& f_view, const PatchView<double, D>& u_view) const override
  {
    *was_called = true;
    for (Side<D> s : Side<D>::getValues()) {
      if (pinfo.hasNbr(s)) {
        auto ghosts = u_view.getSliceOn(s, { -1 });
        auto inner = u_view.getSliceOn(s, { 0 });
        Loop::OverInteriorIndexes<D>(ghosts, [&](const std::array<int, D>& coord) { CHECK((ghosts[coord] + inner[coord]) / 2 == doctest::Approx(schur_fill_value)); });
      }
    }
  }
  bool wasCalled() { return *was_called; }
};
} // namespace
} // namespace ThunderEgg