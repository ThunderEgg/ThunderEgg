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
#include <ThunderEgg/Iterative/Solver.h>
#include <ThunderEgg/PatchOperator.h>
#include <set>

#include <doctest.h>

namespace ThunderEgg {
namespace {
template<int D>
class MockGhostFiller : public GhostFiller<D>
{
public:
  MockGhostFiller<D>* clone() const override { return new MockGhostFiller<D>(*this); }
  void fillGhost(const Vector<D>& u) const override {}
};
template<int D>
class MockPatchOperator : public PatchOperator<D>
{
private:
  std::shared_ptr<bool> rhs_was_modified = std::make_shared<bool>(false);
  std::shared_ptr<bool> bc_enforced = std::make_shared<bool>(false);
  std::shared_ptr<bool> interior_dirichlet = std::make_shared<bool>(false);
  std::shared_ptr<int> num_apply_calls = std::make_shared<int>(0);

public:
  MockPatchOperator(const Domain<D>& domain, const GhostFiller<D>& ghost_filler)
    : PatchOperator<D>(domain, ghost_filler)
  {
  }
  MockPatchOperator<D>* clone() const override { return new MockPatchOperator<D>(*this); }
  void applySinglePatch(const PatchInfo<D>& pinfo, const PatchView<const double, D>& us, const PatchView<double, D>& fs) const override
  {
    *bc_enforced = true;
    (*num_apply_calls)++;
  }
  void applySinglePatchWithInternalBoundaryConditions(const PatchInfo<D>& pinfo, const PatchView<const double, D>& us, const PatchView<double, D>& fs) const override
  {
    *bc_enforced = true;
    *interior_dirichlet = true;
    (*num_apply_calls)++;
  }
  void modifyRHSForInternalBoundaryConditions(const PatchInfo<D>& pinfo, const PatchView<const double, D>& us, const PatchView<double, D>& fs) const override { *rhs_was_modified = true; }
  bool rhsWasModified() { return *rhs_was_modified; }
  bool boundaryConditionsEnforced() { return *bc_enforced; }
  bool internalBoundaryConditionsEnforced() { return *interior_dirichlet; }
  int getNumApplyCalls() { return *num_apply_calls; }
};
template<int D>
class NonLinMockPatchOperator : public PatchOperator<D>
{
private:
  std::shared_ptr<bool> rhs_was_modified = std::make_shared<bool>(false);
  std::shared_ptr<bool> bc_enforced = std::make_shared<bool>(false);
  std::shared_ptr<bool> interior_dirichlet = std::make_shared<bool>(false);

public:
  NonLinMockPatchOperator(const Domain<D>& domain, const GhostFiller<D>& ghost_filler)
    : PatchOperator<D>(domain, ghost_filler)
  {
  }
  NonLinMockPatchOperator<D>* clone() const override { return new NonLinMockPatchOperator<D>(*this); }
  void applySinglePatch(const PatchInfo<D>& pinfo, const PatchView<const double, D>& us, const PatchView<double, D>& fs) const override
  {
    *bc_enforced = true;
    Loop::OverInteriorIndexes<D>(fs, [&](const std::array<int, D + 1>& coord) { fs[coord] += 1; });
  }
  void applySinglePatchWithInternalBoundaryConditions(const PatchInfo<D>& pinfo, const PatchView<const double, D>& us, const PatchView<double, D>& fs) const override
  {
    *bc_enforced = true;
    *interior_dirichlet = true;
    Loop::OverInteriorIndexes<D>(fs, [&](const std::array<int, D + 1>& coord) { fs[coord] += 1; });
  }
  void modifyRHSForInternalBoundaryConditions(const PatchInfo<D>& pinfo, const PatchView<const double, D>& us, const PatchView<double, D>& fs) const override { *rhs_was_modified = true; }
  bool rhsWasModified() { return *rhs_was_modified; }
  bool boundaryConditionsEnforced() { return *bc_enforced; }
  bool internalBoundaryConditionsEnforced() { return *interior_dirichlet; }
};
template<int D>
class MockSolver : public Iterative::Solver<D>
{
private:
  std::function<int(const Operator<D>&, Vector<D>&, const Vector<D>&, const Operator<D>*)> callback;

public:
  MockSolver<D>* clone() const override { return new MockSolver<D>(*this); }
  MockSolver(std::function<int(const Operator<D>&, Vector<D>&, const Vector<D>&, const Operator<D>*)> callback)
    : callback(callback)
  {
  }
  int solve(const Operator<D>& A, Vector<D>& x, const Vector<D>& b, const Operator<D>* Mr = nullptr, bool output = false, std::ostream& os = std::cout) const override { return callback(A, x, b, Mr); }
};

} // namespace
} // namespace ThunderEgg