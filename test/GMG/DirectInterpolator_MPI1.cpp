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

#include "../utils/DomainReader.h"
#include <ThunderEgg/BiLinearGhostFiller.h>
#include <ThunderEgg/DomainTools.h>
#include <ThunderEgg/GMG/DirectInterpolator.h>

#include <doctest.h>

using namespace std;
using namespace ThunderEgg;

const string mesh_file = "mesh_inputs/2d_uniform_4x4_mpi1.json";
const string refined_mesh_file = "mesh_inputs/2d_uniform_2x2_refined_nw_mpi1.json";
TEST_CASE("Test DirectInterpolator on uniform 4x4")
{
  for (auto num_components : { 1, 2, 3 }) {
    for (auto nx : { 2, 10 }) {
      for (auto ny : { 2, 10 }) {
        int num_ghost = 1;
        DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
        Domain<2> d_fine = domain_reader.getFinerDomain();
        Domain<2> d_coarse = domain_reader.getCoarserDomain();

        Vector<2> coarse_vec(d_coarse, num_components);
        Vector<2> fine_vec(d_fine, num_components);
        Vector<2> fine_expected(d_fine, num_components);

        // set coarse vector
        for (auto pinfo : d_coarse.getPatchInfoVector()) {
          PatchView<double, 2> view = coarse_vec.getPatchView(pinfo.local_index);
          Loop::OverInteriorIndexes<3>(view, [&](const array<int, 3>& coord) { view[coord] = 1 + pinfo.id * nx * ny + coord[0] + coord[1] * nx + coord[2]; });
        }

        // set expected finer vector vector
        for (auto pinfo : d_fine.getPatchInfoVector()) {
          PatchView<double, 2> view = fine_expected.getPatchView(pinfo.local_index);

          Orthant<2> orth = pinfo.orth_on_parent;
          std::array<int, 2> starts;
          for (size_t i = 0; i < 2; i++) {
            starts[i] = orth.isOnSide(Side<2>(2 * i)) ? 0 : (view.getEnd()[i] + 1);
          }

          Loop::OverInteriorIndexes<3>(view, [&](const array<int, 3>& coord) { view[coord] = 1 + pinfo.parent_id * nx * ny + (coord[0] + starts[0]) / 2 + (coord[1] + starts[1]) / 2 * nx + coord[2]; });
        }

        GMG::DirectInterpolator<2> interpolator(d_coarse, d_fine);

        interpolator.interpolate(coarse_vec, fine_vec);

        for (auto pinfo : d_fine.getPatchInfoVector()) {
          PatchView<double, 2> vec_view = fine_vec.getPatchView(pinfo.local_index);
          PatchView<double, 2> expected_view = fine_expected.getPatchView(pinfo.local_index);
          Loop::OverInteriorIndexes<3>(vec_view, [&](const array<int, 3>& coord) { REQUIRE_EQ(vec_view[coord], doctest::Approx(expected_view[coord])); });
        }
      }
    }
  }
}
TEST_CASE("Linear Test DirectInterpolator with values already set on uniform 4x4")
{
  for (auto num_components : { 1, 2, 3 }) {
    for (auto nx : { 2, 10 }) {
      for (auto ny : { 2, 10 }) {
        int num_ghost = 1;
        DomainReader<2> domain_reader(mesh_file, { nx, ny }, num_ghost);
        Domain<2> d_fine = domain_reader.getFinerDomain();
        Domain<2> d_coarse = domain_reader.getCoarserDomain();

        Vector<2> coarse_vec(d_coarse, num_components);
        Vector<2> fine_vec(d_fine, num_components);
        Vector<2> fine_expected(d_fine, num_components);

        // set coarse vector
        for (auto pinfo : d_coarse.getPatchInfoVector()) {
          PatchView<double, 2> view = coarse_vec.getPatchView(pinfo.local_index);
          Loop::OverInteriorIndexes<3>(view, [&](const array<int, 3>& coord) { view[coord] = 1 + pinfo.id * nx * ny + coord[0] + coord[1] * nx + coord[2]; });
        }

        // set expected finer vector vector
        for (auto pinfo : d_fine.getPatchInfoVector()) {
          PatchView<double, 2> view = fine_expected.getPatchView(pinfo.local_index);

          Orthant<2> orth = pinfo.orth_on_parent;
          std::array<int, 2> starts;
          for (size_t i = 0; i < 2; i++) {
            starts[i] = orth.isOnSide(Side<2>(2 * i)) ? 0 : (view.getEnd()[i] + 1);
          }

          Loop::OverInteriorIndexes<3>(view, [&](const array<int, 3>& coord) { view[coord] = 2 + pinfo.parent_id * nx * ny + (coord[0] + starts[0]) / 2 + (coord[1] + starts[1]) / 2 * nx + coord[2]; });
        }

        fine_vec.set(1.0);

        GMG::DirectInterpolator<2> interpolator(d_coarse, d_fine);

        interpolator.interpolate(coarse_vec, fine_vec);

        for (auto pinfo : d_fine.getPatchInfoVector()) {
          PatchView<double, 2> vec_view = fine_vec.getPatchView(pinfo.local_index);
          PatchView<double, 2> expected_view = fine_expected.getPatchView(pinfo.local_index);
          Loop::OverInteriorIndexes<3>(vec_view, [&](const array<int, 3>& coord) { REQUIRE_EQ(vec_view[coord], doctest::Approx(expected_view[coord])); });
        }
      }
    }
  }
}
TEST_CASE("Test DirectInterpolator on refined 2x2")
{
  for (auto num_components : { 1, 2, 3 }) {
    for (auto nx : { 2, 10 }) {
      for (auto ny : { 2, 10 }) {
        int num_ghost = 1;
        DomainReader<2> domain_reader(refined_mesh_file, { nx, ny }, num_ghost);
        Domain<2> d_fine = domain_reader.getFinerDomain();
        Domain<2> d_coarse = domain_reader.getCoarserDomain();

        Vector<2> coarse_vec(d_coarse, num_components);
        Vector<2> fine_vec(d_fine, num_components);
        Vector<2> fine_expected(d_fine, num_components);

        // set coarse vector
        for (auto pinfo : d_coarse.getPatchInfoVector()) {
          PatchView<double, 2> view = coarse_vec.getPatchView(pinfo.local_index);
          Loop::OverInteriorIndexes<3>(view, [&](const array<int, 3>& coord) { view[coord] = 1 + pinfo.id * nx * ny + coord[0] + coord[1] * nx + coord[2]; });
        }

        // set expected finer vector vector
        for (auto pinfo : d_fine.getPatchInfoVector()) {
          PatchView<double, 2> view = fine_expected.getPatchView(pinfo.local_index);

          if (pinfo.hasCoarseParent()) {
            Orthant<2> orth = pinfo.orth_on_parent;
            std::array<int, 2> starts;
            for (size_t i = 0; i < 2; i++) {
              starts[i] = orth.isOnSide(Side<2>(2 * i)) ? 0 : (view.getEnd()[i] + 1);
            }

            Loop::OverInteriorIndexes<3>(view, [&](const array<int, 3>& coord) { view[coord] = 1 + pinfo.parent_id * nx * ny + (coord[0] + starts[0]) / 2 + (coord[1] + starts[1]) / 2 * nx + coord[2]; });
          } else {
            Loop::OverInteriorIndexes<3>(view, [&](const array<int, 3>& coord) { view[coord] = 1 + pinfo.id * nx * ny + coord[0] + coord[1] * nx + coord[2]; });
          }
        }

        GMG::DirectInterpolator<2> interpolator(d_coarse, d_fine);

        interpolator.interpolate(coarse_vec, fine_vec);

        for (auto pinfo : d_fine.getPatchInfoVector()) {
          PatchView<double, 2> vec_view = fine_vec.getPatchView(pinfo.local_index);
          PatchView<double, 2> expected_view = fine_expected.getPatchView(pinfo.local_index);
          Loop::OverInteriorIndexes<3>(vec_view, [&](const array<int, 3>& coord) { REQUIRE_EQ(vec_view[coord], doctest::Approx(expected_view[coord])); });
        }
      }
    }
  }
}
