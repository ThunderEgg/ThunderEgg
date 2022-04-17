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

#ifndef MPIGHOSTFILLER_MOCKS_H
#define MPIGHOSTFILLER_MOCKS_H

#include <ThunderEgg/DomainTools.h>
#include <ThunderEgg/MPIGhostFiller.h>

#include <list>
#include <tuple>

#include <doctest.h>

namespace ThunderEgg {
template<int D>
class CallMockMPIGhostFiller : public MPIGhostFiller<D>
{
private:
  int num_components;

  mutable std::list<std::tuple<const PatchInfo<D>*, Side<D>, const NbrType, const Orthant<D - 1>, int, int>> nbr_calls;

  mutable std::list<std::tuple<const PatchInfo<D>*, Edge, const NbrType, const Orthant<1>, int, int>> edge_nbr_calls;

  mutable std::list<std::tuple<const PatchInfo<D>*, Corner<D>, const NbrType, int, int>> corner_nbr_calls;

  mutable std::list<std::tuple<const PatchInfo<D>*, int>> local_calls;

public:
  void fillGhostCellsForNbrPatch(const PatchInfo<D>& pinfo, const PatchView<const double, D>& local_view, const PatchView<const double, D>& nbr_view, Side<D> side, NbrType nbr_type, Orthant<D - 1> orth_on_coarse) const override
  {
    INFO("Side: " << side);
    INFO("NbrType: " << nbr_type);
    called = true;
    bool remote_nbr = false;
    switch (nbr_type) {
      case NbrType::Normal:
        remote_nbr = pinfo.rank != pinfo.getNormalNbrInfo(side).rank;
        break;
      case NbrType::Fine:
        remote_nbr = pinfo.rank != pinfo.getFineNbrInfo(side).ranks[orth_on_coarse.getIndex()];
        break;
      case NbrType::Coarse:
        remote_nbr = pinfo.rank != pinfo.getCoarseNbrInfo(side).rank;
        break;
    }
    if (remote_nbr) {
      int num_ghost_cells = this->getDomain().getNumGhostCells();
      std::array<int, D + 1> expected_ghost_start = local_view.getStart();
      std::array<int, D + 1> expected_start = local_view.getStart();
      std::array<int, D + 1> expected_end = local_view.getEnd();
      std::array<int, D + 1> expected_ghost_end = local_view.getEnd();
      if (side.isHigherOnAxis()) {
        expected_ghost_start[side.getAxisIndex()] = -num_ghost_cells;
        expected_ghost_end[side.getAxisIndex()] = -1;
        expected_end[side.getAxisIndex()] = -1;
      } else {
        expected_start[side.getAxisIndex()] = expected_end[side.getAxisIndex()] + 1;
        expected_ghost_start[side.getAxisIndex()] = expected_end[side.getAxisIndex()] + 1;
        expected_ghost_end[side.getAxisIndex()] = expected_end[side.getAxisIndex()] + num_ghost_cells;
      }
      CHECK(nbr_view.getGhostStart() == expected_ghost_start);
      CHECK(nbr_view.getStart() == expected_start);
      CHECK(nbr_view.getEnd() == expected_end);
      CHECK(nbr_view.getGhostEnd() == expected_ghost_end);
    } else {
      CHECK(local_view.getGhostStart() == nbr_view.getGhostStart());
      CHECK(local_view.getStart() == nbr_view.getStart());
      CHECK(local_view.getEnd() == nbr_view.getEnd());
      CHECK(local_view.getGhostEnd() == nbr_view.getGhostEnd());
    }
    nbr_calls.emplace_back(&pinfo, side, nbr_type, orth_on_coarse, local_view.getEnd()[D] + 1, nbr_view.getEnd()[D] + 1);
  }

  void fillGhostCellsForEdgeNbrPatch(const PatchInfo<D>& pinfo, const PatchView<const double, D>& local_view, const PatchView<const double, D>& nbr_view, Edge edge, NbrType nbr_type, Orthant<1> orth_on_coarse) const override
  {
    if constexpr (D == 3) {
      INFO("Edge: " << edge);
      INFO("NbrType: " << nbr_type);
      called = true;
      bool remote_nbr = false;
      switch (nbr_type) {
        case NbrType::Normal:
          remote_nbr = pinfo.rank != pinfo.getNormalNbrInfo(edge).rank;
          break;
        case NbrType::Fine:
          remote_nbr = pinfo.rank != pinfo.getFineNbrInfo(edge).ranks[orth_on_coarse.getIndex()];
          break;
        case NbrType::Coarse:
          remote_nbr = pinfo.rank != pinfo.getCoarseNbrInfo(edge).rank;
          break;
      }
      if (remote_nbr) {
        int num_ghost_cells = this->getDomain().getNumGhostCells();
        std::array<int, D + 1> expected_ghost_start = local_view.getStart();
        std::array<int, D + 1> expected_start = local_view.getStart();
        std::array<int, D + 1> expected_end = local_view.getEnd();
        std::array<int, D + 1> expected_ghost_end = local_view.getEnd();
        for (Side<D> side : edge.getSides()) {
          if (side.isHigherOnAxis()) {
            expected_ghost_start[side.getAxisIndex()] = -num_ghost_cells;
            expected_ghost_end[side.getAxisIndex()] = -1;
            expected_end[side.getAxisIndex()] = -1;
          } else {
            expected_start[side.getAxisIndex()] = expected_end[side.getAxisIndex()] + 1;
            expected_ghost_start[side.getAxisIndex()] = expected_end[side.getAxisIndex()] + 1;
            expected_ghost_end[side.getAxisIndex()] = expected_end[side.getAxisIndex()] + num_ghost_cells;
          }
        }
        CHECK(nbr_view.getGhostStart() == expected_ghost_start);
        CHECK(nbr_view.getStart() == expected_start);
        CHECK(nbr_view.getEnd() == expected_end);
        CHECK(nbr_view.getGhostEnd() == expected_ghost_end);
      } else {
        CHECK(local_view.getGhostStart() == nbr_view.getGhostStart());
        CHECK(local_view.getStart() == nbr_view.getStart());
        CHECK(local_view.getEnd() == nbr_view.getEnd());
        CHECK(local_view.getGhostEnd() == nbr_view.getGhostEnd());
      }
      edge_nbr_calls.emplace_back(&pinfo, edge, nbr_type, orth_on_coarse, local_view.getEnd()[D] + 1, local_view.getEnd()[D] + 1);
    }
  }

  void fillGhostCellsForCornerNbrPatch(const PatchInfo<D>& pinfo, const PatchView<const double, D>& local_view, const PatchView<const double, D>& nbr_view, Corner<D> corner, NbrType nbr_type) const override
  {
    INFO("Corner: " << corner);
    INFO("NbrType: " << nbr_type);
    called = true;
    bool remote_nbr = false;
    switch (nbr_type) {
      case NbrType::Normal:
        remote_nbr = pinfo.rank != pinfo.getNormalNbrInfo(corner).rank;
        break;
      case NbrType::Fine:
        remote_nbr = pinfo.rank != pinfo.getFineNbrInfo(corner).ranks[0];
        break;
      case NbrType::Coarse:
        remote_nbr = pinfo.rank != pinfo.getCoarseNbrInfo(corner).rank;
        break;
    }
    if (remote_nbr) {
      int num_ghost_cells = this->getDomain().getNumGhostCells();
      std::array<int, D + 1> expected_ghost_start = local_view.getStart();
      std::array<int, D + 1> expected_start = local_view.getStart();
      std::array<int, D + 1> expected_end = local_view.getEnd();
      std::array<int, D + 1> expected_ghost_end = local_view.getEnd();
      for (Side<D> side : corner.getSides()) {
        if (side.isHigherOnAxis()) {
          expected_ghost_start[side.getAxisIndex()] = -num_ghost_cells;
          expected_ghost_end[side.getAxisIndex()] = -1;
          expected_end[side.getAxisIndex()] = -1;
        } else {
          expected_start[side.getAxisIndex()] = expected_end[side.getAxisIndex()] + 1;
          expected_ghost_start[side.getAxisIndex()] = expected_end[side.getAxisIndex()] + 1;
          expected_ghost_end[side.getAxisIndex()] = expected_end[side.getAxisIndex()] + num_ghost_cells;
        }
      }
      CHECK(nbr_view.getGhostStart() == expected_ghost_start);
      CHECK(nbr_view.getStart() == expected_start);
      CHECK(nbr_view.getEnd() == expected_end);
      CHECK(nbr_view.getGhostEnd() == expected_ghost_end);
    } else {
      CHECK(local_view.getGhostStart() == nbr_view.getGhostStart());
      CHECK(local_view.getStart() == nbr_view.getStart());
      CHECK(local_view.getEnd() == nbr_view.getEnd());
      CHECK(local_view.getGhostEnd() == nbr_view.getGhostEnd());
    }
    corner_nbr_calls.emplace_back(&pinfo, corner, nbr_type, local_view.getEnd()[D] + 1, nbr_view.getEnd()[D] + 1);
  }

  void fillGhostCellsForLocalPatch(const PatchInfo<D>& pinfo, const PatchView<const double, D>& local_view) const override
  {
    called = true;
    local_calls.emplace_back(&pinfo, local_view.getEnd()[D] + 1);
  }

  CallMockMPIGhostFiller(const Domain<D>& domain, int num_components, GhostFillingType fill_type)
    : MPIGhostFiller<D>(domain, fill_type)
    , num_components(num_components)
  {
  }
  CallMockMPIGhostFiller<D>* clone() const override { throw 3; }

  /**
   * @brief was any of the functions called
   */
  mutable bool called = false;

  void checkNbrCalls()
  {
    // remove from this collection the calls
    auto remaining_nbr_calls = nbr_calls;

    auto check_for_nbr_call = [&](const std::tuple<const PatchInfo<D>*, Side<D>, const NbrType, const Orthant<D - 1>, int, int>& call) {
      // check if call was made for neighbor
      auto found_call = std::find(remaining_nbr_calls.begin(), remaining_nbr_calls.end(), call);

      CHECK(found_call != remaining_nbr_calls.end());

      if (found_call != remaining_nbr_calls.end()) {
        remaining_nbr_calls.erase(found_call);
      }
    };
    for (const PatchInfo<D>& patch : this->getDomain().getPatchInfoVector()) {
      INFO("id: " << patch.id);
      std::string starts = "starts: ";
      for (size_t i = 0; i < D; i++) {
        starts = starts + " " + std::to_string(patch.starts[i]);
      }
      INFO(starts);

      for (auto side : Side<D>::getValues()) {
        INFO("side: " << side);
        if (patch.hasNbr(side)) {
          switch (patch.getNbrType(side)) {
            case NbrType::Normal: {
              INFO("NbrType: Normal");

              auto call = std::make_tuple(&patch, side, NbrType::Normal, Orthant<D - 1>::null(), num_components, num_components);
              check_for_nbr_call(call);

            } break;
            case NbrType::Fine: {
              INFO("NbrType: Fine");

              for (auto orthant : Orthant<D - 1>::getValues()) {
                INFO("Orthant: " << orthant.getIndex());
                auto call = std::make_tuple(&patch, side, NbrType::Fine, orthant, num_components, num_components);
                check_for_nbr_call(call);
              }
            } break;
            case NbrType::Coarse: {
              INFO("NbrType: Coarse");

              auto orthant = patch.getCoarseNbrInfo(side).orth_on_coarse;

              INFO("Orthant: " << orthant.getIndex());

              auto call = std::make_tuple(&patch, side, NbrType::Coarse, orthant, num_components, num_components);
              check_for_nbr_call(call);

            } break;
            default:
              CHECK(false);
          }
        }
      }
    }
    CHECK(remaining_nbr_calls.empty());
    INFO("REMAINING CALL");
    for (auto call : remaining_nbr_calls) {
      auto patch = *std::get<0>(call);
      INFO("id: " << patch.id);
      std::string starts = "starts: ";
      for (size_t i = 0; i < D; i++) {
        starts = starts + " " + std::to_string(patch.starts[i]);
      }
      INFO(starts);
      INFO("side: " << std::get<1>(call));
      std::string nbrtype = "";
      switch (std::get<2>(call)) {
        case NbrType::Normal:
          nbrtype = "normal";
          break;
        case NbrType::Fine:
          nbrtype = "fine";
          break;
        case NbrType::Coarse:
          nbrtype = "coarse";
          break;
        default:
          nbrtype = "???";
      }
      INFO("nbrtype: " << nbrtype);
      INFO("orthant: " << std::get<3>(call).getIndex());
      CHECK(false);
    }
  }
  void checkEdgeNbrCalls()
  {
    if constexpr (D == 3) {
      // remove from this collection the calls
      auto remaining_nbr_calls = edge_nbr_calls;

      auto check_for_nbr_call = [&](const std::tuple<const PatchInfo<D>*, Edge, const NbrType, const Orthant<1>, int, int>& call) {
        // check if call was made for neighbor
        auto found_call = std::find(remaining_nbr_calls.begin(), remaining_nbr_calls.end(), call);

        CHECK(found_call != remaining_nbr_calls.end());

        if (found_call != remaining_nbr_calls.end()) {
          remaining_nbr_calls.erase(found_call);
        }
      };
      for (const PatchInfo<D>& patch : this->getDomain().getPatchInfoVector()) {
        INFO("id: " << patch.id);
        std::string starts = "starts: ";
        for (size_t i = 0; i < D; i++) {
          starts = starts + " " + std::to_string(patch.starts[i]);
        }
        INFO(starts);

        for (auto edge : Edge::getValues()) {
          INFO("side: " << edge);
          if (patch.hasNbr(edge)) {
            switch (patch.getNbrType(edge)) {
              case NbrType::Normal: {
                INFO("NbrType: Normal");

                auto call = std::make_tuple(&patch, edge, NbrType::Normal, Orthant<1>::null(), num_components, num_components);
                check_for_nbr_call(call);

              } break;
              case NbrType::Fine: {
                INFO("NbrType: Fine");

                for (auto orthant : Orthant<1>::getValues()) {
                  INFO("Orthant: " << orthant.getIndex());
                  auto call = std::make_tuple(&patch, edge, NbrType::Fine, orthant, num_components, num_components);
                  check_for_nbr_call(call);
                }
              } break;
              case NbrType::Coarse: {
                INFO("NbrType: Coarse");

                auto orthant = patch.getCoarseNbrInfo(edge).orth_on_coarse;

                INFO("Orthant: " << orthant.getIndex());

                auto call = std::make_tuple(&patch, edge, NbrType::Coarse, orthant, num_components, num_components);
                check_for_nbr_call(call);

              } break;
              default:
                CHECK(false);
            }
          }
        }
      }
      CHECK(remaining_nbr_calls.empty());
      INFO("REMAINING CALL");
      for (auto call : remaining_nbr_calls) {
        auto patch = *std::get<0>(call);
        INFO("id: " << patch.id);
        std::string starts = "starts: ";
        for (size_t i = 0; i < D; i++) {
          starts = starts + " " + std::to_string(patch.starts[i]);
        }
        INFO(starts);
        INFO("side: " << std::get<1>(call));
        std::string nbrtype = "";
        switch (std::get<2>(call)) {
          case NbrType::Normal:
            nbrtype = "normal";
            break;
          case NbrType::Fine:
            nbrtype = "fine";
            break;
          case NbrType::Coarse:
            nbrtype = "coarse";
            break;
          default:
            nbrtype = "???";
        }
        INFO("nbrtype: " << nbrtype);
        INFO("orthant: " << std::get<3>(call).getIndex());
        CHECK(false);
      }
    }
  }
  void checkCornerNbrCalls()
  {
    // remove from this collection the calls
    auto remaining_nbr_calls = corner_nbr_calls;

    auto check_for_nbr_call = [&](const std::tuple<const PatchInfo<D>*, Corner<D>, const NbrType, int, int>& call) {
      // check if call was made for neighbor
      auto found_call = std::find(remaining_nbr_calls.begin(), remaining_nbr_calls.end(), call);

      CHECK(found_call != remaining_nbr_calls.end());

      if (found_call != remaining_nbr_calls.end()) {
        remaining_nbr_calls.erase(found_call);
      }
    };
    for (const PatchInfo<D>& patch : this->getDomain().getPatchInfoVector()) {
      INFO("id: " << patch.id);
      std::string starts = "starts: ";
      for (size_t i = 0; i < D; i++) {
        starts = starts + " " + std::to_string(patch.starts[i]);
      }
      INFO(starts);

      for (auto corner : Corner<D>::getValues()) {
        INFO("side: " << corner);
        if (patch.hasNbr(corner)) {
          switch (patch.getNbrType(corner)) {
            case NbrType::Normal: {
              INFO("NbrType: Normal");

              auto call = std::make_tuple(&patch, corner, NbrType::Normal, num_components, num_components);
              check_for_nbr_call(call);

            } break;
            case NbrType::Fine: {
              INFO("NbrType: Fine");

              auto call = std::make_tuple(&patch, corner, NbrType::Fine, num_components, num_components);
              check_for_nbr_call(call);
            } break;
            case NbrType::Coarse: {
              INFO("NbrType: Coarse");

              auto call = std::make_tuple(&patch, corner, NbrType::Coarse, num_components, num_components);
              check_for_nbr_call(call);

            } break;
            default:
              CHECK(false);
          }
        }
      }
    }
    CHECK(remaining_nbr_calls.empty());
    INFO("REMAINING CALL");
    for (auto call : remaining_nbr_calls) {
      auto patch = *std::get<0>(call);
      INFO("id: " << patch.id);
      std::string starts = "starts: ";
      for (size_t i = 0; i < D; i++) {
        starts = starts + " " + std::to_string(patch.starts[i]);
      }
      INFO(starts);
      INFO("side: " << std::get<1>(call));
      std::string nbrtype = "";
      switch (std::get<2>(call)) {
        case NbrType::Normal:
          nbrtype = "normal";
          break;
        case NbrType::Fine:
          nbrtype = "fine";
          break;
        case NbrType::Coarse:
          nbrtype = "coarse";
          break;
        default:
          nbrtype = "???";
      }
      INFO("nbrtype: " << nbrtype);
      CHECK(false);
    }
  }
  void checkLocalCalls()
  {
    // remove from this collection the calls
    auto remaining_local_calls = local_calls;

    for (const PatchInfo<D>& patch : this->getDomain().getPatchInfoVector()) {
      INFO("id: " << patch.id);
      std::string starts = "starts: ";
      for (size_t i = 0; i < D; i++) {
        starts = starts + " " + std::to_string(patch.starts[i]);
      }
      INFO(starts);

      // check for local call
      auto found_call = std::find(remaining_local_calls.begin(), remaining_local_calls.end(), std::make_tuple(&patch, num_components));

      CHECK(found_call != remaining_local_calls.end());

      if (found_call != remaining_local_calls.end()) {
        remaining_local_calls.erase(found_call);
      }
    }
    CHECK(remaining_local_calls.empty());
  }
  void checkCalls()
  {
    INFO("NumLocalPatches: " << this->getDomain().getNumLocalPatches());
    INFO("NumGlobalPatches: " << this->getDomain().getNumGlobalPatches());

    switch (this->getFillType()) {
      case GhostFillingType::Corners:
        checkCornerNbrCalls();
      case GhostFillingType::Edges:
        checkEdgeNbrCalls();
      case GhostFillingType::Faces:
        checkNbrCalls();
    }
    checkLocalCalls();
  }
};
template<int D>
class ExchangeMockMPIGhostFiller : public MPIGhostFiller<D>
{
public:
  void fillGhostCellsForNbrPatch(const PatchInfo<D>& pinfo, const PatchView<const double, D>& local_view, const PatchView<const double, D>& nbr_view, Side<D> side, NbrType nbr_type, Orthant<D - 1> orthant) const override
  {
    int index = 0;
    for (int i = 0; i < pinfo.num_ghost_cells; i++) {
      View<double, D> slice = nbr_view.getGhostSliceOn(side.opposite(), { (size_t)i });
      Loop::Nested<D>(slice.getStart(), slice.getEnd(), [&](const std::array<int, D>& coord) {
        INFO(i);
        INFO(coord[0]);
        INFO(index);
        slice[coord] += pinfo.id + index;
        index++;
      });
    }
  }

  void fillGhostCellsForEdgeNbrPatch(const PatchInfo<D>& pinfo, const PatchView<const double, D>& local_view, const PatchView<const double, D>& nbr_view, Edge edge, NbrType nbr_type, Orthant<1> orthant) const override
  {
    if constexpr (D == 3) {
      int index = 0;
      for (int j = 0; j < pinfo.num_ghost_cells; j++) {
        for (int i = 0; i < pinfo.num_ghost_cells; i++) {
          View<double, 2> slice = nbr_view.getGhostSliceOn(edge.opposite(), { (size_t)i, (size_t)j });
          Loop::Nested<2>(slice.getStart(), slice.getEnd(), [&](const std::array<int, 2>& coord) {
            slice[coord] += pinfo.id + index;
            index++;
          });
        }
      }
    }
  }

  void fillGhostCellsForCornerNbrPatch(const PatchInfo<D>& pinfo, const PatchView<const double, D>& local_view, const PatchView<const double, D>& nbr_view, Corner<D> corner, NbrType nbr_type) const override
  {
    int index = 0;
    std::array<size_t, D> start;
    start.fill(0);
    std::array<size_t, D> end;
    end.fill(pinfo.num_ghost_cells - 1);
    Loop::Nested<D>(start, end, [&](const std::array<size_t, D>& offset) {
      View<double, 1> ghosts = nbr_view.getGhostSliceOn(corner.opposite(), offset);
      for (int c = 0; c <= ghosts.getEnd()[0]; c++) {
        ghosts[{ c }] += pinfo.id + index;
        index++;
      }
    });
  }

  void fillGhostCellsForLocalPatch(const PatchInfo<D>& pinfo, const PatchView<const double, D>& local_data) const override {}

  ExchangeMockMPIGhostFiller(const Domain<D>& domain_in, GhostFillingType fill_type)
    : MPIGhostFiller<D>(domain_in, fill_type)
  {
  }
  ExchangeMockMPIGhostFiller<D>* clone() const override { return new ExchangeMockMPIGhostFiller<D>(*this); }

  void checkInterior(const Vector<D>& vec)
  {
    for (auto pinfo : this->getDomain().getPatchInfoVector()) {
      INFO("id: " << pinfo.id);
      std::string starts = "start: ";
      for (size_t i = 0; i < D; i++) {
        starts += " " + std::to_string(pinfo.starts[i]);
      }
      INFO(starts);
      std::string ns = "ns: ";
      for (size_t i = 0; i < D; i++) {
        ns += " " + std::to_string(pinfo.ns[i]);
      }
      INFO(ns);
      INFO("num_ghost_cells: " << pinfo.num_ghost_cells);
      for (int c = 0; c < vec.getNumComponents(); c++) {
        INFO("c: " << c);
        auto data = vec.getComponentView(c, pinfo.local_index);
        // check that vector was not modified on the interior
        Loop::Nested<D>(data.getStart(), data.getEnd(), [&](const std::array<int, D>& coord) {
          std::string coord_str = "coord: ";
          for (size_t i = 0; i < D; i++) {
            coord_str += " " + std::to_string(coord[i]);
          }
          INFO(coord_str);
          CHECK(data[coord] == pinfo.id);
        });
      }
    }
  }
  void checkCorners(const Vector<D>& vec)
  {
    for (auto pinfo : this->getDomain().getPatchInfoVector()) {
      INFO("id: " << pinfo.id);
      std::string starts = "start: ";
      for (size_t i = 0; i < D; i++) {
        starts += " " + std::to_string(pinfo.starts[i]);
      }
      INFO(starts);
      std::string ns = "ns: ";
      for (size_t i = 0; i < D; i++) {
        ns += " " + std::to_string(pinfo.ns[i]);
      }
      INFO(ns);
      INFO("num_ghost_cells: " << pinfo.num_ghost_cells);
      auto data = vec.getPatchView(pinfo.local_index);
      // check the ghost cells
      for (Corner<D> corner : Corner<D>::getValues()) {
        INFO("Corner: " << corner);

        if (pinfo.hasNbr(corner)) {
          switch (pinfo.getNbrType(corner)) {
            case NbrType::Normal: {
              INFO("NbrType: Normal");

              // value should be id of neighbor + index +c
              int index = 0;
              auto nbrinfo = pinfo.getNormalNbrInfo(corner);

              std::array<size_t, D> start;
              start.fill(0);
              std::array<size_t, D> end;
              end.fill(pinfo.num_ghost_cells - 1);
              Loop::Nested<D>(start, end, [&](std::array<size_t, D>& coord) {
                std::string coord_str = "coord: ";
                for (size_t i = 0; i < D; i++) {
                  coord_str += " " + std::to_string(coord[i]);
                }
                INFO(coord_str);
                for (int c = 0; c < vec.getNumComponents(); c++) {
                  CHECK(data.getGhostSliceOn(corner, coord)[{ c }] == nbrinfo.id + index);
                  index++;
                }
              });
            } break;
            case NbrType::Fine: {
              INFO("NbrType: Fine");

              auto nbrinfo = pinfo.getFineNbrInfo(corner);

              // value should be id of neighbors + index
              int index = 0;

              std::array<size_t, D> start;
              start.fill(0);
              std::array<size_t, D> end;
              end.fill(pinfo.num_ghost_cells - 1);
              Loop::Nested<D>(start, end, [&](std::array<size_t, D>& coord) {
                std::string coord_str = "coord: ";
                for (size_t i = 0; i < D; i++) {
                  coord_str += " " + std::to_string(coord[i]);
                }
                INFO(coord_str);
                for (int c = 0; c < vec.getNumComponents(); c++) {
                  CHECK(data.getGhostSliceOn(corner, coord)[{ c }] == nbrinfo.ids[0] + index);
                  index++;
                }
              });
            } break;
            case NbrType::Coarse: {
              INFO("NbrType: Coarse");

              // value should be id of neighbor + index
              int index = 0;
              auto nbrinfo = pinfo.getCoarseNbrInfo(corner);

              std::array<size_t, D> start;
              start.fill(0);
              std::array<size_t, D> end;
              end.fill(pinfo.num_ghost_cells - 1);
              Loop::Nested<D>(start, end, [&](std::array<size_t, D>& coord) {
                std::string coord_str = "coord: ";
                for (size_t i = 0; i < D; i++) {
                  coord_str += " " + std::to_string(coord[i]);
                }
                INFO(coord_str);
                for (int c = 0; c < vec.getNumComponents(); c++) {
                  CHECK(data.getGhostSliceOn(corner, coord)[{ c }] == nbrinfo.id + index);
                  index++;
                }
              });

            } break;
          }
        } else {
          // values should be zero on ghost cells
          INFO("Physical Boundary");

          std::array<int, D> start;
          for (int i = 0; i < D; i++) {
            start[i] = data.getGhostStart()[i];
          }
          std::array<int, D> end;
          end.fill(-1);
          Loop::Nested<D>(start, end, [&](std::array<int, D>& coord) {
            std::string coord_str = "coord: ";
            for (size_t i = 0; i < D; i++) {
              coord_str += " " + std::to_string(coord[i]);
            }
            INFO(coord_str);
            for (int c = 0; c < vec.getNumComponents(); c++) {
              CHECK(data.getSliceOn(corner, coord)[{ c }] == 0);
            }
          });
        }
      }
    }
  }

  void checkEdges(const Vector<D>& vec)
  {
    if constexpr (D == 3) {
      for (auto pinfo : this->getDomain().getPatchInfoVector()) {
        INFO("id: " << pinfo.id);
        std::string starts = "start: ";
        for (size_t i = 0; i < D; i++) {
          starts += " " + std::to_string(pinfo.starts[i]);
        }
        INFO(starts);
        std::string ns = "ns: ";
        for (size_t i = 0; i < D; i++) {
          ns += " " + std::to_string(pinfo.ns[i]);
        }
        INFO(ns);
        INFO("num_ghost_cells: " << pinfo.num_ghost_cells);
        auto data = vec.getPatchView(pinfo.local_index);
        // check the ghost cells
        for (Edge e : Edge::getValues()) {
          INFO("Edge: " << e);

          if (pinfo.hasNbr(e)) {
            switch (pinfo.getNbrType(e)) {
              case NbrType::Normal: {
                INFO("NbrType: Normal");

                // value should be id of neighbor + index + c
                int index = 0;
                auto nbrinfo = pinfo.getNormalNbrInfo(e);
                for (int j = 0; j < pinfo.num_ghost_cells; j++) {
                  for (int i = 0; i < pinfo.num_ghost_cells; i++) {
                    auto slice = data.getGhostSliceOn(e, { (size_t)i, (size_t)j });

                    Loop::Nested<2>(slice.getStart(), slice.getEnd(), [&](std::array<int, 2>& coord) {
                      std::string coord_str = "coord: ";
                      for (size_t i = 0; i < 2; i++) {
                        coord_str += " " + std::to_string(coord[i]);
                      }
                      INFO(coord_str);
                      CHECK(slice[coord] == nbrinfo.id + index);
                      index++;
                    });
                  }
                }

              } break;
              case NbrType::Fine: {
                INFO("NbrType: Fine");

                int ids = 0;
                auto nbrinfo = pinfo.getFineNbrInfo(e);
                for (size_t i = 0; i < nbrinfo.ids.size(); i++) {
                  ids += nbrinfo.ids[i];
                }

                // value should be id of neighbors + index
                int index = 0;
                for (int j = 0; j < pinfo.num_ghost_cells; j++) {
                  for (int i = 0; i < pinfo.num_ghost_cells; i++) {
                    auto slice = data.getGhostSliceOn(e, { (size_t)i, (size_t)j });
                    Loop::Nested<2>(slice.getStart(), slice.getEnd(), [&](std::array<int, 2>& coord) {
                      std::string coord_str = "coord: ";
                      for (size_t i = 0; i < 2; i++) {
                        coord_str += " " + std::to_string(coord[i]);
                      }
                      INFO(coord_str);
                      CHECK(slice[coord] == ids + 2 * (index));
                      index++;
                    });
                  }
                }
              } break;
              case NbrType::Coarse: {
                INFO("NbrType: Coarse");

                // value should be id of neighbor + index
                int index = 0;
                auto nbrinfo = pinfo.getCoarseNbrInfo(e);
                for (int j = 0; j < pinfo.num_ghost_cells; j++) {
                  for (int i = 0; i < pinfo.num_ghost_cells; i++) {
                    auto slice = data.getGhostSliceOn(e, { (size_t)i, (size_t)j });
                    Loop::Nested<2>(slice.getStart(), slice.getEnd(), [&](std::array<int, 2>& coord) {
                      std::string coord_str = "coord: ";
                      for (size_t i = 0; i < 2; i++) {
                        coord_str += " " + std::to_string(coord[i]);
                      }
                      INFO(coord_str);
                      CHECK(slice[coord] == (double)(nbrinfo.id + index));
                      index++;
                    });
                  }
                }

              } break;
            }
          } else {
            // values should be zero on ghost cells
            INFO("Physical Boundary");

            for (int j = 0; j < pinfo.num_ghost_cells; j++) {
              for (int i = 0; i < pinfo.num_ghost_cells; i++) {
                auto slice = data.getSliceOn(e, { -1 - i, -1 - j });
                Loop::Nested<2>(slice.getStart(), slice.getEnd(), [&](std::array<int, 2>& coord) {
                  INFO("coord: ");
                  for (size_t i = 0; i < 2; i++) {
                    INFO(coord[i]);
                  }
                  CHECK(slice[coord] == 0.0);
                });
              }
            }
          }
        }
      }
    }
  }

  void checkFaces(const Vector<D>& vec)
  {
    if constexpr (D == 3) {
      for (auto pinfo : this->getDomain().getPatchInfoVector()) {
        INFO("id: " << pinfo.id);
        std::string starts = "start: ";
        for (size_t i = 0; i < D; i++) {
          starts += " " + std::to_string(pinfo.starts[i]);
        }
        INFO(starts);
        std::string ns = "ns: ";
        for (size_t i = 0; i < D; i++) {
          ns += " " + std::to_string(pinfo.ns[i]);
        }
        INFO(ns);
        INFO("num_ghost_cells: " << pinfo.num_ghost_cells);
        auto data = vec.getPatchView(pinfo.local_index);
        // check the ghost cells
        for (Side<D> s : Side<D>::getValues()) {
          INFO("Side: " << s);

          if (pinfo.hasNbr(s)) {
            switch (pinfo.getNbrType(s)) {
              case NbrType::Normal: {
                INFO("NbrType: Normal");

                // value should be id of neighbor + index +c
                int index = 0;
                auto nbrinfo = pinfo.getNormalNbrInfo(s);
                for (int i = 0; i < pinfo.num_ghost_cells; i++) {
                  auto slice = data.getGhostSliceOn(s, { (size_t)i });

                  Loop::OverInteriorIndexes<D>(slice, [&](std::array<int, D>& coord) {
                    std::string coord_str = "coord: ";
                    for (size_t i = 0; i < D; i++) {
                      coord_str += " " + std::to_string(coord[i]);
                    }
                    INFO(coord_str);
                    CHECK(slice[coord] == nbrinfo.id + index);
                    index++;
                  });
                }

              } break;
              case NbrType::Fine: {
                INFO("NbrType: Fine");

                int ids = 0;
                auto nbrinfo = pinfo.getFineNbrInfo(s);
                for (size_t i = 0; i < nbrinfo.ids.size(); i++) {
                  ids += nbrinfo.ids[i];
                }

                // value should be id of neighbors + index
                int index = 0;
                for (int i = 0; i < pinfo.num_ghost_cells; i++) {
                  auto slice = data.getGhostSliceOn(s, { (size_t)i });

                  Loop::OverInteriorIndexes<D>(slice, [&](std::array<int, D>& coord) {
                    std::string coord_str = "coord: ";
                    for (size_t i = 0; i < D; i++) {
                      coord_str += " " + std::to_string(coord[i]);
                    }
                    INFO(coord_str);
                    CHECK(slice[coord] == ids + (1 << (D - 1)) * (index));
                    index++;
                  });
                }
              } break;
              case NbrType::Coarse: {
                INFO("NbrType: Coarse");

                // value should be id of neighbor + index
                int index = 0;
                auto nbrinfo = pinfo.getCoarseNbrInfo(s);
                for (int i = 0; i < pinfo.num_ghost_cells; i++) {
                  auto slice = data.getGhostSliceOn(s, { (size_t)i });

                  Loop::OverInteriorIndexes<D>(slice, [&](std::array<int, D>& coord) {
                    std::string coord_str = "coord: ";
                    for (size_t i = 0; i < D; i++) {
                      coord_str += " " + std::to_string(coord[i]);
                    }
                    INFO(coord_str);
                    CHECK(slice[coord] == nbrinfo.id + index);
                    index++;
                  });
                }

              } break;
            }
          } else {
            // values should be zero on ghost cells
            INFO("Physical Boundary");

            for (int i = 0; i < pinfo.num_ghost_cells; i++) {
              auto slice = data.getSliceOn(s, { -i - 1 });

              Loop::Nested<D>(slice.getStart(), slice.getEnd(), [&](std::array<int, D>& coord) {
                INFO("coord: ");
                for (size_t i = 0; i < D; i++) {
                  INFO(coord[i]);
                }
                CHECK(slice[coord] == 0.0);
              });
            }
          }
        }
      }
    }
  }

  void checkVector(const Vector<D>& vec)
  {
    INFO("NumLocalPatches: " << this->getDomain().getNumLocalPatches());
    INFO("NumGlobalPatches: " << this->getDomain().getNumGlobalPatches());
    checkInterior(vec);
    switch (this->getFillType()) {
      case GhostFillingType::Corners:
        checkCorners(vec);
      case GhostFillingType::Edges:
        checkEdges(vec);
      case GhostFillingType::Faces:
        checkFaces(vec);
    }
  }
};
} // namespace ThunderEgg
#endif