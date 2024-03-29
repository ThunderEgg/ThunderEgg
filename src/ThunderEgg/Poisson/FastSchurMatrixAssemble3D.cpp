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

#include "FastSchurMatrixAssemble3D.h"
#include <ThunderEgg/MPIGhostFiller.h>
#include <ThunderEgg/TriLinearGhostFiller.h>
using namespace std;
using namespace ThunderEgg;
using namespace ThunderEgg::Schur;
namespace {
enum class Rotation : char
{
  x_cw,
  x_ccw,
  y_cw,
  y_ccw,
  z_cw,
  z_ccw
};
bool
sideIsLeftOriented(const Side<3> s)
{
  return (s == Side<3>::north() || s == Side<3>::west() || s == Side<3>::bottom());
}
class Block
{
public:
  static const Side<3> side_table[6][6];
  static const char rots_table[6][6];
  static const vector<Rotation> main_rot_plan[6];
  static const vector<Rotation> aux_rot_plan_dirichlet[6];
  static const vector<Rotation> aux_rot_plan_neumann[16];
  static const char rot_quad_lookup_left[4][4];
  static const char rot_quad_lookup_right[4][4];
  static const char quad_flip_lookup[4];
  IfaceType<3> type;
  Side<3> main;
  Side<3> aux;
  int j;
  int i;
  bitset<6> non_dirichlet_boundary;
  bool orig_main_is_left_oriented;
  bool orig_aux_is_left_oriented;
  unsigned char main_rotation = 0;
  unsigned char aux_rotation = 0;
  Block(Side<3> main,
        int j,
        Side<3> aux,
        int i,
        bitset<6> non_dirichlet_boundary,
        IfaceType<3> type)
    : type(type)
    , main(main)
    , aux(aux)
    , j(j)
    , i(i)
    , non_dirichlet_boundary(non_dirichlet_boundary)
    , orig_main_is_left_oriented(sideIsLeftOriented(main))
    , orig_aux_is_left_oriented(sideIsLeftOriented(aux))
  {
    rotate();
  }
  void applyRotation(const Rotation rot)
  {
    // main rotation
    main_rotation = (main_rotation + rots_table[static_cast<int>(rot)][main.getIndex()]) % 4;
    // aux rotation
    aux_rotation = (aux_rotation + rots_table[static_cast<int>(rot)][aux.getIndex()]) % 4;
    main = side_table[static_cast<int>(rot)][main.getIndex()];
    aux = side_table[static_cast<int>(rot)][aux.getIndex()];
    bitset<6> old_neumann = non_dirichlet_boundary;
    for (int idx = 0; idx < 6; idx++) {
      non_dirichlet_boundary[side_table[static_cast<int>(rot)][idx].getIndex()] = old_neumann[idx];
    }
  }
  void rotate()
  {
    for (Rotation rot : main_rot_plan[main.getIndex()]) {
      applyRotation(rot);
    }
    if (non_dirichlet_boundary.to_ulong() == 0) {
      for (Rotation rot : aux_rot_plan_dirichlet[aux.getIndex()]) {
        applyRotation(rot);
      }
    } else {
      for (Rotation rot : aux_rot_plan_neumann[non_dirichlet_boundary.to_ulong() >> 2]) {
        applyRotation(rot);
      }
    }
    // updated iface type
    auto rotateQuad = [&](int quad) {
      if (orig_aux_is_left_oriented) {
        quad = rot_quad_lookup_left[aux_rotation][quad];
      } else {
        quad = rot_quad_lookup_right[aux_rotation][quad];
      }
      if (auxFlipped()) {
        quad = quad_flip_lookup[quad];
      }
      return quad;
    };
    if (type.isFineToCoarse() || type.isFineToFine() || type.isCoarseToFine()) {
      int quad = (int)type.getOrthant().getIndex();
      quad = rotateQuad(quad);
      type.setOrthant(Orthant<2>((unsigned char)quad));
    }
  }
  bool operator<(const Block& b) const
  {
    return std::tie(i, j, main_rotation, orig_main_is_left_oriented, aux_rotation) <
           std::tie(b.i, b.j, b.main_rotation, b.orig_main_is_left_oriented, b.aux_rotation);
  }
  bool mainFlipped() const { return sideIsLeftOriented(main) != orig_main_is_left_oriented; }
  bool auxFlipped() const { return sideIsLeftOriented(aux) != orig_aux_is_left_oriented; }
};

const Side<3> Block::side_table[6][6] = { { Side<3>::west(),
                                            Side<3>::east(),
                                            Side<3>::top(),
                                            Side<3>::bottom(),
                                            Side<3>::south(),
                                            Side<3>::north() },
                                          { Side<3>::west(),
                                            Side<3>::east(),
                                            Side<3>::bottom(),
                                            Side<3>::top(),
                                            Side<3>::north(),
                                            Side<3>::south() },
                                          { Side<3>::bottom(),
                                            Side<3>::top(),
                                            Side<3>::south(),
                                            Side<3>::north(),
                                            Side<3>::east(),
                                            Side<3>::west() },
                                          { Side<3>::top(),
                                            Side<3>::bottom(),
                                            Side<3>::south(),
                                            Side<3>::north(),
                                            Side<3>::west(),
                                            Side<3>::east() },
                                          { Side<3>::north(),
                                            Side<3>::south(),
                                            Side<3>::west(),
                                            Side<3>::east(),
                                            Side<3>::bottom(),
                                            Side<3>::top() },
                                          { Side<3>::south(),
                                            Side<3>::north(),
                                            Side<3>::east(),
                                            Side<3>::west(),
                                            Side<3>::bottom(),
                                            Side<3>::top() } };
const char Block::rots_table[6][6] = { { 3, 1, 0, 0, 2, 2 }, { 1, 3, 2, 2, 0, 0 },
                                       { 1, 3, 3, 1, 1, 3 }, { 1, 3, 1, 3, 3, 1 },
                                       { 0, 0, 0, 0, 3, 1 }, { 0, 0, 0, 0, 1, 3 } };
const vector<Rotation> Block::main_rot_plan[6] = { {},
                                                   { Rotation::z_cw, Rotation::z_cw },
                                                   { Rotation::z_cw },
                                                   { Rotation::z_ccw },
                                                   { Rotation::y_ccw },
                                                   { Rotation::y_cw } };
const vector<Rotation> Block::aux_rot_plan_dirichlet[6] = {
  {}, {}, {}, { Rotation::x_cw, Rotation::x_cw }, { Rotation::x_cw }, { Rotation::x_ccw }
};
const vector<Rotation> Block::aux_rot_plan_neumann[16] = { {},
                                                           {},
                                                           { Rotation::x_cw, Rotation::x_cw },
                                                           {},
                                                           { Rotation::x_cw },
                                                           { Rotation::x_cw },
                                                           { Rotation::x_cw, Rotation::x_cw },
                                                           { Rotation::x_cw, Rotation::x_cw },
                                                           { Rotation::x_ccw },
                                                           {},
                                                           { Rotation::x_ccw },
                                                           {},
                                                           { Rotation::x_ccw },
                                                           { Rotation::x_cw },
                                                           { Rotation::x_ccw },
                                                           {} };
const char Block::rot_quad_lookup_left[4][4] = { { 0, 1, 2, 3 },
                                                 { 1, 3, 0, 2 },
                                                 { 3, 2, 1, 0 },
                                                 { 2, 0, 3, 1 } };
const char Block::rot_quad_lookup_right[4][4] = { { 0, 1, 2, 3 },
                                                  { 2, 0, 3, 1 },
                                                  { 3, 2, 1, 0 },
                                                  { 1, 3, 0, 2 } };
const char Block::quad_flip_lookup[4] = { 1, 0, 3, 2 };

/**
 * @brief Get the view for the buffer
 *
 * @param buffer_ptr pointer to the ghost cells position in the buffer
 * @param pinfo  the PatchInfo object
 * @param side  the side that the ghost cells are on
 * @return PatchView<D> the view
 */
PatchView<const double, 3>
getPatchViewForBuffer(double* buffer_ptr, const PatchInfo<3>& pinfo, const Side<3> side)
{
  auto ns = pinfo.ns;
  int num_ghost_cells = pinfo.num_ghost_cells;
  // determine striding
  std::array<int, 4> strides;
  strides[0] = 1;
  for (size_t i = 1; i < 3; i++) {
    if (i == side.getAxisIndex() + 1) {
      strides[i] = num_ghost_cells * strides[i - 1];
    } else {
      strides[i] = ns[i - 1] * strides[i - 1];
    }
  }
  strides[3] = ns[2] * strides[2];
  // determine bounds of buffer
  std::array<int, 4> ghost_start;
  ghost_start.fill(0);

  std::array<int, 4> start;
  start.fill(0);

  std::array<int, 4> end;
  for (int i = 0; i < 3; i++) {
    end[i] = ns[i] - 1;
  }
  end[3] = 0;

  std::array<int, 4> ghost_end = end;

  if (side.isLowerOnAxis()) {
    ghost_start[side.getAxisIndex()] = -1;
    ghost_end[side.getAxisIndex()] = -1;
    end[side.getAxisIndex()] = -1;
  } else {
    ghost_start[side.getAxisIndex()] = ns[side.getAxisIndex()];
    ghost_end[side.getAxisIndex()] = ns[side.getAxisIndex()];
    start[side.getAxisIndex()] = ns[side.getAxisIndex()];
  }

  return PatchView<double, 3>(buffer_ptr, strides, ghost_start, start, end, ghost_end);
}

/**
 * @brief Fill a block column for a normal interface
 *
 * @param j the column of the block to fill
 * @param u the patch data
 * @param s the side of the patch that the block is on
 * @param block
 */
void
FillBlockColumnForNormalInterface(int j,
                                  const PatchView<const double, 3>& u,
                                  Side<3> s,
                                  std::vector<double>& block)
{
  int n = u.getEnd()[0] + 1;
  View<const double, 3> slice = u.getSliceOn(s, { 0 });
  for (int yi = 0; yi < n; yi++) {
    for (int xi = 0; xi < n; xi++) {
      block[(xi + n * yi) * n * n + j] = -slice(xi, yi, 0) / 2;
    }
  }
}
/**
 * @brief Fill a block column for a coarse to coarse interface
 *
 * @param j the column of the block to fill
 * @param u the patch data
 * @param s the side of the patch that the block is on
 * @param ghost_filler the GhostFiller
 * @param pinfo the PatchInfo
 * @param block
 */
void
FillBlockColumnForCoarseToCoarseInterface(int j,
                                          const PatchView<const double, 3>& u,
                                          Side<3> s,
                                          const MPIGhostFiller<3>& ghost_filler,
                                          const PatchInfo<3>& pinfo,
                                          std::vector<double>& block)
{
  int n = pinfo.ns[0];
  PatchInfo<3> new_pinfo = pinfo;
  new_pinfo.setNbrInfo(Side<3>::west(), nullptr);
  new_pinfo.setNbrInfo(s, new FineNbrInfo<2>());
  ghost_filler.fillGhostCellsForLocalPatch(new_pinfo, u);
  View<const double, 3> slice = u.getSliceOn(s, { 0 });
  View<double, 3> ghost_slice = u.getGhostSliceOn(s, { 0 });
  for (int yi = 0; yi < n; yi++) {
    for (int xi = 0; xi < n; xi++) {
      block[(xi + yi * n) * n * n + j] = -(slice(xi, yi, 0) + ghost_slice(xi, yi, 0)) / 2;
      ghost_slice(xi, yi, 0) = 0;
    }
  }
}
/**
 * @brief Fill a block column for a fine to fine interface
 *
 * @param j the column of the block to fill
 * @param u the patch data
 * @param s the side of the patch that the block is on
 * @param ghost_filler the GhostFiller
 * @param pinfo the PatchInfo
 * @param type the IfaceType
 * @param block
 */
void
FillBlockColumnForFineToFineInterface(int j,
                                      const PatchView<const double, 3>& u,
                                      Side<3> s,
                                      const MPIGhostFiller<3>& ghost_filler,
                                      const PatchInfo<3>& pinfo,
                                      IfaceType<3> type,
                                      std::vector<double>& block)
{
  int n = pinfo.ns[0];
  PatchInfo<3> new_pinfo = pinfo;
  new_pinfo.setNbrInfo(Side<3>::west(), nullptr);
  new_pinfo.setNbrInfo(s, new CoarseNbrInfo<2>(100, type.getOrthant()));
  ghost_filler.fillGhostCellsForLocalPatch(new_pinfo, u);
  View<const double, 3> slice = u.getSliceOn(s, { 0 });
  View<double, 3> ghost_slice = u.getGhostSliceOn(s, { 0 });
  for (int yi = 0; yi < n; yi++) {
    for (int xi = 0; xi < n; xi++) {
      block[(xi + yi * n) * n * n + j] = -(slice(xi, yi, 0) + ghost_slice(xi, yi, 0)) / 2;
      ghost_slice(xi, yi, 0) = 0;
    }
  }
}
/**
 * @brief Fill a block column for a coarse to fine interface
 *
 * @param j the column of the block to fill
 * @param u the patch data
 * @param s the side of the patch that the block is on
 * @param ghost_filler the GhostFiller
 * @param pinfo the PatchInfo
 * @param type the IfaceType
 * @param block
 */
void
FillBlockColumnForCoarseToFineInterface(int j,
                                        const PatchView<const double, 3>& u,
                                        Side<3> s,
                                        const MPIGhostFiller<3>& ghost_filler,
                                        const PatchInfo<3>& pinfo,
                                        IfaceType<3> type,
                                        std::vector<double>& block)
{
  int n = pinfo.ns[0];
  PatchInfo<3> new_pinfo = pinfo;
  new_pinfo.setNbrInfo(Side<3>::west(), nullptr);
  new_pinfo.setNbrInfo(s, new FineNbrInfo<2>());
  vector<double> ghosts(n * n);
  PatchView<const double, 3> nbr_view = getPatchViewForBuffer(ghosts.data(), pinfo, s.opposite());
  ghost_filler.fillGhostCellsForNbrPatch(
    new_pinfo, u, nbr_view, s, NbrType::Fine, type.getOrthant());
  for (int yi = 0; yi < n; yi++) {
    for (int xi = 0; xi < n; xi++) {
      block[(xi + yi * n) * n * n + j] = -ghosts[xi + yi * n] / 2;
    }
  }
}
/**
 * @brief Fill a block column for a fine to coarse interface
 *
 * @param j the column of the block to fill
 * @param u the patch data
 * @param s the side of the patch that the block is on
 * @param ghost_filler the GhostFiller
 * @param pinfo the PatchInfo
 * @param type the IfaceType
 * @param block
 */
void
FillBlockColumnForFineToCoarseInterface(int j,
                                        const PatchView<const double, 3>& u,
                                        Side<3> s,
                                        const MPIGhostFiller<3>& ghost_filler,
                                        const PatchInfo<3>& pinfo,
                                        IfaceType<3> type,
                                        std::vector<double>& block)
{
  int n = pinfo.ns[0];
  PatchInfo<3> new_pinfo = pinfo;
  new_pinfo.setNbrInfo(Side<3>::west(), nullptr);
  new_pinfo.setNbrInfo(s, new CoarseNbrInfo<2>(100, type.getOrthant()));
  vector<double> ghosts(n * n);
  PatchView<const double, 3> nbr_view = getPatchViewForBuffer(ghosts.data(), pinfo, s.opposite());
  ghost_filler.fillGhostCellsForNbrPatch(
    new_pinfo, u, nbr_view, s, NbrType::Coarse, type.getOrthant());
  for (int yi = 0; yi < n; yi++) {
    for (int xi = 0; xi < n; xi++) {
      block[(xi + yi * n) * n * n + j] = -ghosts[xi + yi * n] / 2;
    }
  }
}
/**
 * @brief Get a vector of set<Block>, each set of blocks have the same boundary conditions
 *
 * @param iface_domain the InterfaceDomain
 * @return vector<set<Block>> the vector of blocks
 */
vector<set<Block>>
GetBlocks(const InterfaceDomain<3>& iface_domain, std::bitset<6> neumann)
{
  map<unsigned long, set<Block>> bc_to_blocks;
  for (auto iface : iface_domain.getInterfaces()) {
    int i = iface->global_index;
    for (auto patch : iface->patches) {
      Side<3> aux = patch.side;
      const PatchIfaceInfo<3>& sinfo = *patch.piinfo;
      IfaceType<3> type = patch.type;
      for (Side<3> s : Side<3>::getValues()) {
        if (sinfo.pinfo.hasNbr(s)) {
          int j = sinfo.getIfaceInfo(s)->global_index;
          std::bitset<6> patch_neumann;
          for (Side<3> s : Side<3>::getValues()) {
            patch_neumann[s.getIndex()] = neumann[s.getIndex()] && !sinfo.pinfo.hasNbr(s);
          }
          Block block(s, j, aux, i, patch_neumann, type);
          bc_to_blocks[block.non_dirichlet_boundary.to_ulong()].insert(block);
        }
      }
    }
  }
  vector<set<Block>> blocks_vector;
  blocks_vector.reserve(bc_to_blocks.size());
  for (auto& pair : bc_to_blocks) {
    blocks_vector.push_back(std::move(pair.second));
  }
  return blocks_vector;
}
/**
 * @brief Fill the coefficients for the blocks
 *
 * @tparam CoeffMap The map from block to coefficients
 * @param coeffs The map from block to coefficients
 * @param pinfo the patchinfo that has to be solved on
 * @param solver the patch solver
 */
template<class CoeffMap>
void
FillBlockCoeffs(CoeffMap coeffs, const PatchInfo<3>& pinfo, Poisson::FFTWPatchSolver<3>& solver)
{
  auto ns = solver.getDomain().getNs();
  int n = ns[0];
  const MPIGhostFiller<3>& ghost_filler =
    dynamic_cast<const MPIGhostFiller<3>&>(solver.getGhostFiller());
  for (int yi = 0; yi < n; yi++) {
    for (int xi = 0; xi < n; xi++) {
      int j = xi + yi * n;
      // create some work vectors
      auto u_vec = make_shared<Vector<3>>(Communicator(MPI_COMM_SELF), ns, 1, 1, 1);
      auto f_vec = make_shared<Vector<3>>(Communicator(MPI_COMM_SELF), ns, 1, 1, 1);
      PatchView<double, 3> u_view = u_vec->getPatchView(0);
      View<double, 3> u_west_ghosts = u_view.getGhostSliceOn(Side<3>::west(), { 0 });
      PatchView<double, 3> f_view = f_vec->getPatchView(0);

      u_west_ghosts(xi, yi, 0) = 2;

      solver.solveSinglePatch(pinfo, f_view, u_view);

      u_west_ghosts(xi, yi, 0) = 0;

      for (const auto& pair : coeffs) {
        Side<3> s = pair.first.aux;
        IfaceType<3> type = pair.first.type;
        vector<double>& block = *pair.second;
        if (type.isNormal()) {
          FillBlockColumnForNormalInterface(j, u_view, s, block);
        } else if (type.isCoarseToCoarse()) {
          FillBlockColumnForCoarseToCoarseInterface(j, u_view, s, ghost_filler, pinfo, block);

        } else if (type.isFineToFine()) {
          FillBlockColumnForFineToFineInterface(j, u_view, s, ghost_filler, pinfo, type, block);

        } else if (type.isCoarseToFine()) {
          FillBlockColumnForCoarseToFineInterface(j, u_view, s, ghost_filler, pinfo, type, block);

        } else if (type.isFineToCoarse()) {
          FillBlockColumnForFineToCoarseInterface(j, u_view, s, ghost_filler, pinfo, type, block);
        }

        if (s == Side<3>::west()) {
          if (type.isNormal()) {
            block[n * n * j + j] += 0.5;
          } else if (type.isFineToFine() || type.isCoarseToCoarse()) {
            block[n * n * j + j] += 1;
          }
        }
      }
    }
  }
}
template<class Inserter>
void
AssembleMatrix(const Schur::InterfaceDomain<3>& iface_domain,
               Poisson::FFTWPatchSolver<3>& solver,
               Inserter insertBlock)
{
  auto ns = iface_domain.getDomain().getNs();
  int n = ns[0];

  for (const set<Block>& blocks : GetBlocks(iface_domain, solver.getNeumann())) {
    // create domain representing curr_type
    PatchInfo<3> pinfo;
    pinfo.setNbrInfo(Side<3>::west(), new NormalNbrInfo<2>());
    pinfo.ns.fill(n);
    pinfo.spacings.fill(1.0 / n);
    pinfo.num_ghost_cells = 1;
    for (Side<3> s : Side<3>::getValues()) {
      if (!blocks.begin()->non_dirichlet_boundary[s.getIndex()]) {
        pinfo.setNbrInfo(s, new NormalNbrInfo<2>());
      }
    }
    solver.addPatch(pinfo);

    map<Block, shared_ptr<vector<double>>, std::function<bool(const Block& a, const Block& b)>>
      coeffs([](const Block& a, const Block& b) {
        return std::tie(a.aux, a.type) < std::tie(b.aux, b.type);
      });
    // allocate blocks of coefficients
    for (const Block& b : blocks) {
      shared_ptr<vector<double>> ptr = coeffs[b];
      if (ptr.get() == nullptr) {
        coeffs[b] = shared_ptr<vector<double>>(new vector<double>(n * n * n * n));
      }
    }

    FillBlockCoeffs(coeffs, pinfo, solver);

    // now insert these results into the matrix for each interface
    for (Block block : blocks) {
      insertBlock(block, coeffs[block]);
    }
  }
}
const function<int(int, int, int)> transforms_left[4] = {
  [](int n, int xi, int yi) { return xi + yi * n; },
  [](int n, int xi, int yi) { return n - yi - 1 + xi * n; },
  [](int n, int xi, int yi) { return n - xi - 1 + (n - yi - 1) * n; },
  [](int n, int xi, int yi) { return yi + (n - xi - 1) * n; }
};
const function<int(int, int, int)> transforms_right[4] = {
  [](int n, int xi, int yi) { return xi + yi * n; },
  [](int n, int xi, int yi) { return yi + (n - xi - 1) * n; },
  [](int n, int xi, int yi) { return n - xi - 1 + (n - yi - 1) * n; },
  [](int n, int xi, int yi) { return n - yi - 1 + xi * n; }
};
const function<int(int, int, int)> transforms_left_inv[4] = { [](int n, int xi, int yi) {
                                                               xi = n - xi - 1;
                                                               return xi + yi * n;
                                                             },
                                                              [](int n, int xi, int yi) {
                                                                xi = n - xi - 1;
                                                                return n - yi - 1 + xi * n;
                                                              },
                                                              [](int n, int xi, int yi) {
                                                                xi = n - xi - 1;
                                                                return n - xi - 1 +
                                                                       (n - yi - 1) * n;
                                                              },
                                                              [](int n, int xi, int yi) {
                                                                xi = n - xi - 1;
                                                                return yi + (n - xi - 1) * n;
                                                              } };
const function<int(int, int, int)> transforms_right_inv[4] = { [](int n, int xi, int yi) {
                                                                xi = n - xi - 1;
                                                                return xi + yi * n;
                                                              },
                                                               [](int n, int xi, int yi) {
                                                                 xi = n - xi - 1;
                                                                 return yi + (n - xi - 1) * n;
                                                               },
                                                               [](int n, int xi, int yi) {
                                                                 xi = n - xi - 1;
                                                                 return n - xi - 1 +
                                                                        (n - yi - 1) * n;
                                                               },
                                                               [](int n, int xi, int yi) {
                                                                 xi = n - xi - 1;
                                                                 return n - yi - 1 + xi * n;
                                                               } };
std::function<int(int, int, int)>
GetColTransform(const Block& b)
{
  function<int(int, int, int)> col_trans;
  if (sideIsLeftOriented(b.main)) {
    if (b.mainFlipped()) {
      col_trans = transforms_left_inv[b.main_rotation];
    } else {
      col_trans = transforms_left[b.main_rotation];
    }
  } else {
    if (b.mainFlipped()) {
      col_trans = transforms_right_inv[b.main_rotation];
    } else {
      col_trans = transforms_right[b.main_rotation];
    }
  }
  return col_trans;
}
std::function<int(int, int, int)>
GetRowTransform(const Block& b)
{
  function<int(int, int, int)> row_trans;
  if (sideIsLeftOriented(b.aux)) {
    if (b.auxFlipped()) {
      row_trans = transforms_left_inv[b.aux_rotation];
    } else {
      row_trans = transforms_left[b.aux_rotation];
    }
  } else {
    if (b.auxFlipped()) {
      row_trans = transforms_right_inv[b.aux_rotation];
    } else {
      row_trans = transforms_right[b.aux_rotation];
    }
  }
  return row_trans;
}
vector<double>
FlipBlock(int n, const Block& b, const vector<double>& orig)
{
  vector<double> copy(n * n * n * n);
  function<int(int, int, int)> col_trans = GetColTransform(b);
  function<int(int, int, int)> row_trans = GetRowTransform(b);

  for (int row_yi = 0; row_yi < n; row_yi++) {
    for (int row_xi = 0; row_xi < n; row_xi++) {
      int i_dest = row_xi + row_yi * n;
      int i_orig = row_trans(n, row_xi, row_yi);
      for (int col_yi = 0; col_yi < n; col_yi++) {
        for (int col_xi = 0; col_xi < n; col_xi++) {
          int j_dest = col_xi + col_yi * n;
          int j_orig = col_trans(n, col_xi, col_yi);

          copy[i_dest * n * n + j_dest] = orig[i_orig * n * n + j_orig];
        }
      }
    }
  }
  return copy;
}
} // namespace
Mat
ThunderEgg::Poisson::FastSchurMatrixAssemble3D(const Schur::InterfaceDomain<3>& iface_domain,
                                               Poisson::FFTWPatchSolver<3>& solver)
{
  auto ns = iface_domain.getDomain().getNs();
  if (ns[0] != ns[1] && ns[0] != ns[2]) {
    throw RuntimeError("FastSchurMatrixAssemble3D only supports cube shaped patches");
  }
  const GhostFiller<3>& gf = solver.getGhostFiller();
  if (typeid(gf) != typeid(TriLinearGhostFiller)) {
    throw RuntimeError("FastSchurMatrixAssemble3D only supports TriLinearGhostFiller");
  }
  Mat A;
  MatCreate(MPI_COMM_WORLD, &A);
  int n = ns[0];
  int local_size = iface_domain.getNumLocalInterfaces() * n * n;
  int global_size = iface_domain.getNumGlobalInterfaces() * n * n;
  MatSetSizes(A, local_size, local_size, global_size, global_size);
  MatSetType(A, MATMPIAIJ);
  MatMPIAIJSetPreallocation(A, 26 * n * n, nullptr, 26 * n * n, nullptr);

  auto insertBlock = [&](const Block& b, shared_ptr<vector<double>> coeffs) {
    int global_i = b.i * n * n;
    int global_j = b.j * n * n;

    vector<double> copy = FlipBlock(n, b, *coeffs);

    vector<int> inds_i(n * n);
    iota(inds_i.begin(), inds_i.end(), global_i);
    vector<int> inds_j(n * n);
    iota(inds_j.begin(), inds_j.end(), global_j);

    MatSetValues(A, n * n, &inds_i[0], n * n, &inds_j[0], &copy[0], ADD_VALUES);
  };

  AssembleMatrix(iface_domain, solver, insertBlock);
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  return A;
}