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

#include "FastSchurMatrixAssemble3D.h"
#include <ThunderEgg/MPIGhostFiller.h>
#include <ThunderEgg/PETSc/VecWrapper.h>
#include <ThunderEgg/TriLinearGhostFiller.h>
using namespace std;
using namespace ThunderEgg;
using namespace ThunderEgg::Schur;
namespace
{
enum class Rotation : char { x_cw, x_ccw, y_cw, y_ccw, z_cw, z_ccw };
bool sideIsLeftOriented(const Side<3> s)
{
	return (s == Side<3>::north() || s == Side<3>::west() || s == Side<3>::bottom());
}
class Block
{
	public:
	static const Side<3>          side_table[6][6];
	static const char             rots_table[6][6];
	static const vector<Rotation> main_rot_plan[6];
	static const vector<Rotation> aux_rot_plan_dirichlet[6];
	static const vector<Rotation> aux_rot_plan_neumann[16];
	static const char             rot_quad_lookup_left[4][4];
	static const char             rot_quad_lookup_right[4][4];
	static const char             quad_flip_lookup[4];
	IfaceType<3>                  type;
	Side<3>                       main;
	Side<3>                       aux;
	int                           j;
	int                           i;
	bitset<6>                     non_dirichlet_boundary;
	bool                          orig_main_is_left_oriented;
	bool                          orig_aux_is_left_oriented;
	unsigned char                 main_rotation = 0;
	unsigned char                 aux_rotation  = 0;
	Block(Side<3> main, int j, Side<3> aux, int i, bitset<6> non_dirichlet_boundary, IfaceType<3> type)
	: type(type),
	  main(main),
	  aux(aux),
	  j(j),
	  i(i),
	  non_dirichlet_boundary(non_dirichlet_boundary),
	  orig_main_is_left_oriented(sideIsLeftOriented(main)),
	  orig_aux_is_left_oriented(sideIsLeftOriented(aux))
	{
		rotate();
	}
	void applyRotation(const Rotation rot)
	{
		// main rotation
		main_rotation = (main_rotation + rots_table[static_cast<int>(rot)][main.getIndex()]) % 4;
		// aux rotation
		aux_rotation          = (aux_rotation + rots_table[static_cast<int>(rot)][aux.getIndex()]) % 4;
		main                  = side_table[static_cast<int>(rot)][main.getIndex()];
		aux                   = side_table[static_cast<int>(rot)][aux.getIndex()];
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
			int quad = (int) type.getOrthant().getIndex();
			quad     = rotateQuad(quad);
			type.setOrthant(Orthant<2>((unsigned char) quad));
		}
	}
	bool operator<(const Block &b) const
	{
		return std::tie(i, j, main_rotation, orig_main_is_left_oriented, aux_rotation)
		       < std::tie(b.i, b.j, b.main_rotation, b.orig_main_is_left_oriented, b.aux_rotation);
	}
	bool mainFlipped() const
	{
		return sideIsLeftOriented(main) != orig_main_is_left_oriented;
	}
	bool auxFlipped() const
	{
		return sideIsLeftOriented(aux) != orig_aux_is_left_oriented;
	}
};

const Side<3> Block::side_table[6][6] = {{Side<3>::west(), Side<3>::east(), Side<3>::top(), Side<3>::bottom(), Side<3>::south(), Side<3>::north()},
                                         {Side<3>::west(), Side<3>::east(), Side<3>::bottom(), Side<3>::top(), Side<3>::north(), Side<3>::south()},
                                         {Side<3>::bottom(), Side<3>::top(), Side<3>::south(), Side<3>::north(), Side<3>::east(), Side<3>::west()},
                                         {Side<3>::top(), Side<3>::bottom(), Side<3>::south(), Side<3>::north(), Side<3>::west(), Side<3>::east()},
                                         {Side<3>::north(), Side<3>::south(), Side<3>::west(), Side<3>::east(), Side<3>::bottom(), Side<3>::top()},
                                         {Side<3>::south(), Side<3>::north(), Side<3>::east(), Side<3>::west(), Side<3>::bottom(), Side<3>::top()}};
const char    Block::rots_table[6][6]
= {{3, 1, 0, 0, 2, 2}, {1, 3, 2, 2, 0, 0}, {1, 3, 3, 1, 1, 3}, {1, 3, 1, 3, 3, 1}, {0, 0, 0, 0, 3, 1}, {0, 0, 0, 0, 1, 3}};
const vector<Rotation> Block::main_rot_plan[6]
= {{}, {Rotation::z_cw, Rotation::z_cw}, {Rotation::z_cw}, {Rotation::z_ccw}, {Rotation::y_ccw}, {Rotation::y_cw}};
const vector<Rotation> Block::aux_rot_plan_dirichlet[6]   = {{}, {}, {}, {Rotation::x_cw, Rotation::x_cw}, {Rotation::x_cw}, {Rotation::x_ccw}};
const vector<Rotation> Block::aux_rot_plan_neumann[16]    = {{},
                                                          {},
                                                          {Rotation::x_cw, Rotation::x_cw},
                                                          {},
                                                          {Rotation::x_cw},
                                                          {Rotation::x_cw},
                                                          {Rotation::x_cw, Rotation::x_cw},
                                                          {Rotation::x_cw, Rotation::x_cw},
                                                          {Rotation::x_ccw},
                                                          {},
                                                          {Rotation::x_ccw},
                                                          {},
                                                          {Rotation::x_ccw},
                                                          {Rotation::x_cw},
                                                          {Rotation::x_ccw},
                                                          {}};
const char             Block::rot_quad_lookup_left[4][4]  = {{0, 1, 2, 3}, {1, 3, 0, 2}, {3, 2, 1, 0}, {2, 0, 3, 1}};
const char             Block::rot_quad_lookup_right[4][4] = {{0, 1, 2, 3}, {2, 0, 3, 1}, {3, 2, 1, 0}, {1, 3, 0, 2}};
const char             Block::quad_flip_lookup[4]         = {1, 0, 3, 2};

/**
 * @brief Get the LocalData object for the buffer
 *
 * @param buffer_ptr pointer to the ghost cells position in the buffer
 * @param pinfo  the PatchInfo object
 * @param side  the side that the ghost cells are on
 * @return LocalData<D> the LocalData object
 */
LocalData<3> getLocalDataForBuffer(double *buffer_ptr, shared_ptr<const PatchInfo<3>> pinfo, const Side<3> side)
{
	auto ns              = pinfo->ns;
	int  num_ghost_cells = pinfo->num_ghost_cells;
	// determine striding
	std::array<int, 3> strides;
	strides[0] = 1;
	for (size_t i = 1; i < 3; i++) {
		if (i == side.getAxisIndex() + 1) {
			strides[i] = num_ghost_cells * strides[i - 1];
		} else {
			strides[i] = ns[i - 1] * strides[i - 1];
		}
	}
	// transform buffer ptr so that it points to first non-ghost cell
	double *transformed_buffer_ptr;
	if (side.isLowerOnAxis()) {
		transformed_buffer_ptr = buffer_ptr - (-num_ghost_cells) * strides[side.getAxisIndex()];
	} else {
		transformed_buffer_ptr = buffer_ptr - ns[side.getAxisIndex()] * strides[side.getAxisIndex()];
	}

	LocalData<3> buffer_data(transformed_buffer_ptr, strides, ns, num_ghost_cells);
	return buffer_data;
}
/**
 * @brief Fill a block column for a normal interface
 *
 * @param j the column of the block to fill
 * @param u the patch data
 * @param s the side of the patch that the block is on
 * @param block
 */
void FillBlockColumnForNormalInterface(int j, const LocalData<3> &u, Side<3> s, std::vector<double> &block)
{
	int  n     = u.getLengths()[0];
	auto slice = u.getSliceOnSide(s);
	for (int yi = 0; yi < n; yi++) {
		for (int xi = 0; xi < n; xi++) {
			block[(xi + n * yi) * n * n + j] = -slice[{xi, yi}] / 2;
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
void FillBlockColumnForCoarseToCoarseInterface(int                                      j,
                                               const LocalData<3> &                     u,
                                               Side<3>                                  s,
                                               std::shared_ptr<const MPIGhostFiller<3>> ghost_filler,
                                               std::shared_ptr<const PatchInfo<3>>      pinfo,
                                               std::vector<double> &                    block)
{
	int  n                            = pinfo->ns[0];
	auto new_pinfo                    = make_shared<PatchInfo<3>>(*pinfo);
	new_pinfo->nbr_info[0]            = nullptr;
	new_pinfo->nbr_info[s.getIndex()] = make_unique<FineNbrInfo<3>>();
	std::vector<LocalData<3>> us      = {u};
	ghost_filler->fillGhostCellsForLocalPatch(new_pinfo, us);
	auto slice       = u.getSliceOnSide(s);
	auto ghost_slice = u.getGhostSliceOnSide(s, 1);
	for (int yi = 0; yi < n; yi++) {
		for (int xi = 0; xi < n; xi++) {
			block[(xi + yi * n) * n * n + j] = -(slice[{xi, yi}] + ghost_slice[{xi, yi}]) / 2;
			ghost_slice[{xi, yi}]            = 0;
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
void FillBlockColumnForFineToFineInterface(int                                      j,
                                           const LocalData<3> &                     u,
                                           Side<3>                                  s,
                                           std::shared_ptr<const MPIGhostFiller<3>> ghost_filler,
                                           std::shared_ptr<const PatchInfo<3>>      pinfo,
                                           IfaceType<3>                             type,
                                           std::vector<double> &                    block)
{
	int  n                            = pinfo->ns[0];
	auto new_pinfo                    = make_shared<PatchInfo<3>>(*pinfo);
	new_pinfo->nbr_info[0]            = nullptr;
	new_pinfo->nbr_info[s.getIndex()] = make_unique<CoarseNbrInfo<3>>(100, type.getOrthant());
	std::vector<LocalData<3>> us      = {u};
	ghost_filler->fillGhostCellsForLocalPatch(new_pinfo, us);
	auto slice       = u.getSliceOnSide(s);
	auto ghost_slice = u.getGhostSliceOnSide(s, 1);
	for (int yi = 0; yi < n; yi++) {
		for (int xi = 0; xi < n; xi++) {
			block[(xi + yi * n) * n * n + j] = -(slice[{xi, yi}] + ghost_slice[{xi, yi}]) / 2;
			ghost_slice[{xi, yi}]            = 0;
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
void FillBlockColumnForCoarseToFineInterface(int                                      j,
                                             const LocalData<3> &                     u,
                                             Side<3>                                  s,
                                             std::shared_ptr<const MPIGhostFiller<3>> ghost_filler,
                                             std::shared_ptr<const PatchInfo<3>>      pinfo,
                                             IfaceType<3>                             type,
                                             std::vector<double> &                    block)
{
	int  n                            = pinfo->ns[0];
	auto new_pinfo                    = make_shared<PatchInfo<3>>(*pinfo);
	new_pinfo->nbr_info[0]            = nullptr;
	new_pinfo->nbr_info[s.getIndex()] = make_unique<FineNbrInfo<3>>();
	vector<double>            ghosts(n * n);
	std::vector<LocalData<3>> us        = {u};
	std::vector<LocalData<3>> nbr_datas = {getLocalDataForBuffer(ghosts.data(), pinfo, s.opposite())};
	std::vector<Side<3>>      sides     = {s};
	ghost_filler->fillGhostCellsForNbrPatch(
	new_pinfo, us, nbr_datas, sides, NbrType::Fine, Orthant<3>::getValuesOnSide(s)[type.getOrthant().getIndex()]);
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
void FillBlockColumnForFineToCoarseInterface(int                                      j,
                                             const LocalData<3> &                     u,
                                             Side<3>                                  s,
                                             std::shared_ptr<const MPIGhostFiller<3>> ghost_filler,
                                             std::shared_ptr<const PatchInfo<3>>      pinfo,
                                             IfaceType<3>                             type,
                                             std::vector<double> &                    block)
{
	int  n                            = pinfo->ns[0];
	auto new_pinfo                    = make_shared<PatchInfo<3>>(*pinfo);
	new_pinfo->nbr_info[0]            = nullptr;
	new_pinfo->nbr_info[s.getIndex()] = make_unique<CoarseNbrInfo<3>>(100, type.getOrthant());
	vector<double>            ghosts(n * n);
	std::vector<LocalData<3>> us        = {u};
	std::vector<LocalData<3>> nbr_datas = {getLocalDataForBuffer(ghosts.data(), pinfo, s.opposite())};
	std::vector<Side<3>>      sides     = {s};
	ghost_filler->fillGhostCellsForNbrPatch(
	new_pinfo, us, nbr_datas, sides, NbrType::Coarse, Orthant<3>::getValuesOnSide(s.opposite())[type.getOrthant().getIndex()]);
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
vector<set<Block>> GetBlocks(shared_ptr<const InterfaceDomain<3>> iface_domain, std::bitset<6> neumann)
{
	map<unsigned long, set<Block>> bc_to_blocks;
	for (auto iface : iface_domain->getInterfaces()) {
		int i = iface->global_index;
		for (auto patch : iface->patches) {
			Side<3>                  aux   = patch.side;
			const PatchIfaceInfo<3> &sinfo = *patch.piinfo;
			IfaceType<3>             type  = patch.type;
			for (Side<3> s : Side<3>::getValues()) {
				if (sinfo.pinfo->hasNbr(s)) {
					int            j = sinfo.getIfaceInfo(s)->global_index;
					std::bitset<6> patch_neumann;
					for (Side<3> s : Side<3>::getValues()) {
						patch_neumann[s.getIndex()] = neumann[s.getIndex()] && !sinfo.pinfo->hasNbr(s);
					}
					Block block(s, j, aux, i, patch_neumann, type);
					bc_to_blocks[block.non_dirichlet_boundary.to_ulong()].insert(block);
				}
			}
		}
	}
	vector<set<Block>> blocks_vector;
	blocks_vector.reserve(bc_to_blocks.size());
	for (auto &pair : bc_to_blocks) {
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
template <class CoeffMap>
void FillBlockCoeffs(CoeffMap coeffs, std::shared_ptr<const PatchInfo<3>> pinfo, std::shared_ptr<Poisson::FFTWPatchSolver<3>> solver)
{
	auto ns           = solver->getDomain()->getNs();
	int  n            = ns[0];
	auto ghost_filler = dynamic_pointer_cast<const MPIGhostFiller<3>>(solver->getGhostFiller());
	for (int yi = 0; yi < n; yi++) {
		for (int xi = 0; xi < n; xi++) {
			int j = xi + yi * n;
			// create some work vectors
			auto         u_vec         = make_shared<ValVector<3>>(MPI_COMM_SELF, ns, 1, 1, 1);
			auto         f_vec         = make_shared<ValVector<3>>(MPI_COMM_SELF, ns, 1, 1, 1);
			LocalData<3> u_local_data  = u_vec->getLocalData(0, 0);
			auto         u_local_datas = u_vec->getLocalDatas(0);
			LocalData<2> u_west_ghosts = u_local_data.getGhostSliceOnSide(Side<3>::west(), 1);
			LocalData<3> f_local_data  = f_vec->getLocalData(0, 0);
			auto         f_local_datas = f_vec->getLocalDatas(0);

			u_west_ghosts[{xi, yi}] = 2;

			solver->solveSinglePatch(pinfo, f_local_datas, u_local_datas);

			for (const auto &pair : coeffs) {
				Side<3>         s     = pair.first.aux;
				IfaceType<3>    type  = pair.first.type;
				vector<double> &block = *pair.second;
				if (type.isNormal()) {
					FillBlockColumnForNormalInterface(j, u_local_data, s, block);
				} else if (type.isCoarseToCoarse()) {
					FillBlockColumnForCoarseToCoarseInterface(j, u_local_data, s, ghost_filler, pinfo, block);

				} else if (type.isFineToFine()) {
					FillBlockColumnForFineToFineInterface(j, u_local_data, s, ghost_filler, pinfo, type, block);

				} else if (type.isCoarseToFine()) {
					FillBlockColumnForCoarseToFineInterface(j, u_local_data, s, ghost_filler, pinfo, type, block);

				} else if (type.isFineToCoarse()) {
					FillBlockColumnForFineToCoarseInterface(j, u_local_data, s, ghost_filler, pinfo, type, block);
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
template <class Inserter>
void AssembleMatrix(std::shared_ptr<const Schur::InterfaceDomain<3>> iface_domain,
                    std::shared_ptr<Poisson::FFTWPatchSolver<3>>     solver,
                    Inserter                                         insertBlock)
{
	auto ns = iface_domain->getDomain()->getNs();
	int  n  = ns[0];

	for (const set<Block> &blocks : GetBlocks(iface_domain, solver->getNeumann())) {
		// create domain representing curr_type
		std::shared_ptr<PatchInfo<3>> pinfo(new PatchInfo<3>());
		pinfo->nbr_info[0] = make_unique<NormalNbrInfo<3>>();
		pinfo->ns.fill(n);
		pinfo->spacings.fill(1.0 / n);
		pinfo->num_ghost_cells = 1;
		for (Side<3> s : Side<3>::getValues()) {
			if (!blocks.begin()->non_dirichlet_boundary[s.getIndex()]) {
				pinfo->nbr_info[s.getIndex()] = make_unique<NormalNbrInfo<3>>();
			}
		}
		solver->addPatch(*pinfo);

		map<Block, shared_ptr<vector<double>>, std::function<bool(const Block &a, const Block &b)>> coeffs(
		[](const Block &a, const Block &b) { return std::tie(a.aux, a.type) < std::tie(b.aux, b.type); });
		// allocate blocks of coefficients
		for (const Block &b : blocks) {
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
const function<int(int, int, int)> transforms_left[4]      = {[](int n, int xi, int yi) { return xi + yi * n; },
                                                         [](int n, int xi, int yi) { return n - yi - 1 + xi * n; },
                                                         [](int n, int xi, int yi) { return n - xi - 1 + (n - yi - 1) * n; },
                                                         [](int n, int xi, int yi) { return yi + (n - xi - 1) * n; }};
const function<int(int, int, int)> transforms_right[4]     = {[](int n, int xi, int yi) { return xi + yi * n; },
                                                          [](int n, int xi, int yi) { return yi + (n - xi - 1) * n; },
                                                          [](int n, int xi, int yi) { return n - xi - 1 + (n - yi - 1) * n; },
                                                          [](int n, int xi, int yi) { return n - yi - 1 + xi * n; }};
const function<int(int, int, int)> transforms_left_inv[4]  = {[](int n, int xi, int yi) {
                                                                 xi = n - xi - 1;
                                                                 return xi + yi * n;
                                                             },
                                                             [](int n, int xi, int yi) {
                                                                 xi = n - xi - 1;
                                                                 return n - yi - 1 + xi * n;
                                                             },
                                                             [](int n, int xi, int yi) {
                                                                 xi = n - xi - 1;
                                                                 return n - xi - 1 + (n - yi - 1) * n;
                                                             },
                                                             [](int n, int xi, int yi) {
                                                                 xi = n - xi - 1;
                                                                 return yi + (n - xi - 1) * n;
                                                             }};
const function<int(int, int, int)> transforms_right_inv[4] = {[](int n, int xi, int yi) {
	                                                              xi = n - xi - 1;
	                                                              return xi + yi * n;
                                                              },
                                                              [](int n, int xi, int yi) {
	                                                              xi = n - xi - 1;
	                                                              return yi + (n - xi - 1) * n;
                                                              },
                                                              [](int n, int xi, int yi) {
	                                                              xi = n - xi - 1;
	                                                              return n - xi - 1 + (n - yi - 1) * n;
                                                              },
                                                              [](int n, int xi, int yi) {
	                                                              xi = n - xi - 1;
	                                                              return n - yi - 1 + xi * n;
                                                              }};
std::function<int(int, int, int)>  GetColTransform(const Block &b)
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
std::function<int(int, int, int)> GetRowTransform(const Block &b)
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
vector<double> FlipBlock(int n, const Block &b, const vector<double> &orig)
{
	vector<double>               copy(n * n * n * n);
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
Mat ThunderEgg::Poisson::FastSchurMatrixAssemble3D(std::shared_ptr<const Schur::InterfaceDomain<3>> iface_domain,
                                                   std::shared_ptr<Poisson::FFTWPatchSolver<3>>     solver)
{
	auto ns = iface_domain->getDomain()->getNs();
	if (ns[0] != ns[1] && ns[0] != ns[2]) {
		throw RuntimeError("FastSchurMatrixAssemble3D only supports cube shaped patches");
	}
	if (dynamic_pointer_cast<const TriLinearGhostFiller>(solver->getGhostFiller()) == nullptr) {
		throw RuntimeError("FastSchurMatrixAssemble3D only supports TriLinearGhostFiller");
	}
	Mat A;
	MatCreate(MPI_COMM_WORLD, &A);
	int n           = ns[0];
	int local_size  = iface_domain->getNumLocalInterfaces() * n * n;
	int global_size = iface_domain->getNumGlobalInterfaces() * n * n;
	MatSetSizes(A, local_size, local_size, global_size, global_size);
	MatSetType(A, MATMPIAIJ);
	MatMPIAIJSetPreallocation(A, 26 * n * n, nullptr, 26 * n * n, nullptr);

	auto insertBlock = [&](const Block &b, shared_ptr<vector<double>> coeffs) {
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