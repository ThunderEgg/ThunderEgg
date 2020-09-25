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

#include "FastSchurMatrixAssemble2D.h"
#include <ThunderEgg/BiLinearGhostFiller.h>
#include <ThunderEgg/BiQuadraticGhostFiller.h>
#include <ThunderEgg/PETSc/VecWrapper.h>
#include <numeric>
#include <typeinfo>
using namespace std;
using namespace ThunderEgg;
using namespace ThunderEgg::Schur;
namespace
{
/**
 * @brief Given The interface side, and an auxilary side,
 */
const Side<2> rot_table[4][4]
= {{Side<2>::west(), Side<2>::east(), Side<2>::south(), Side<2>::north()},
   {Side<2>::east(), Side<2>::west(), Side<2>::north(), Side<2>::south()},
   {Side<2>::north(), Side<2>::south(), Side<2>::west(), Side<2>::east()},
   {Side<2>::south(), Side<2>::north(), Side<2>::east(), Side<2>::west()}};
/**
 * @brief true of the the j indexes have to be flipped for an interfaces side
 */
const bitset<4> flip_j_table = 0b0110;
/**
 * @brief true of the the i indexes have to be flipped for an interfaces side and auxilary side
 */
const array<bitset<4>, 4> flip_i_table = {{0b0000, 0b1111, 0b0011, 0b1100}};
/**
 * @brief Represents a block in the Schur compliment matrix
 */
class Block
{
	public:
	/**
	 * @brief The interface type
	 */
	IfaceType<2> type;
	/**
	 * @brief The side of the patch that the block is on
	 */
	Side<2> s;
	/**
	 * @brief True if boundary on a given side is not Dirichlet
	 */
	bitset<4> non_dirichlet_boundary;
	/**
	 * @brief i block index in matrix
	 */
	int i;
	/**
	 * @brief j block index in matrix
	 */
	int j;
	/**
	 * @brief true if j indexes are flipped
	 */
	bool flip_i;
	/**
	 * @brief true if j indexes are flipped
	 */
	bool flip_j;

	/**
	 * @brief Construct a new Block object
	 *
	 * @param main the side of the patch that the interface being set is on
	 * @param j the block j index in the matrix
	 * @param aux the side of the patch that the interface being affected is on
	 * @param i the block i index in the matrix
	 * @param non_dirichlet_boundary true if side of patch has non dirichlet boundary conditions
	 * @param type the type of interface
	 */
	Block(Side<2> main, int j, Side<2> aux, int i, bitset<4> non_dirichlet_boundary,
	      IfaceType<2> type)
	: type(type), i(i), j(j)
	{
		s = rot_table[main.getIndex()][aux.getIndex()];
		for (int side_index = 0; side_index < 4; side_index++) {
			this->non_dirichlet_boundary[rot_table[main.getIndex()][side_index].getIndex()]
			= non_dirichlet_boundary[side_index];
		}
		flip_j = flip_j_table[main.getIndex()];
		flip_i = flip_i_table[main.getIndex()][s.getIndex()];
		if (flip_i) {
			this->type.setOrthant(type.getOrthant().getNbrOnSide(Side<1>::west()));
		}
	}
	bool operator==(const Block &b) const
	{
		return non_dirichlet_boundary.to_ulong() == b.non_dirichlet_boundary.to_ulong();
	}
	bool operator<(const Block &b) const
	{
		return std::tie(type, i, j, flip_j) < std::tie(b.type, b.i, b.j, b.flip_j);
	}
};
/**
 * @brief Get the LocalData object for the buffer
 *
 * @param buffer_ptr pointer to the ghost cells position in the buffer
 * @param pinfo  the PatchInfo object
 * @param side  the side that the ghost cells are on
 * @return LocalData<D> the LocalData object
 */
LocalData<2> getLocalDataForBuffer(double *buffer_ptr, shared_ptr<const PatchInfo<2>> pinfo,
                                   const Side<2> side)
{
	auto ns              = pinfo->ns;
	int  num_ghost_cells = pinfo->num_ghost_cells;
	// determine striding
	std::array<int, 2> strides;
	strides[0] = 1;
	for (size_t i = 1; i < 2; i++) {
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
		transformed_buffer_ptr
		= buffer_ptr - ns[side.getAxisIndex()] * strides[side.getAxisIndex()];
	}

	LocalData<2> buffer_data(transformed_buffer_ptr, strides, ns, num_ghost_cells);
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
void FillBlockColumnForNormalInterface(int j, const LocalData<2> &u, Side<2> s,
                                       std::vector<double> &block)
{
	int  n     = u.getLengths()[0];
	auto slice = u.getSliceOnSide(s);
	for (int i = 0; i < n; i++) {
		block[i * n + j] = -slice[{i}] / 2;
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
void FillBlockColumnForCoarseToCoarseInterface(
int j, const LocalData<2> &u, Side<2> s, std::shared_ptr<const MPIGhostFiller<2>> ghost_filler,
std::shared_ptr<const PatchInfo<2>> pinfo, std::vector<double> &block)
{
	int  n                            = pinfo->ns[0];
	auto new_pinfo                    = make_shared<PatchInfo<2>>(*pinfo);
	new_pinfo->nbr_info[0]            = nullptr;
	new_pinfo->nbr_info[s.getIndex()] = make_shared<FineNbrInfo<2>>();
	ghost_filler->fillGhostCellsForLocalPatch(new_pinfo, u);
	auto slice       = u.getSliceOnSide(s);
	auto ghost_slice = u.getGhostSliceOnSide(s, 1);
	for (int i = 0; i < n; i++) {
		block[i * n + j] = -(slice[{i}] + ghost_slice[{i}]) / 2;
		ghost_slice[{i}] = 0;
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
void FillBlockColumnForFineToFineInterface(int j, const LocalData<2> &u, Side<2> s,
                                           std::shared_ptr<const MPIGhostFiller<2>> ghost_filler,
                                           std::shared_ptr<const PatchInfo<2>>      pinfo,
                                           IfaceType<2> type, std::vector<double> &block)
{
	int  n                            = pinfo->ns[0];
	auto new_pinfo                    = make_shared<PatchInfo<2>>(*pinfo);
	new_pinfo->nbr_info[0]            = nullptr;
	new_pinfo->nbr_info[s.getIndex()] = make_shared<CoarseNbrInfo<2>>(100, type.getOrthant());
	ghost_filler->fillGhostCellsForLocalPatch(new_pinfo, u);
	auto slice       = u.getSliceOnSide(s);
	auto ghost_slice = u.getGhostSliceOnSide(s, 1);
	for (int i = 0; i < n; i++) {
		block[i * n + j] = -(slice[{i}] + ghost_slice[{i}]) / 2;
		ghost_slice[{i}] = 0;
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
void FillBlockColumnForCoarseToFineInterface(int j, const LocalData<2> &u, Side<2> s,
                                             std::shared_ptr<const MPIGhostFiller<2>> ghost_filler,
                                             std::shared_ptr<const PatchInfo<2>>      pinfo,
                                             IfaceType<2> type, std::vector<double> &block)
{
	int  n                            = pinfo->ns[0];
	auto new_pinfo                    = make_shared<PatchInfo<2>>(*pinfo);
	new_pinfo->nbr_info[0]            = nullptr;
	new_pinfo->nbr_info[s.getIndex()] = make_shared<FineNbrInfo<2>>();
	vector<double> ghosts(n);
	LocalData<2>   nbr_data = getLocalDataForBuffer(ghosts.data(), pinfo, s.opposite());
	ghost_filler->fillGhostCellsForNbrPatch(
	new_pinfo, u, nbr_data, s, NbrType::Fine,
	Orthant<2>::getValuesOnSide(s)[type.getOrthant().getIndex()]);
	for (int i = 0; i < n; i++) {
		block[i * n + j] = -ghosts[i] / 2;
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
void FillBlockColumnForFineToCoarseInterface(int j, const LocalData<2> &u, Side<2> s,
                                             std::shared_ptr<const MPIGhostFiller<2>> ghost_filler,
                                             std::shared_ptr<const PatchInfo<2>>      pinfo,
                                             IfaceType<2> type, std::vector<double> &block)
{
	int  n                            = pinfo->ns[0];
	auto new_pinfo                    = make_shared<PatchInfo<2>>(*pinfo);
	new_pinfo->nbr_info[0]            = nullptr;
	new_pinfo->nbr_info[s.getIndex()] = make_shared<CoarseNbrInfo<2>>(100, type.getOrthant());
	vector<double> ghosts(n);
	LocalData<2>   nbr_data = getLocalDataForBuffer(ghosts.data(), pinfo, s.opposite());
	ghost_filler->fillGhostCellsForNbrPatch(
	new_pinfo, u, nbr_data, s, NbrType::Coarse,
	Orthant<2>::getValuesOnSide(s.opposite())[type.getOrthant().getIndex()]);
	for (int i = 0; i < n; i++) {
		block[i * n + j] = -ghosts[i] / 2;
	}
}
/**
 * @brief Get a vector of set<Block>, each set of blocks have the same boundary conditions
 *
 * @param iface_domain the InterfaceDomain
 * @return vector<set<Block>> the vector of blocks
 */
vector<set<Block>> GetBlocks(shared_ptr<const InterfaceDomain<2>> iface_domain)
{
	map<unsigned long, set<Block>> bc_to_blocks;
	for (auto iface : iface_domain->getInterfaces()) {
		int i = iface->global_index;
		for (auto patch : iface->patches) {
			Side<2>                  aux   = patch.side;
			const PatchIfaceInfo<2> &sinfo = *patch.piinfo;
			IfaceType<2>             type  = patch.type;
			for (Side<2> s : Side<2>::getValues()) {
				if (sinfo.pinfo->hasNbr(s)) {
					int   j = sinfo.getIfaceInfo(s)->global_index;
					Block block(s, j, aux, i, sinfo.pinfo->neumann, type);
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
void FillBlockCoeffs(CoeffMap coeffs, std::shared_ptr<const PatchInfo<2>> pinfo,
                     std::shared_ptr<Poisson::FFTWPatchSolver<2>> solver)
{
	auto ns           = solver->getDomain()->getNs();
	int  n            = ns[0];
	auto ghost_filler = dynamic_pointer_cast<const MPIGhostFiller<2>>(solver->getGhostFiller());
	for (int j = 0; j < n; j++) {
		// create some work vectors
		auto         u_vec         = make_shared<ValVector<2>>(MPI_COMM_SELF, ns, 1, 1);
		auto         f_vec         = make_shared<ValVector<2>>(MPI_COMM_SELF, ns, 1, 1);
		LocalData<2> u_local_data  = u_vec->getLocalData(0);
		LocalData<1> u_west_ghosts = u_local_data.getGhostSliceOnSide(Side<2>::west(), 1);
		LocalData<2> f_local_data  = f_vec->getLocalData(0);

		u_west_ghosts[{j}] = 2;

		solver->solveSinglePatch(pinfo, u_local_data, f_local_data);

		for (const auto &pair : coeffs) {
			Side<2>         s    = pair.first.s;
			IfaceType<2>    type = pair.first.type;
			vector<double>  filled_ghosts(n);
			vector<double> &block = *pair.second;
			if (type.isNormal()) {
				FillBlockColumnForNormalInterface(j, u_local_data, s, block);
			} else if (type.isCoarseToCoarse()) {
				FillBlockColumnForCoarseToCoarseInterface(j, u_local_data, s, ghost_filler, pinfo,
				                                          block);

			} else if (type.isFineToFine()) {
				FillBlockColumnForFineToFineInterface(j, u_local_data, s, ghost_filler, pinfo, type,
				                                      block);

			} else if (type.isCoarseToFine()) {
				FillBlockColumnForCoarseToFineInterface(j, u_local_data, s, ghost_filler, pinfo,
				                                        type, block);

			} else if (type.isFineToCoarse()) {
				FillBlockColumnForFineToCoarseInterface(j, u_local_data, s, ghost_filler, pinfo,
				                                        type, block);
			}

			if (s == Side<2>::west()) {
				if (type.isNormal()) {
					block[n * j + j] += 0.5;
				} else if (type.isFineToFine() || type.isCoarseToCoarse()) {
					block[n * j + j] += 1;
				}
			}
		}
	}
}
/**
 * @brief Assemble the matrix
 *
 * @tparam Inserter has the follow arguments
 *  (int block_i, int block_j, shared_ptr<vector<double>> block, bool flip_i, bool flip_j)
 * @param iface_domain the InterfaceDomain
 * @param solver the PatchSolver
 * @param insertBlock the Inserter
 */
template <class Inserter>
void assembleMatrix(std::shared_ptr<const InterfaceDomain<2>>    iface_domain,
                    std::shared_ptr<Poisson::FFTWPatchSolver<2>> solver, Inserter insertBlock)
{
	auto ns = iface_domain->getDomain()->getNs();
	int  n  = ns[0];

	for (const set<Block> &blocks : GetBlocks(iface_domain)) {
		// create domain representing curr_type
		auto pinfo             = make_shared<PatchInfo<2>>();
		pinfo->nbr_info[0]     = make_shared<NormalNbrInfo<2>>();
		pinfo->num_ghost_cells = 1;
		auto piinfo            = make_shared<PatchIfaceInfo<2>>(pinfo);
		piinfo->pinfo          = pinfo;
		pinfo->ns.fill(n);
		pinfo->spacings.fill(1.0 / n);
		pinfo->neumann = blocks.begin()->non_dirichlet_boundary;
		piinfo->setIfaceInfo(Side<2>::west(),
		                     make_shared<NormalIfaceInfo<2>>(pinfo, Side<2>::west()));

		solver->addPatch(pinfo);

		std::vector<std::shared_ptr<PatchIfaceInfo<2>>> single_domain;
		single_domain.push_back(piinfo);

		// coefficients are grouped by block's side and type
		map<Block, shared_ptr<vector<double>>, std::function<bool(const Block &a, const Block &b)>>
		coeffs([](const Block &a, const Block &b) {
			return std::tie(a.s, a.type) < std::tie(b.s, b.type);
		});

		// allocate blocks of coefficients
		for (const Block &b : blocks) {
			shared_ptr<vector<double>> ptr = coeffs[b];
			if (ptr.get() == nullptr) {
				coeffs[b] = shared_ptr<vector<double>>(new vector<double>(n * n));
			}
		}

		FillBlockCoeffs(coeffs, pinfo, solver);

		// now insert these results into the matrix for each interface
		for (Block block : blocks) {
			insertBlock(block.i, block.j, coeffs[block], block.flip_i, block.flip_j);
		}
	}
}
} // namespace
Mat ThunderEgg::Poisson::FastSchurMatrixAssemble2D(
std::shared_ptr<const InterfaceDomain<2>>    iface_domain,
std::shared_ptr<Poisson::FFTWPatchSolver<2>> solver)
{
	array<int, 2> ns = iface_domain->getDomain()->getNs();
	if (ns[0] != ns[1]) {
		throw RuntimeError("FastSchurMatrixAssembler2D does not support non-square patches");
	}
	if (dynamic_pointer_cast<const BiLinearGhostFiller>(solver->getGhostFiller()) == nullptr
	    && dynamic_pointer_cast<const BiQuadraticGhostFiller>(solver->getGhostFiller())
	       == nullptr) {
		throw RuntimeError(
		"FastSchurMatrixAssembler2D only supports BiLinearGhostFiller and BiQuadraticGhostFiller");
	}
	int n = ns[0];
	Mat A;
	MatCreate(MPI_COMM_WORLD, &A);
	int local_size  = iface_domain->getNumLocalInterfaces() * n;
	int global_size = iface_domain->getNumGlobalInterfaces() * n;
	MatSetSizes(A, local_size, local_size, global_size, global_size);
	MatSetType(A, MATMPIAIJ);
	MatMPIAIJSetPreallocation(A, 19 * n, nullptr, 19 * n, nullptr);

	auto insertBlock
	= [&](int block_i, int block_j, shared_ptr<vector<double>> block, bool flip_i, bool flip_j) {
		  int matrix_i = block_i * n;
		  int matrix_j = block_j * n;

		  vector<double> &orig = *block;
		  vector<double>  copy(n * n);
		  for (int i = 0; i < n; i++) {
			  int orig_i = i;
			  if (flip_i) {
				  orig_i = n - i - 1;
			  }
			  for (int j = 0; j < n; j++) {
				  int orig_j = j;
				  if (flip_j) {
					  orig_j = n - j - 1;
				  }
				  copy[i * n + j] = orig[orig_i * n + orig_j];
			  }
		  }
		  vector<int> inds_i(n);
		  iota(inds_i.begin(), inds_i.end(), matrix_i);
		  vector<int> inds_j(n);
		  iota(inds_j.begin(), inds_j.end(), matrix_j);

		  MatSetValues(A, n, &inds_i[0], n, &inds_j[0], &copy[0], ADD_VALUES);
	  };

	assembleMatrix(iface_domain, solver, insertBlock);
	MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
	return A;
}
