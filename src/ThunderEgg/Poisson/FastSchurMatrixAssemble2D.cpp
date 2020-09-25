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
Side<2> rot_table[4][4] = {{Side<2>::west(), Side<2>::east(), Side<2>::south(), Side<2>::north()},
                           {Side<2>::east(), Side<2>::west(), Side<2>::north(), Side<2>::south()},
                           {Side<2>::north(), Side<2>::south(), Side<2>::west(), Side<2>::east()},
                           {Side<2>::south(), Side<2>::north(), Side<2>::east(), Side<2>::west()}};
/**
 * @brief Represents a block in the Schur compliment matrix
 */
struct Block {
	IfaceType<2> type;
	Side<2>      s;
	/**
	 * @brief True if boundary on a given side is not dirichlet
	 */
	bitset<4>           non_dirichlet_boundary;
	int                 i;
	int                 j;
	bool                flip_i;
	bool                flip_j;
	bitset<4>           flip_j_table = 0b0110;
	array<bitset<4>, 4> flip_i_table = {{0b0000, 0b1111, 0b0011, 0b1100}};

	Block(Side<2> main, int j, Side<2> aux, int i, bitset<4> non_dirichlet_boundary,
	      IfaceType<2> type)
	: type(type)
	{
		s       = rot_table[main.getIndex()][aux.getIndex()];
		this->i = i;
		this->j = j;
		for (int s = 0; s < 4; s++) {
			this->non_dirichlet_boundary[rot_table[main.getIndex()][s].getIndex()]
			= non_dirichlet_boundary[s];
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
struct BlockKey {
	IfaceType<2> type;
	Side<2>      s;

	BlockKey(const Block &b) : type(b.type)
	{
		s = b.s;
	}
	friend bool operator<(const BlockKey &l, const BlockKey &r)
	{
		return std::tie(l.s, l.type) < std::tie(r.s, r.type);
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
void assembleMatrix(
std::shared_ptr<const InterfaceDomain<2>>                                       iface_domain,
std::shared_ptr<Poisson::FFTWPatchSolver<2>>                                    solver,
std::function<void(int, int, std::shared_ptr<std::vector<double>>, bool, bool)> insertBlock)
{
	auto ns = iface_domain->getDomain()->getNs();
	int  n  = ns[0];

	auto ghost_filler = dynamic_pointer_cast<const MPIGhostFiller<2>>(solver->getGhostFiller());

	// get block type
	set<Block> blocks;
	for (auto iface : iface_domain->getInterfaces()) {
		int i = iface->global_index;
		for (auto patch : iface->patches) {
			Side<2>                  aux   = patch.side;
			const PatchIfaceInfo<2> &sinfo = *patch.piinfo;
			IfaceType<2>             type  = patch.type;
			for (int s = 0; s < 4; s++) {
				if (sinfo.iface_info[s] != nullptr) {
					int     j    = sinfo.iface_info[s]->global_index;
					Side<2> main = static_cast<Side<2>>(s);
					blocks.insert(Block(main, j, aux, i, sinfo.pinfo->neumann, type));
				}
			}
		}
	}
	int          num_types     = 0;
	auto         u_vec         = make_shared<ValVector<2>>(MPI_COMM_SELF, ns, 1, 1);
	auto         f_vec         = make_shared<ValVector<2>>(MPI_COMM_SELF, ns, 1, 1);
	LocalData<2> u_local_data  = u_vec->getLocalData(0);
	LocalData<1> u_west_ghosts = u_local_data.getGhostSliceOnSide(Side<2>::west(), 1);
	LocalData<2> f_local_data  = f_vec->getLocalData(0);
	while (!blocks.empty()) {
		num_types++;
		// the first in the set is the type of interface that we are going to solve for
		set<Block> todo;
		Block      curr_type = *blocks.begin();
		blocks.erase(blocks.begin());
		todo.insert(curr_type);
		set<Block> to_be_deleted;
		for (auto iter = blocks.begin(); iter != blocks.end(); iter++) {
			if (*iter == curr_type) {
				todo.insert(*iter);
				to_be_deleted.insert(*iter);
			}
		}
		for (Block i : to_be_deleted) {
			blocks.erase(i);
		}

		// create domain representing curr_type
		auto pinfo             = make_shared<PatchInfo<2>>();
		pinfo->nbr_info[0]     = make_shared<NormalNbrInfo<2>>();
		pinfo->num_ghost_cells = 1;
		auto piinfo            = make_shared<PatchIfaceInfo<2>>(pinfo);
		piinfo->pinfo          = pinfo;
		pinfo->ns.fill(n);
		pinfo->spacings.fill(1.0 / n);
		pinfo->neumann = curr_type.non_dirichlet_boundary;
		piinfo->setIfaceInfo(Side<2>::west(),
		                     make_shared<NormalIfaceInfo<2>>(pinfo, Side<2>::west()));

		solver->addPatch(pinfo);

		std::vector<std::shared_ptr<PatchIfaceInfo<2>>> single_domain;
		single_domain.push_back(piinfo);

		map<BlockKey, shared_ptr<vector<double>>> coeffs;
		// allocate blocks of coefficients
		for (const Block &b : todo) {
			shared_ptr<vector<double>> ptr = coeffs[b];
			if (ptr.get() == nullptr) {
				coeffs[b] = shared_ptr<vector<double>>(new vector<double>(n * n));
			}
		}

		for (int j = 0; j < n; j++) {
			u_vec->set(0);
			u_west_ghosts[{j}] = 2;
			solver->solveSinglePatch(pinfo, u_local_data, f_local_data);
			u_west_ghosts[{j}] = 0;

			// fill the blocks
			for (auto &p : coeffs) {
				Side<2>         s    = p.first.s;
				IfaceType<2>    type = p.first.type;
				vector<double>  filled_ghosts(n);
				vector<double> &block = *p.second;
				if (type.isNormal()) {
					auto slice = u_local_data.getSliceOnSide(s);
					for (int i = 0; i < n; i++) {
						block[i * n + j] = slice[{i}] / 2;
					}
				} else if (type.isCoarseToCoarse()) {
					auto new_pinfo                    = make_shared<PatchInfo<2>>(*pinfo);
					new_pinfo->nbr_info[0]            = nullptr;
					new_pinfo->nbr_info[s.getIndex()] = make_shared<FineNbrInfo<2>>();
					ghost_filler->fillGhostCellsForLocalPatch(new_pinfo, u_local_data);
					auto slice       = u_local_data.getSliceOnSide(s);
					auto ghost_slice = u_local_data.getGhostSliceOnSide(s, 1);
					for (int i = 0; i < n; i++) {
						block[i * n + j] = (slice[{i}] + ghost_slice[{i}]) / 2;
						ghost_slice[{i}] = 0;
					}
				} else if (type.isFineToFine()) {
					auto new_pinfo         = make_shared<PatchInfo<2>>(*pinfo);
					new_pinfo->nbr_info[0] = nullptr;
					new_pinfo->nbr_info[s.getIndex()]
					= make_shared<CoarseNbrInfo<2>>(100, type.getOrthant());
					ghost_filler->fillGhostCellsForLocalPatch(new_pinfo, u_local_data);
					auto slice       = u_local_data.getSliceOnSide(s);
					auto ghost_slice = u_local_data.getGhostSliceOnSide(s, 1);
					for (int i = 0; i < n; i++) {
						block[i * n + j] = (slice[{i}] + ghost_slice[{i}]) / 2;
						ghost_slice[{i}] = 0;
					}
				} else if (type.isCoarseToFine()) {
					auto new_pinfo                    = make_shared<PatchInfo<2>>(*pinfo);
					new_pinfo->nbr_info[0]            = nullptr;
					new_pinfo->nbr_info[s.getIndex()] = make_shared<FineNbrInfo<2>>();
					vector<double> ghosts(n);
					LocalData<2>   nbr_data
					= getLocalDataForBuffer(ghosts.data(), pinfo, s.opposite());
					ghost_filler->fillGhostCellsForNbrPatch(
					new_pinfo, u_local_data, nbr_data, s, NbrType::Fine,
					Orthant<2>::getValuesOnSide(s)[type.getOrthant().getIndex()]);
					for (int i = 0; i < n; i++) {
						block[i * n + j] = ghosts[i] / 2;
					}
				} else if (type.isFineToCoarse()) {
					auto new_pinfo         = make_shared<PatchInfo<2>>(*pinfo);
					new_pinfo->nbr_info[0] = nullptr;
					new_pinfo->nbr_info[s.getIndex()]
					= make_shared<CoarseNbrInfo<2>>(100, type.getOrthant());
					vector<double> ghosts(n);
					LocalData<2>   nbr_data
					= getLocalDataForBuffer(ghosts.data(), pinfo, s.opposite());
					ghost_filler->fillGhostCellsForNbrPatch(
					new_pinfo, u_local_data, nbr_data, s, NbrType::Coarse,
					Orthant<2>::getValuesOnSide(s.opposite())[type.getOrthant().getIndex()]);
					for (int i = 0; i < n; i++) {
						block[i * n + j] = ghosts[i] / 2;
					}
				}
				// interpolator->interpolate(*piinfo, s, 0, type, u_vec, interp_vec);

				if (s == Side<2>::west()) {
					if (type.isNormal()) {
						block[n * j + j] -= 0.5;
					} else if (type.isFineToFine() || type.isCoarseToCoarse()) {
						block[n * j + j] -= 1;
					}
				}
			}
		}

		// now insert these results into the matrix for each interface
		for (Block block : todo) {
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
	int local_size  = iface_domain->getInterfaces().size() * n;
	int global_size = iface_domain->getNumGlobalInterfaces() * n;
	MatSetSizes(A, local_size, local_size, global_size, global_size);
	MatSetType(A, MATMPIAIJ);
	MatMPIAIJSetPreallocation(A, 19 * n, nullptr, 19 * n, nullptr);

	auto insertBlock
	= [&](int i, int j, shared_ptr<vector<double>> block, bool flip_i, bool flip_j) {
		  int local_i = i * n;
		  int local_j = j * n;

		  vector<double> &orig = *block;
		  vector<double>  copy(n * n);
		  for (int i = 0; i < n; i++) {
			  int block_i = i;
			  if (flip_i) {
				  block_i = n - i - 1;
			  }
			  for (int j = 0; j < n; j++) {
				  int block_j = j;
				  if (flip_j) {
					  block_j = n - j - 1;
				  }
				  copy[i * n + j] = orig[block_i * n + block_j];
			  }
		  }
		  vector<int> inds_i(n);
		  iota(inds_i.begin(), inds_i.end(), local_i);
		  vector<int> inds_j(n);
		  iota(inds_j.begin(), inds_j.end(), local_j);

		  MatSetValues(A, n, &inds_i[0], n, &inds_j[0], &copy[0], ADD_VALUES);
	  };

	assembleMatrix(iface_domain, solver, insertBlock);
	MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
	return A;
}
