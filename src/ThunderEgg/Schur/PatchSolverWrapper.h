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

#ifndef THUNDEREGG_SCHUR_PATCHSOLVERWRAPPER_H
#define THUNDEREGG_SCHUR_PATCHSOLVERWRAPPER_H

#include <ThunderEgg/PatchSolver.h>
#include <ThunderEgg/Schur/PatchIfaceScatter.h>
#include <ThunderEgg/ValVectorGenerator.h>

namespace ThunderEgg
{
namespace Schur
{
/**
 * @brief Creates a Schur compliment matrix operator for an InterfaceDomain by using a PatchSolver.
 *
 * @tparam D the number of Cartesian dimensions
 */
template <int D> class PatchSolverWrapper : public Operator<D - 1>
{
	private:
	/**
	 * @brief The InterfaceDomain
	 */
	std::shared_ptr<const InterfaceDomain<D>> iface_domain;
	/**
	 * @brief The PatchSolver that is being wrapped
	 */
	std::shared_ptr<const PatchSolver<D>> solver;
	/**
	 * @brief The scatter object
	 */
	mutable PatchIfaceScatter<D> scatter;
	/**
	 * @brief A VectorGenerator for the domain
	 *
	 */
	ThunderEgg::ValVectorGenerator<D> vg;
	/**
	 * @brief Set of patches that have only local interfaces
	 */
	std::deque<std::shared_ptr<const PatchIfaceInfo<D>>> patches_with_only_local_ifaces;
	/**
	 * @brief Set of patches that have interafce on neighboring ranks
	 */
	std::deque<std::shared_ptr<const PatchIfaceInfo<D>>> patches_with_ifaces_on_neighbor_rank;

	public:
	/**
	 * @brief Construct a new PatchSolverWrapper object
	 *
	 * @param iface_domain the InterfaceDomain for the Schur compliment system
	 * @param solver the PatchSolver to wrap
	 */
	PatchSolverWrapper(std::shared_ptr<const InterfaceDomain<D>> iface_domain, std::shared_ptr<const PatchSolver<D>> solver)
	: iface_domain(iface_domain),
	  solver(solver),
	  scatter(iface_domain),
	  vg(solver->getDomain(), 1)
	{
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		for (auto piinfo : iface_domain->getPatchIfaceInfos()) {
			for (Side<D> s : Side<D>::getValues()) {
				if (piinfo->pinfo.hasNbr(s) && piinfo->getIfaceInfo(s)->rank != rank) {
					patches_with_ifaces_on_neighbor_rank.push_back(piinfo);
					break;
				}
			}
			if (patches_with_ifaces_on_neighbor_rank.empty() || patches_with_ifaces_on_neighbor_rank.back() != piinfo) {
				patches_with_only_local_ifaces.push_back(piinfo);
			}
		}
	}
	/**
	 * @brief Apply Schur matrix
	 *
	 * @param x the input vector.
	 * @param b the output vector.
	 */
	void apply(std::shared_ptr<const Vector<D - 1>> x, std::shared_ptr<Vector<D - 1>> b) const override
	{
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		auto u = vg.getNewVector();
		auto f = vg.getNewVector();
		// scatter local iface vector
		auto local_x = scatter.getNewLocalPatchIfaceVector();
		scatter.scatterStart(x, local_x);

		// each patch has a local interface, go ahead and set the ghost values using those
		// interfaces
		for (auto piinfo : iface_domain->getPatchIfaceInfos()) {
			for (Side<D> s : Side<D>::getValues()) {
				auto local_data = u->getLocalData(0, piinfo->pinfo.local_index);
				if (piinfo->pinfo.hasNbr(s) && piinfo->getIfaceInfo(s)->rank == rank) {
					auto ghosts    = local_data.getGhostSliceOnSide(s, 1);
					auto interface = local_x->getLocalData(0, piinfo->getIfaceInfo(s)->patch_local_index);
					nested_loop<D - 1>(
					interface.getStart(), interface.getEnd(), [&](const std::array<int, D - 1> &coord) { ghosts[coord] = 2 * interface[coord]; });
				}
			}
		}
		// go ahead and solve for patches with only local interfaces
		for (auto piinfo : patches_with_only_local_ifaces) {
			auto us = u->getLocalDatas(piinfo->pinfo.local_index);
			auto fs = f->getLocalDatas(piinfo->pinfo.local_index);
			solver->solveSinglePatch(piinfo->pinfo, fs, us);
		}

		scatter.scatterFinish(x, local_x);

		// set ghosts using interfaces that were on a neighboring rank
		for (auto piinfo : patches_with_ifaces_on_neighbor_rank) {
			for (Side<D> s : Side<D>::getValues()) {
				auto local_data = u->getLocalData(0, piinfo->pinfo.local_index);
				if (piinfo->pinfo.hasNbr(s) && piinfo->getIfaceInfo(s)->rank != rank) {
					auto ghosts    = local_data.getGhostSliceOnSide(s, 1);
					auto interface = local_x->getLocalData(0, piinfo->getIfaceInfo(s)->patch_local_index);
					nested_loop<D - 1>(
					interface.getStart(), interface.getEnd(), [&](const std::array<int, D - 1> &coord) { ghosts[coord] = 2 * interface[coord]; });
				}
			}
		}
		// solve the remaining patches
		for (auto piinfo : patches_with_ifaces_on_neighbor_rank) {
			auto us = u->getLocalDatas(piinfo->pinfo.local_index);
			auto fs = f->getLocalDatas(piinfo->pinfo.local_index);
			solver->solveSinglePatch(piinfo->pinfo, fs, us);
		}

		solver->getGhostFiller()->fillGhost(u);

		for (auto iface : iface_domain->getInterfaces()) {
			for (auto patch : iface->patches) {
				if (patch.piinfo->pinfo.rank == rank && (patch.type.isNormal() || patch.type.isCoarseToCoarse() || patch.type.isFineToFine())) {
					auto local_data = u->getLocalData(0, patch.piinfo->pinfo.local_index);
					auto ghosts     = local_data.getGhostSliceOnSide(patch.side, 1);
					auto inner      = local_data.getSliceOnSide(patch.side);
					auto interface  = b->getLocalData(0, patch.piinfo->getIfaceInfo(patch.side)->patch_local_index);
					nested_loop<D - 1>(interface.getStart(), interface.getEnd(), [&](const std::array<int, D - 1> &coord) {
						interface[coord] = (ghosts[coord] + inner[coord]) / 2;
					});
					break;
				}
			}
		}
		b->scaleThenAdd(-1, x);
	}
	/**
	 * @brief Get the RHS for the Schur system from a given RHS for the domain system
	 *
	 * @param domain_b the domain rhs
	 * @param schur_b the Schur rhs
	 */
	void getSchurRHSFromDomainRHS(std::shared_ptr<const Vector<D>> domain_b, std::shared_ptr<Vector<D - 1>> schur_b) const
	{
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		auto u = vg.getNewVector();

		for (auto piinfo : iface_domain->getPatchIfaceInfos()) {
			for (Side<D> s : Side<D>::getValues()) {
				auto local_data = u->getLocalData(0, piinfo->pinfo.local_index);
				if (piinfo->pinfo.hasNbr(s)) {
					auto ghosts = local_data.getGhostSliceOnSide(s, 1);
					auto inner  = local_data.getSliceOnSide(s);
					nested_loop<D - 1>(
					ghosts.getStart(), ghosts.getEnd(), [&](const std::array<int, D - 1> &coord) { ghosts[coord] = -inner[coord]; });
				}
			}
		}
		for (auto piinfo : iface_domain->getPatchIfaceInfos()) {
			auto us = u->getLocalDatas(piinfo->pinfo.local_index);
			auto fs = domain_b->getLocalDatas(piinfo->pinfo.local_index);
			solver->solveSinglePatch(piinfo->pinfo, fs, us);
		}

		solver->getGhostFiller()->fillGhost(u);

		for (auto iface : iface_domain->getInterfaces()) {
			for (auto patch : iface->patches) {
				if (patch.piinfo->pinfo.rank == rank && (patch.type.isNormal() || patch.type.isCoarseToCoarse() || patch.type.isFineToFine())) {
					auto local_data = u->getLocalData(0, patch.piinfo->pinfo.local_index);
					auto ghosts     = local_data.getGhostSliceOnSide(patch.side, 1);
					auto inner      = local_data.getSliceOnSide(patch.side);
					auto interface  = schur_b->getLocalData(0, patch.piinfo->getIfaceInfo(patch.side)->patch_local_index);
					nested_loop<D - 1>(interface.getStart(), interface.getEnd(), [&](const std::array<int, D - 1> &coord) {
						interface[coord] = (ghosts[coord] + inner[coord]) / 2;
					});
					break;
				}
			}
		}
	}
};
extern template class PatchSolverWrapper<2>;
extern template class PatchSolverWrapper<3>;
} // namespace Schur
} // namespace ThunderEgg
#endif
