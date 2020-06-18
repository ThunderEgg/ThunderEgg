/***************************************************************************
 *  Thunderegg, a library for solving Poisson's equation on adaptively
 *  refined block-structured Cartesian grids
 *
 *  Copyright (C) 2019  Thunderegg Developers. See AUTHORS.md file at the
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

#ifndef THUNDEREGG_MPIGHOSTFILLER_H
#define THUNDEREGG_MPIGHOSTFILLER_H

#include <Thunderegg/GhostFiller.h>
namespace Thunderegg
{
template <size_t D> class MPIGhostFiller : public GhostFiller<D>
{
	private:
	using remote_call = std::tuple<std::shared_ptr<const PatchInfo<D>>, const Side<D>,
	                               const NbrType, const Orthant<D>>;
	using LocalCall = std::tuple<std::shared_ptr<const PatchInfo<D>>, const Side<D>, const NbrType,
	                             const Orthant<D>, int, int>;

	std::deque<remote_call> remote_calls;
	std::deque<LocalCall>   local_calls;
	std::map<int, size_t>   rank_buff_size_map;

	virtual void fillGhostCellsForNbrPatch(std::shared_ptr<const PatchInfo<D>> pinfo,
	                                       const LocalData<D>                  local_data,
	                                       const LocalData<D> nbr_data, const Side<D> side,
	                                       const NbrType    nbr_type,
	                                       const Orthant<D> orthant) const = 0;

	virtual void fillGhostCellsForLocalPatch(std::shared_ptr<const PatchInfo<D>> pinfo,
	                                         const LocalData<D> local_data) const = 0;

	protected:
	std::shared_ptr<const Domain<D>> domain;
	int                              side_cases;

	public:
	MPIGhostFiller(std::shared_ptr<const Domain<D>> domain_in, int side_cases_in)
	: domain(domain_in), side_cases(side_cases_in)
	{
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		using risop = std::tuple<int, int, const Side<D>, const Orthant<D>, size_t>;
		// std::map < riso,
		std::set<remote_call> remote_calls_set;
		for (auto pinfo : domain->getPatchInfoVector()) {
			for (Side<D> s : Side<D>::getValues()) {
				if (pinfo->hasNbr(s)) {
					switch (pinfo->getNbrType(s)) {
						case NbrType::Normal: {
							auto nbrinfo = pinfo->getNormalNbrInfo(s);
							if (nbrinfo.rank == rank) {
								local_calls.emplace_back(pinfo, s, NbrType::Normal,
								                         Orthant<D>::null(), pinfo->local_index,
								                         nbrinfo.local_index);
							} else {
							}
						} break;
						case NbrType::Fine: {
							auto nbrinfo  = pinfo->getFineNbrInfo(s);
							auto orthants = Orthant<D>::getValuesOnSide(s);
							for (size_t i = 0; i < orthants.size(); i++) {
								if (nbrinfo.ranks[i] == rank) {
									local_calls.emplace_back(pinfo, s, NbrType::Fine, orthants[i],
									                         pinfo->local_index,
									                         nbrinfo.local_indexes[i]);
								}
							}
						} break;
						case NbrType::Coarse: {
							auto nbrinfo = pinfo->getCoarseNbrInfo(s);
							auto orthant = Orthant<D>::getValuesOnSide(
							s.opposite())[nbrinfo.orth_on_coarse.toInt()];
							if (nbrinfo.rank == rank) {
								local_calls.emplace_back(pinfo, s, NbrType::Coarse, orthant,
								                         pinfo->local_index, nbrinfo.local_index);
							}
						} break;
					}
				}
			}
		}
	}
	/**
	 * @brief Fill ghost cells on a vector
	 *
	 * @param u  the vector
	 */
	void fillGhost(std::shared_ptr<const Vector<D>> u) const
	{
		for (auto pinfo : domain->getPatchInfoVector()) {
			const LocalData<2> this_patch = u->getLocalData(pinfo->local_index);
			for (Side<2> s : Side<2>::getValues()) {
				if (pinfo->hasNbr(s)) {
					for (int i = 0; i < pinfo->num_ghost_cells; i++) {
						auto this_ghost = this_patch.getGhostSliceOnSide(s, i + 1);
						nested_loop<1>(
						this_ghost.getStart(), this_ghost.getEnd(),
						[&](const std::array<int, 1> &coord) { this_ghost[coord] = 0; });
					}
				}
			}
		}
		for (auto pinfo : domain->getPatchInfoVector()) {
			auto data = u->getLocalData(pinfo->local_index);
			fillGhostCellsForLocalPatch(pinfo, data);
		}
		for (const LocalCall &call : local_calls) {
			auto pinfo      = std::get<0>(call);
			auto side       = std::get<1>(call);
			auto nbr_type   = std::get<2>(call);
			auto orthant    = std::get<3>(call);
			auto local_data = u->getLocalData(std::get<4>(call));
			auto nbr_data   = u->getLocalData(std::get<5>(call));
			fillGhostCellsForNbrPatch(pinfo, local_data, nbr_data, side, nbr_type, orthant);
		}
	}
};
} // namespace Thunderegg
#endif