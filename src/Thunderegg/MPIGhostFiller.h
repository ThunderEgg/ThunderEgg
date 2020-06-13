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
	virtual void fillGhostCellsForNbrPatch(std::shared_ptr<const PatchInfo<D>> pinfo,
	                                       const LocalData<D>                  local_data,
	                                       const LocalData<D> nbr_data, const Side<D> side,
	                                       const NbrType    nbr_type,
	                                       const Orthant<D> orthant) const = 0;

	virtual void fillGhostCellsForLocalPatch(std::shared_ptr<const PatchInfo<D>> pinfo,
	                                         const LocalData<D> local_data) const = 0;

	protected:
	std::shared_ptr<Domain<D>> domain;
	int                        side_cases;

	public:
	MPIGhostFiller(std::shared_ptr<Domain<D>> domain_in, int side_cases_in)
	: domain(domain_in), side_cases(side_cases_in)
	{
	}
	/**
	 * @brief Fill ghost cells on a vector
	 *
	 * @param u  the vector
	 */
	void fillGhost(std::shared_ptr<const Vector<D>> u) const
	{
		for (auto pinfo : domain->getPatchInfoVector()) {
			auto data = u->getLocalData(pinfo->local_index);
			fillGhostCellsForLocalPatch(pinfo, data);

			for (Side<D> s : Side<D>::getValues()) {
				if (pinfo->hasNbr(s)) {
					switch (pinfo->getNbrType(s)) {
						case NbrType::Normal: {
							auto nbrinfo  = pinfo->getNormalNbrInfo(s);
							auto nbr_data = u->getLocalData(nbrinfo.local_index);
							fillGhostCellsForNbrPatch(pinfo, data, nbr_data, s, NbrType::Normal, 0);
						} break;
						case NbrType::Fine: {
							auto nbrinfo  = pinfo->getFineNbrInfo(s);
							auto orthants = Orthant<D>::getValuesOnSide(s);
							for (size_t i = 0; i < orthants.size(); i++) {
								auto nbr_data = u->getLocalData(nbrinfo.local_indexes[i]);
								fillGhostCellsForNbrPatch(pinfo, data, nbr_data, s, NbrType::Fine,
								                          orthants[i]);
							}
						} break;
						case NbrType::Coarse: {
							auto nbrinfo  = pinfo->getCoarseNbrInfo(s);
							auto nbr_data = u->getLocalData(nbrinfo.local_index);
							auto orthant  = Orthant<D>::getValuesOnSide(
                            s.opposite())[nbrinfo.orth_on_coarse.toInt()];
							fillGhostCellsForNbrPatch(pinfo, data, nbr_data, s, NbrType::Coarse,
							                          orthant);
						} break;
					}
				}
			}
		}
	}
};
} // namespace Thunderegg
#endif