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
#ifndef THUNDEREGG_COARSENBRINFO_H
#define THUNDEREGG_COARSENBRINFO_H
#include <ThunderEgg/BufferReader.h>
#include <ThunderEgg/BufferWriter.h>
#include <ThunderEgg/NbrInfo.h>
#include <ThunderEgg/Orthant.h>

namespace ThunderEgg
{
/**
 * @brief Represents a neighbor that is at a coarser refinement level.
 *
 * @tparam D the number of Cartesian dimensions.
 */
template <size_t D> class CoarseNbrInfo : public NbrInfo<D>
{
	public:
	/**
	 * @brief The mpi rank that the neighbor resides on.
	 */
	int rank = 0;
	/**
	 * @brief The id of the neighbor
	 */
	int id = 0;
	/**
	 * @brief The local index of the neighbor
	 */
	int local_index = 0;
	/**
	 * @brief The global index of the neighbor
	 */
	int global_index = 0;
	/**
	 * @brief The orthant that this patch in relation to the coarser patch's interface.
	 */
	Orthant<D - 1> orth_on_coarse;
	/**
	 * @brief Construct a new empty CoarseNbrInfo object
	 */
	CoarseNbrInfo()  = default;
	~CoarseNbrInfo() = default;
	/**
	 * @brief Construct a new CoarseNbrInfo object
	 *
	 * @param id the id of the neighbor
	 * @param orth_on_coarse The orthant that this patch in relation to the coarser
	 * patch's interface.
	 */
	CoarseNbrInfo(int id, Orthant<D - 1> orth_on_coarse)
	{
		this->id             = id;
		this->orth_on_coarse = orth_on_coarse;
	}
	NbrType getNbrType()
	{
		return NbrType::Coarse;
	}
	void getNbrIds(std::deque<int> &nbr_ids)
	{
		nbr_ids.push_back(id);
	};
	void getNbrRanks(std::deque<int> &nbr_ranks)
	{
		nbr_ranks.push_back(rank);
	}
	void setGlobalIndexes(std::map<int, int> &rev_map)
	{
		global_index = rev_map.at(local_index);
	}
	void setLocalIndexes(std::map<int, int> &rev_map)
	{
		local_index = rev_map.at(id);
	}
	int serialize(char *buffer) const
	{
		BufferWriter writer(buffer);
		writer << rank;
		writer << id;
		writer << orth_on_coarse;
		return writer.getPos();
	}
	int deserialize(char *buffer)
	{
		BufferReader reader(buffer);
		reader >> rank;
		reader >> id;
		reader >> orth_on_coarse;
		return reader.getPos();
	}
};
} // namespace ThunderEgg
#endif