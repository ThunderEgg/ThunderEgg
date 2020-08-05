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
#ifndef THUNDEREGG_FINENBRINFO_H
#define THUNDEREGG_FINENBRINFO_H
#include <ThunderEgg/BufferReader.h>
#include <ThunderEgg/BufferWriter.h>
#include <ThunderEgg/NbrInfo.h>
#include <ThunderEgg/Orthant.h>

namespace ThunderEgg
{
/**
 * @brief Represents neighbors that are at a finer refinement level.
 *
 * @tparam D the number of Cartesian dimensions.
 */
template <size_t D> class FineNbrInfo : public NbrInfo<D>
{
	public:
	/**
	 * @brief The mpi rank that the neighbor resides on.
	 */
	std::array<int, Orthant<D - 1>::num_orthants> ranks;
	/**
	 * @brief The ids of the neighbors
	 */
	std::array<int, Orthant<D - 1>::num_orthants> ids;
	/**
	 * @brief The global indexes of the neighbors
	 */
	std::array<int, Orthant<D - 1>::num_orthants> global_indexes;
	/**
	 * @brief The local indexes of the neighbors
	 */
	std::array<int, Orthant<D - 1>::num_orthants> local_indexes;
	/**
	 * @brief Construct a new empty FineNbrInfo object
	 */
	FineNbrInfo()
	{
		ranks.fill(0);
		ids.fill(0);
		global_indexes.fill(0);
		local_indexes.fill(0);
	}
	~FineNbrInfo() = default;
	/**
	 * @brief Construct a new FineNbrInfo object
	 *
	 * @param ids the ids of the neighbors
	 */
	FineNbrInfo(std::array<int, Orthant<D - 1>::num_orthants> ids)
	{
		ranks.fill(0);
		this->ids = ids;
	}
	NbrType getNbrType()
	{
		return NbrType::Fine;
	}
	void getNbrIds(std::deque<int> &nbr_ids)
	{
		for (size_t i = 0; i < ids.size(); i++) {
			nbr_ids.push_back(ids[i]);
		}
	};
	void getNbrRanks(std::deque<int> &nbr_ranks)
	{
		for (size_t i = 0; i < ranks.size(); i++) {
			nbr_ranks.push_back(ranks[i]);
		}
	}
	void setGlobalIndexes(std::map<int, int> &rev_map)
	{
		for (size_t i = 0; i < global_indexes.size(); i++) {
			global_indexes[i] = rev_map.at(local_indexes[i]);
		}
	}
	void setLocalIndexes(std::map<int, int> &rev_map)
	{
		for (size_t i = 0; i < local_indexes.size(); i++) {
			local_indexes[i] = rev_map.at(ids[i]);
		}
	}
	int serialize(char *buffer) const
	{
		BufferWriter writer(buffer);
		writer << ranks;
		writer << ids;
		return writer.getPos();
	}
	int deserialize(char *buffer)
	{
		BufferReader reader(buffer);
		reader >> ranks;
		reader >> ids;
		return reader.getPos();
	}
};
} // namespace ThunderEgg
#endif