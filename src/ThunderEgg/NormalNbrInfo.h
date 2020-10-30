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
#ifndef THUNDEREGG_NORMALNBRINFO_H
#define THUNDEREGG_NORMALNBRINFO_H
#include <ThunderEgg/BufferReader.h>
#include <ThunderEgg/BufferWriter.h>
#include <ThunderEgg/NbrInfo.h>

namespace ThunderEgg
{
/**
 * @brief Represents a neighbor that is at the same refinement level.
 *
 * @tparam D the number of Cartesian dimensions on a patch.
 */
template <int D> class NormalNbrInfo : public NbrInfo<D>
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
	 * @brief Construct a new empty NormalNbrInfo object
	 */
	NormalNbrInfo() {}
	~NormalNbrInfo() = default;
	/**
	 * @brief Construct a new NormalNbrInfo object
	 *
	 * @param id the id of the neighbor.
	 */
	NormalNbrInfo(int id)
	{
		this->id = id;
	}
	NbrType getNbrType()
	{
		return NbrType::Normal;
	}
	void getNbrIds(std::deque<int> &nbr_ids)
	{
		nbr_ids.push_back(id);
	}
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
		return writer.getPos();
	}
	int deserialize(char *buffer)
	{
		BufferReader reader(buffer);
		reader >> rank;
		reader >> id;
		return reader.getPos();
	}
};
template <int D> void to_json(nlohmann::json &j, const NormalNbrInfo<D> &n)
{
	j["type"]  = NbrType::Normal;
	j["ids"]   = {n.id};
	j["ranks"] = {n.rank};
}
template <int D> void from_json(const nlohmann::json &j, NormalNbrInfo<D> &n)
{
	n.id   = j["ids"][0];
	n.rank = j["ranks"][0];
}
} // namespace ThunderEgg
#endif