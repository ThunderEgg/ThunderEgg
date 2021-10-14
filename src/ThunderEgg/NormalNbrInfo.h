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
#ifndef THUNDEREGG_NORMALNBRINFO_H
#define THUNDEREGG_NORMALNBRINFO_H
/**
 * @file
 *
 * @brief NormalNbrType class
 */
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
	int local_index = -1;
	/**
	 * @brief The global index of the neighbor
	 */
	int global_index = -1;
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
	NbrType getNbrType() const override
	{
		return NbrType::Normal;
	}
	void getNbrIds(std::deque<int> &nbr_ids) const override
	{
		nbr_ids.push_back(id);
	}
	void getNbrRanks(std::deque<int> &nbr_ranks) const override
	{
		nbr_ranks.push_back(rank);
	}
	void setGlobalIndexes(const std::map<int, int> &id_to_global_index_map) override
	{
		global_index = id_to_global_index_map.at(id);
	}
	void setLocalIndexes(const std::map<int, int> &id_to_local_index_map) override
	{
		auto iter = id_to_local_index_map.find(id);
		if (iter != id_to_local_index_map.end()) {
			local_index = iter->second;
		}
	}
	int serialize(char *buffer) const override
	{
		BufferWriter writer(buffer);
		writer << rank;
		writer << id;
		return writer.getPos();
	}
	int deserialize(char *buffer) override
	{
		BufferReader reader(buffer);
		reader >> rank;
		reader >> id;
		return reader.getPos();
	}
	std::unique_ptr<NbrInfoBase> clone() const override
	{
		return std::make_unique<NormalNbrInfo<D>>(*this);
	}
};
template <int D> void to_json(tpl::nlohmann::json &j, const NormalNbrInfo<D> &n)
{
	j["type"]  = NbrType::Normal;
	j["ids"]   = {n.id};
	j["ranks"] = {n.rank};
}
template <int D> void from_json(const tpl::nlohmann::json &j, NormalNbrInfo<D> &n)
{
	n.id   = j["ids"][0];
	n.rank = j["ranks"][0];
}
} // namespace ThunderEgg
#endif