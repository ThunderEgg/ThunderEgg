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
template <int D> class CoarseNbrInfo : public NbrInfo<D>
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
	 * @brief The orthant that this patch in relation to the coarser patch's interface.
	 */
	Orthant<D> orth_on_coarse;
	/**
	 * @brief Construct a new empty CoarseNbrInfo object
	 */
	CoarseNbrInfo()  = default;
	~CoarseNbrInfo() = default;
	/**
	 * @brief Construct a new CoarseNbrInfo object
	 *
	 * @param id the id of the neighbor
	 * @param orth_on_coarse The orthant of the neighboring patch's interface that the this patch
	 * lies along.
	 */
	CoarseNbrInfo(int id, Orthant<D> orth_on_coarse)
	{
		this->id             = id;
		this->orth_on_coarse = orth_on_coarse;
	}
	NbrType getNbrType() const override
	{
		return NbrType::Coarse;
	}
	void getNbrIds(std::deque<int> &nbr_ids) const override
	{
		nbr_ids.push_back(id);
	};
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
		writer << orth_on_coarse;
		return writer.getPos();
	}
	int deserialize(char *buffer) override
	{
		BufferReader reader(buffer);
		reader >> rank;
		reader >> id;
		reader >> orth_on_coarse;
		return reader.getPos();
	}
	std::unique_ptr<NbrInfoBase> clone() const override
	{
		return std::make_unique<CoarseNbrInfo<D>>(*this);
	}
};
template <int D> void to_json(nlohmann::json &j, const CoarseNbrInfo<D> &n)
{
	j["type"]           = NbrType::Coarse;
	j["ids"]            = {n.id};
	j["ranks"]          = {n.rank};
	j["orth_on_coarse"] = n.orth_on_coarse;
}
template <int D> void from_json(const nlohmann::json &j, CoarseNbrInfo<D> &n)
{
	n.id   = j["ids"][0];
	n.rank = j["ranks"][0];
	if (j.contains("orth_on_coarse")) {
		n.orth_on_coarse = j["orth_on_coarse"].get<Orthant<D>>();
	} else {
		n.orth_on_coarse = Orthant<D>::null();
	}
}
} // namespace ThunderEgg
#endif