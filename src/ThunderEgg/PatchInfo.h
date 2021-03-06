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
#ifndef THUNDEREGG_PATCHINFO_H
#define THUNDEREGG_PATCHINFO_H
#include <ThunderEgg/CoarseNbrInfo.h>
#include <ThunderEgg/FineNbrInfo.h>
#include <ThunderEgg/NormalNbrInfo.h>
#include <ThunderEgg/Orthant.h>
#include <ThunderEgg/RuntimeError.h>
#include <ThunderEgg/Serializable.h>
#include <ThunderEgg/TypeDefs.h>
#include <bitset>
#include <deque>
#include <map>
#include <memory>

namespace ThunderEgg
{
/**
 * @brief Contains metadata for a patch
 *
 * This contains metadata for a specific patch. Information like:
 * * The globally unique id of this patch
 * * The local and global indexes of this patch in the Domain
 * * The parent patch in the tree (if there is one)
 *
 * It also contains information for a patch's neighbor:
 * * What are the neighbors id?
 * * Are the neighbors at the same refinement level? Are coarser or finer?
 *
 * @tparam D the number of cartesian dimensions in the patch
 */
template <int D> struct PatchInfo : public Serializable {
	/**
	 * @brief The globally unique ID of the patch
	 * This ID only needs to be unique within a Domain.
	 */
	int id = 0;
	/**
	 * @brief The local index of the patch in the Domain.
	 */
	int local_index = 0;
	/**
	 * @brief The global index of the patch in the Domain.
	 */
	int global_index = 0;
	/**
	 * @brief The refinement level
	 */
	int refine_level = -1;
	/**
	 * @brief The id of the parent patch.
	 *
	 * Set to -1 if there is no parent.
	 */
	int parent_id = -1;
	/**
	 * @brief the rank that the parent patch resides on
	 */
	int parent_rank = -1;
	/**
	 * @brief The id's of the children.
	 *
	 * Set to -1 if there are no children
	 */
	std::array<int, Orthant<D>::num_orthants> child_ids;
	/**
	 * @brief The ranks of the children
	 *
	 * Set to -1 if there are no children
	 */
	std::array<int, Orthant<D>::num_orthants> child_ranks;
	/**
	 * @brief Number of ghost cells on each side of the patch.
	 */
	int num_ghost_cells = 0;
	/**
	 * @brief MPI rank of this patch
	 */
	int rank = -1;
	/**
	 * @brief The orthant of the parent that this parent resides on.
	 *
	 * If the parent is the same size, it should be set to Orthant::null
	 */
	Orthant<D> orth_on_parent = Orthant<D>::null();
	/**
	 * @brief Whether the patch has neumann boundary conditions on one side.
	 */
	std::bitset<Side<D>::num_sides> neumann;
	/**
	 * @brief The number of cells in each direction
	 */
	std::array<int, D> ns;
	/**
	 * @brief The lower-left-bottom index of the patch
	 */
	std::array<double, D> starts;
	/**
	 * @brief The cell spacings in each direction
	 */
	std::array<double, D> spacings;
	/**
	 * @brief Nbr info objects for each side.
	 * If there is no neighbor, it should be set to nullptr.
	 */
	std::array<std::shared_ptr<NbrInfo<D>>, Side<D>::num_sides> nbr_info;
	/**
	 * @brief local index in the boundary conditions vector
	 */
	std::array<int, Side<D>::num_sides> bc_local_index;
	/**
	 * @brief global index in the boundary conditions vector
	 */
	std::array<int, Side<D>::num_sides> bc_global_index;

	/**
	 * @brief Construct a new Patch Info object
	 * starts, ns, and spacings are all set to 0
	 */
	PatchInfo()
	{
		starts.fill(0);
		nbr_info.fill(nullptr);
		ns.fill(1);
		spacings.fill(1);
		bc_local_index.fill(-1);
		bc_global_index.fill(-1);
		child_ids.fill(-1);
		child_ranks.fill(-1);
	}
	/**
	 * @brief Compare the ids of the patches
	 *
	 * @param l left operand
	 * @param r right operand
	 * @return true if r's id is lower
	 * @return false if r's id is not lower
	 */
	friend bool operator<(const PatchInfo &l, const PatchInfo &r)
	{
		return l.id < r.id;
	}
	/**
	 * @brief Get the NbrType for a side
	 *
	 * @param s the side
	 * @return The NbrType
	 */
	NbrType getNbrType(Side<D> s) const
	{
		return nbr_info[s.getIndex()]->getNbrType();
	}
	/**
	 * @brief Get the NormalNbrInfo object for a side
	 *
	 * Neighbor must be of Normal type, otherwise behavior is undefined.
	 *
	 * @param s the side
	 * @return NormalNbrInfo<D>& the object
	 */
	NormalNbrInfo<D> &getNormalNbrInfo(Side<D> s) const
	{
		return *std::dynamic_pointer_cast<NormalNbrInfo<D>>(nbr_info[s.getIndex()]);
	}
	/**
	 * @brief Get the NormalNbrInfo pointer
	 *
	 * Neighbor must be of Normal type, otherwise a nullptr will be returned.
	 *
	 * @param s the side
	 * @return std::shared_ptr<NormalNbrInfo<D>> the pointer
	 */
	std::shared_ptr<NormalNbrInfo<D>> getNormalNbrInfoPtr(Side<D> s) const
	{
		return std::dynamic_pointer_cast<NormalNbrInfo<D>>(nbr_info[s.getIndex()]);
	}
	/**
	 * @brief Get the CoarseNbrInfo object
	 *

	 * @param s the side
	 * @return CoarseNbrInfo<D>& the object
	*/
	CoarseNbrInfo<D> &getCoarseNbrInfo(Side<D> s) const
	{
		return *std::dynamic_pointer_cast<CoarseNbrInfo<D>>(nbr_info[s.getIndex()]);
	}
	/**
	 * @brief Get the CoarseNbrInfo pointer
	 *
	 * Neighbor must be of Coarse type, otherwise a nullptr will be returned.
	 *
	 * @param s the side
	 * @return std::shared_ptr<CoarseNbrInfo<D>> the pointer
	 */
	std::shared_ptr<CoarseNbrInfo<D>> getCoarseNbrInfoPtr(Side<D> s) const
	{
		return std::dynamic_pointer_cast<CoarseNbrInfo<D>>(nbr_info[s.getIndex()]);
	}
	/**
	 * @brief Get the FineNbrInfo object
	 *
	 * Neighbor must be of Fine type, otherwise behavior is undefined.
	 *
	 * @param s the side
	 * @return FineNbrInfo<D>& the object
	 */
	FineNbrInfo<D> &getFineNbrInfo(Side<D> s) const
	{
		return *std::dynamic_pointer_cast<FineNbrInfo<D>>(nbr_info[s.getIndex()]);
	}
	/**
	 * @brief Get the FineNbrInfo pointer
	 *
	 * Neighbor must be of Coarse type, otherwise a nullptr will be returned.
	 *
	 * @param s the side
	 * @return std::shared_ptr<FineNbrInfo<D>> the pointer
	 */
	std::shared_ptr<FineNbrInfo<D>> getFineNbrInfoPtr(Side<D> s) const
	{
		return std::dynamic_pointer_cast<FineNbrInfo<D>>(nbr_info[s.getIndex()]);
	}
	/**
	 * @brief Return whether the patch has a neighbor
	 *
	 * @param s the side
	 * @return true if the is neighbor
	 * @return false if at domain boundary
	 */
	inline bool hasNbr(Side<D> s) const
	{
		return nbr_info[s.getIndex()] != nullptr;
	}
	/**
	 * @brief Return whether the patch has a coarser parent
	 *
	 * @return true if there is a parent
	 */
	inline bool hasCoarseParent() const
	{
		return orth_on_parent != Orthant<D>::null();
	}
	/**
	 * @brief Return wether the boundary conditions are neumann
	 *
	 * @param s the side
	 */
	inline bool isNeumann(Side<D> s) const
	{
		return neumann[s.getIndex()];
	}
	/**
	 * @brief Set the local indexes in the NbrInfo objects
	 *
	 * @param rev_map map from id to local_index
	 */
	void setLocalNeighborIndexes(std::map<int, int> &rev_map)
	{
		local_index = rev_map.at(id);
		for (Side<D> s : Side<D>::getValues()) {
			if (hasNbr(s)) {
				nbr_info[s.getIndex()]->setLocalIndexes(rev_map);
			}
		}
	}
	/**
	 * @brief Set the global indexes in the NbrInfo objects
	 *
	 * @param rev_map map form local_index to global_index
	 */
	void setGlobalNeighborIndexes(std::map<int, int> &rev_map)
	{
		global_index = rev_map.at(local_index);
		for (Side<D> s : Side<D>::getValues()) {
			if (hasNbr(s)) {
				nbr_info[s.getIndex()]->setGlobalIndexes(rev_map);
			}
		}
	}
	/**
	 * @brief Set the neumann boundary conditions
	 *
	 * @param inf the function for determining boundary conditions
	 */
	void setNeumann(IsNeumannFunc<D> inf)
	{
		for (Side<D> s : Side<D>::getValues()) {
			if (!hasNbr(s)) {
				std::array<double, D> bound_start = starts;
				if (!s.isLowerOnAxis()) {
					bound_start[s.getAxisIndex()]
					+= spacings[s.getAxisIndex()] * ns[s.getAxisIndex()];
				}
				std::array<double, D> bound_end = bound_start;
				for (size_t dir = 0; dir < D; dir++) {
					if (dir != s.getAxisIndex()) {
						bound_end[dir] += spacings[dir] * ns[dir];
					}
				}
				neumann[s.getIndex()] = inf(s, bound_end, bound_start);
			}
		}
	}
	/**
	 * @brief return a vector of neighbor ids
	 */
	std::deque<int> getNbrIds()
	{
		std::deque<int> retval;
		for (Side<D> s : Side<D>::getValues()) {
			if (hasNbr(s)) {
				nbr_info[s.getIndex()]->getNbrIds(retval);
			}
		}
		return retval;
	}
	/**
	 * @brief return a vector of neighbor ranks
	 */
	std::deque<int> getNbrRanks()
	{
		std::deque<int> retval;
		for (Side<D> s : Side<D>::getValues()) {
			if (hasNbr(s)) {
				nbr_info[s.getIndex()]->getNbrRanks(retval);
			}
		}
		return retval;
	}
	/**
	 * @brief set the local index in the boundary condition vector
	 *
	 * @param s the side that the boundary is on
	 * @param local_index the index
	 */
	void setBCLocalIndex(Side<D> s, int local_index)
	{
		bc_local_index[s.getIndex()] = local_index;
	}
	/**
	 * @brief get the local index in the boundnary condition vector
	 *
	 * @param s the side that the boundary is on
	 * @return int the index
	 */
	int getBCLocalIndex(Side<D> s)
	{
		return bc_local_index[s.getIndex()];
	}
	/**
	 * @brief set the global index in the boundary condition vector
	 *
	 * @param s the side that the boundary is on
	 * @param local_index the index
	 */
	void setBCGlobalIndex(Side<D> s, int global_index)
	{
		bc_global_index[s.getIndex()] = global_index;
	}
	/**
	 * @brief get the global index in the boundnary condition vector
	 *
	 * @param s the side that the boundary is on
	 * @return int the index
	 */
	int getBCGlobalIndex(Side<D> s)
	{
		return bc_global_index[s.getIndex()];
	}
	int serialize(char *buffer) const
	{
		BufferWriter writer(buffer);
		writer << id;
		writer << ns;
		writer << refine_level;
		writer << parent_id;
		writer << parent_rank;
		writer << child_ids;
		writer << child_ranks;
		writer << rank;
		writer << orth_on_parent;
		writer << neumann;
		writer << starts;
		writer << spacings;
		std::bitset<Side<D>::num_sides> has_nbr;
		for (size_t i = 0; i < Side<D>::num_sides; i++) {
			has_nbr[i] = nbr_info[i] != nullptr;
		}
		writer << has_nbr;
		for (Side<D> s : Side<D>::getValues()) {
			if (hasNbr(s)) {
				NbrType type = getNbrType(s);
				writer << type;
				switch (type) {
					case NbrType::Normal: {
						NormalNbrInfo<D> tmp = getNormalNbrInfo(s);
						writer << tmp;
					} break;
					case NbrType::Fine: {
						FineNbrInfo<D> tmp = getFineNbrInfo(s);
						writer << tmp;
					} break;
					case NbrType::Coarse: {
						CoarseNbrInfo<D> tmp = getCoarseNbrInfo(s);
						writer << tmp;
					} break;
				}
			}
		}
		return writer.getPos();
	}
	int deserialize(char *buffer)
	{
		BufferReader reader(buffer);
		reader >> id;
		reader >> ns;
		reader >> refine_level;
		reader >> parent_id;
		reader >> parent_rank;
		reader >> child_ids;
		reader >> child_ranks;
		reader >> rank;
		reader >> orth_on_parent;
		reader >> neumann;
		reader >> starts;
		reader >> spacings;
		std::bitset<Side<D>::num_sides> has_nbr;
		reader >> has_nbr;
		for (size_t i = 0; i < Side<D>::num_sides; i++) {
			if (has_nbr[i]) {
				NbrType type;
				reader >> type;
				std::shared_ptr<NbrInfo<D>> info;
				switch (type) {
					case NbrType::Normal:
						info.reset(new NormalNbrInfo<D>());
						reader >> *std::static_pointer_cast<NormalNbrInfo<D>>(info);
						break;
					case NbrType::Fine:
						info.reset(new FineNbrInfo<D>());
						reader >> *std::static_pointer_cast<FineNbrInfo<D>>(info);
						break;
					case NbrType::Coarse:
						info.reset(new CoarseNbrInfo<D>());
						reader >> *std::static_pointer_cast<CoarseNbrInfo<D>>(info);
						break;
				}
				nbr_info[i] = info;
			}
		}
		return reader.getPos();
	}
};
template <int D> void to_json(nlohmann::json &j, const PatchInfo<D> &pinfo)
{
	j["id"]             = pinfo.id;
	j["parent_id"]      = pinfo.parent_id;
	j["parent_rank"]    = pinfo.parent_rank;
	j["orth_on_parent"] = pinfo.orth_on_parent;
	j["rank"]           = pinfo.rank;
	j["starts"]         = pinfo.starts;
	j["lengths"]        = pinfo.spacings;
	for (int i = 0; i < D; i++) {
		j["lengths"][i] = pinfo.spacings[i] * pinfo.ns[i];
	}
	j["nbrs"] = nlohmann::json::array();
	for (Side<D> s : Side<D>::getValues()) {
		if (pinfo.hasNbr(s)) {
			switch (pinfo.getNbrType(s)) {
				case NbrType::Normal:
					j["nbrs"].push_back(pinfo.getNormalNbrInfo(s));
					break;
				case NbrType::Coarse:
					j["nbrs"].push_back(pinfo.getCoarseNbrInfo(s));
					break;
				case NbrType::Fine:
					j["nbrs"].push_back(pinfo.getFineNbrInfo(s));
					break;
				default:
					throw RuntimeError("Unsupported NbrType");
			}
			j["nbrs"].back()["type"] = pinfo.getNbrType(s);
			j["nbrs"].back()["side"] = s;
		}
	}
	if (pinfo.child_ids[0] != -1) {
		j["child_ids"]   = pinfo.child_ids;
		j["child_ranks"] = pinfo.child_ranks;
	}
}
template <int D> void from_json(const nlohmann::json &j, PatchInfo<D> &pinfo)
{
	pinfo.id          = j["id"];
	pinfo.parent_id   = j["parent_id"];
	pinfo.parent_rank = j["parent_rank"];
	if (j.contains("orth_on_parent")) {
		j["orth_on_parent"].get_to(pinfo.orth_on_parent);
	}
	pinfo.rank     = j["rank"];
	pinfo.starts   = j["starts"].get<std::array<double, D>>();
	pinfo.spacings = j["lengths"].get<std::array<double, D>>();
	pinfo.ns.fill(1);
	for (const auto &nbr_j : j["nbrs"]) {
		Side<D> s = nbr_j["side"].get<Side<D>>();
		switch (nbr_j["type"].get<NbrType>()) {
			case NbrType::Normal:
				pinfo.nbr_info[s.getIndex()] = std::make_shared<NormalNbrInfo<D>>();
				pinfo.getNormalNbrInfo(s)    = nbr_j.get<NormalNbrInfo<D>>();
				break;
			case NbrType::Coarse:
				pinfo.nbr_info[s.getIndex()] = std::make_shared<CoarseNbrInfo<D>>();
				pinfo.getCoarseNbrInfo(s)    = nbr_j.get<CoarseNbrInfo<D>>();
				break;
			case NbrType::Fine:
				pinfo.nbr_info[s.getIndex()] = std::make_shared<FineNbrInfo<D>>();
				pinfo.getFineNbrInfo(s)      = nbr_j.get<FineNbrInfo<D>>();
				break;
			default:
				throw RuntimeError("Unsupported NbrType");
		}
	}
	if (j.contains("child_ids")) {
		pinfo.child_ids   = j["child_ids"].get<std::array<int, Orthant<D>::num_orthants>>();
		pinfo.child_ranks = j["child_ranks"].get<std::array<int, Orthant<D>::num_orthants>>();
	}
}
extern template struct PatchInfo<2>;
extern template struct PatchInfo<3>;
} // namespace ThunderEgg
#endif