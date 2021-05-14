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
template <int D> class PatchInfo : public Serializable
{
	private:
	/**
	 * @brief Nbr info objects for each side.
	 * If there is no neighbor, it should be set to nullptr.
	 */
	std::array<std::unique_ptr<NbrInfoBase>, Face<D, D>::sum_of_faces> nbr_infos;

	public:
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
	 * @brief Construct a new Patch Info object
	 * starts, ns, and spacings are all set to 0
	 */
	PatchInfo()
	{
		starts.fill(0);
		ns.fill(1);
		spacings.fill(1);
		child_ids.fill(-1);
		child_ranks.fill(-1);
	}

	/**
	 * @brief Copy constructor
	 *
	 * @param other_pinfo  object to copy
	 */
	PatchInfo(const PatchInfo<D> &other_pinfo)
	: id(other_pinfo.id),
	  local_index(other_pinfo.local_index),
	  global_index(other_pinfo.global_index),
	  refine_level(other_pinfo.refine_level),
	  parent_id(other_pinfo.parent_id),
	  parent_rank(other_pinfo.parent_rank),
	  child_ids(other_pinfo.child_ids),
	  child_ranks(other_pinfo.child_ranks),
	  num_ghost_cells(other_pinfo.num_ghost_cells),
	  rank(other_pinfo.rank),
	  orth_on_parent(other_pinfo.orth_on_parent),
	  ns(other_pinfo.ns),
	  starts(other_pinfo.starts),
	  spacings(other_pinfo.spacings)
	{
		for (size_t i = 0; i < other_pinfo.nbr_infos.size(); i++) {
			if (other_pinfo.nbr_infos[i] != nullptr) {
				nbr_infos[i] = other_pinfo.nbr_infos[i]->clone();
			}
		}
	}
	/**
	 * @brief Copy asisgnment
	 *
	 * @param other_pinfo the object to copy
	 * @return PatchInfo<D>& this object
	 */
	PatchInfo<D> &operator=(const PatchInfo<D> &other_pinfo)
	{
		id              = other_pinfo.id;
		local_index     = other_pinfo.local_index;
		global_index    = other_pinfo.global_index;
		refine_level    = other_pinfo.refine_level;
		parent_id       = other_pinfo.parent_id;
		parent_rank     = other_pinfo.parent_rank;
		child_ids       = other_pinfo.child_ids;
		child_ranks     = other_pinfo.child_ranks;
		num_ghost_cells = other_pinfo.num_ghost_cells;
		rank            = other_pinfo.rank;
		orth_on_parent  = other_pinfo.orth_on_parent;
		ns              = other_pinfo.ns;
		starts          = other_pinfo.starts;
		spacings        = other_pinfo.spacings;

		for (size_t i = 0; i < other_pinfo.nbr_infos.size(); i++) {
			if (other_pinfo.nbr_infos[i] != nullptr) {
				nbr_infos[i] = other_pinfo.nbr_infos[i]->clone();
			} else {
				nbr_infos[i] = nullptr;
			}
		}
		return *this;
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
	template <int M> void setNbrInfo(Face<D, M> f, nullptr_t)
	{
		nbr_infos[Face<D, M>::sum_of_faces + f.getIndex()] = nullptr;
	}
	/**
	 * @brief Set the Nbr Info object on a given face
	 *
	 * @tparam M the dimensionality of the face
	 * @param s the face
	 * @param nbr_info the neighbor info object, this patchinfo will take ownership of it
	 */
	template <int M> void setNbrInfo(Face<D, M> f, NbrInfo<M> *nbr_info)
	{
		nbr_infos[Face<D, M>::sum_of_faces + f.getIndex()].reset(nbr_info);
	}
	/**
	 * @brief Get the NbrType for a side
	 *
	 * @param s the side
	 * @return The NbrType
	 */
	template <int M> NbrType getNbrType(Face<D, M> s) const
	{
		return nbr_infos[Face<D, M>::sum_of_faces + s.getIndex()]->getNbrType();
	}
	/**
	 * @brief Get the NormalNbrInfo object for a side
	 *
	 * Neighbor must be of Normal type, otherwise behavior is undefined.
	 *
	 * @param s the side
	 * @return NormalNbrInfo<D>& the object
	 */
	template <int M> NormalNbrInfo<M> &getNormalNbrInfo(Face<D, M> s) const
	{
		return *dynamic_cast<NormalNbrInfo<M> *>(nbr_infos[Face<D, M>::sum_of_faces + s.getIndex()].get());
	}
	/**
	 * @brief Get the CoarseNbrInfo object
	 *

	 * @param s the side
	 * @return CoarseNbrInfo<D>& the object
	*/
	template <int M> CoarseNbrInfo<M> &getCoarseNbrInfo(Face<D, M> s) const
	{
		return *dynamic_cast<CoarseNbrInfo<M> *>(nbr_infos[Face<D, M>::sum_of_faces + s.getIndex()].get());
	}
	/**
	 * @brief Get the FineNbrInfo object
	 *
	 * Neighbor must be of Fine type, otherwise behavior is undefined.
	 *
	 * @param s the side
	 * @return FineNbrInfo<D>& the object
	 */
	template <int M> FineNbrInfo<M> &getFineNbrInfo(Face<D, M> s) const
	{
		return *dynamic_cast<FineNbrInfo<M> *>(nbr_infos[Face<D, M>::sum_of_faces + s.getIndex()].get());
	}
	/**
	 * @brief Return whether the patch has a neighbor
	 *
	 * @param s the side
	 * @return true if the is neighbor
	 * @return false if at domain boundary
	 */
	template <int M> inline bool hasNbr(Face<D, M> s) const
	{
		return nbr_infos[Face<D, M>::sum_of_faces + s.getIndex()] != nullptr;
	}
	/**
	 * @brief Return if this patch has a neighbor
	 */
	inline bool hasNbr() const
	{
		for (const std::unique_ptr<NbrInfoBase> &nbr_info : nbr_infos) {
			if (nbr_info != nullptr) {
				return true;
			}
		}
		return false;
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
	 * @brief Set the local indexes in the NbrInfo objects
	 *
	 * @param id_to_local_index_map map from id to local_index
	 */
	void setNeighborLocalIndexes(const std::map<int, int> &id_to_local_index_map)
	{
		for (size_t i = 0; i < nbr_infos.size(); i++) {
			if (nbr_infos[i] != nullptr) {
				nbr_infos[i]->setLocalIndexes(id_to_local_index_map);
			}
		}
	}
	/**
	 * @brief Set the global indexes in the NbrInfo objects
	 *
	 * @param id_to_global_index_map map form id to global_index
	 */
	void setNeighborGlobalIndexes(const std::map<int, int> &id_to_global_index_map)
	{
		for (size_t i = 0; i < nbr_infos.size(); i++) {
			if (nbr_infos[i] != nullptr) {
				nbr_infos[i]->setGlobalIndexes(id_to_global_index_map);
			}
		}
	}
	/**
	 * @brief return a vector of neighbor ids
	 */
	std::deque<int> getNbrIds() const
	{
		std::deque<int> retval;
		for (size_t i = 0; i < nbr_infos.size(); i++) {
			if (nbr_infos[i] != nullptr) {
				nbr_infos[i]->getNbrIds(retval);
			}
		}
		return retval;
	}
	/**
	 * @brief return a vector of neighbor ranks
	 */
	std::deque<int> getNbrRanks() const
	{
		std::deque<int> retval;
		for (size_t i = 0; i < nbr_infos.size(); i++) {
			if (nbr_infos[i] != nullptr) {
				nbr_infos[i]->getNbrRanks(retval);
			}
		}
		return retval;
	}
	template <int M> void serializeNeighbors(BufferWriter &writer) const
	{
		std::bitset<Face<D, M>::number_of> has_nbr;
		for (Face<D, M> f : Face<D, M>::getValues()) {
			has_nbr[f.getIndex()] = hasNbr(f);
		}
		writer << has_nbr;
		for (Face<D, M> f : Face<D, M>::getValues()) {
			if (hasNbr(f)) {
				NbrType type = getNbrType(f);
				writer << type;
				switch (type) {
					case NbrType::Normal: {
						NormalNbrInfo<M> tmp = getNormalNbrInfo(f);
						writer << tmp;
					} break;
					case NbrType::Fine: {
						FineNbrInfo<M> tmp = getFineNbrInfo(f);
						writer << tmp;
					} break;
					case NbrType::Coarse: {
						CoarseNbrInfo<M> tmp = getCoarseNbrInfo(f);
						writer << tmp;
					} break;
					default:
						throw RuntimeError("Unsupported NbrType");
				}
			}
		}
		if constexpr (M > 0) {
			serializeNeighbors<M - 1>(writer);
		}
	}
	template <int M> void deserializeNeighbors(BufferReader &reader)
	{
		std::bitset<Face<D, M>::number_of> has_nbr;
		reader >> has_nbr;
		for (Face<D, M> f : Face<D, M>::getValues()) {
			if (has_nbr[f.getIndex()]) {
				NbrType type;
				reader >> type;
				NbrInfo<M> *info;
				switch (type) {
					case NbrType::Normal:
						info = new NormalNbrInfo<M>();
						reader >> *static_cast<NormalNbrInfo<M> *>(info);
						break;
					case NbrType::Fine:
						info = new FineNbrInfo<M>();
						reader >> *static_cast<FineNbrInfo<M> *>(info);
						break;
					case NbrType::Coarse:
						info = new CoarseNbrInfo<M>();
						reader >> *static_cast<CoarseNbrInfo<M> *>(info);
						break;
					default:
						throw RuntimeError("Unsupported NbrType");
				}

				setNbrInfo(f, info);
			}
		}
		if constexpr (M > 0) {
			deserializeNeighbors<M - 1>(reader);
		}
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
		writer << starts;
		writer << spacings;

		serializeNeighbors<D - 1>(writer);

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
		reader >> starts;
		reader >> spacings;

		deserializeNeighbors<D - 1>(reader);

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
	j["refine_level"]   = pinfo.refine_level;
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
	if constexpr (D == 3) {
		j["edge_nbrs"] = nlohmann::json::array();
		for (Edge e : Edge::getValues()) {
			if (pinfo.hasNbr(e)) {
				switch (pinfo.getNbrType(e)) {
					case NbrType::Normal:
						j["edge_nbrs"].push_back(pinfo.getNormalNbrInfo(e));
						break;
					case NbrType::Coarse:
						j["edge_nbrs"].push_back(pinfo.getCoarseNbrInfo(e));
						break;
					case NbrType::Fine:
						j["edge_nbrs"].push_back(pinfo.getFineNbrInfo(e));
						break;
					default:
						throw RuntimeError("Unsupported NbrType");
				}
				j["edge_nbrs"].back()["type"] = pinfo.getNbrType(e);
				j["edge_nbrs"].back()["edge"] = e;
			}
		}
	}
	if constexpr (D >= 2) {
		j["corner_nbrs"] = nlohmann::json::array();
		for (Corner<D> c : Corner<D>::getValues()) {
			if (pinfo.hasNbr(c)) {
				switch (pinfo.getNbrType(c)) {
					case NbrType::Normal:
						j["corner_nbrs"].push_back(pinfo.getNormalNbrInfo(c));
						break;
					case NbrType::Coarse:
						j["corner_nbrs"].push_back(pinfo.getCoarseNbrInfo(c));
						break;
					case NbrType::Fine:
						j["corner_nbrs"].push_back(pinfo.getFineNbrInfo(c));
						break;
					default:
						throw RuntimeError("Unsupported NbrType");
				}
				j["corner_nbrs"].back()["type"]   = pinfo.getNbrType(c);
				j["corner_nbrs"].back()["corner"] = c;
			}
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
	if (j.contains("refine_level")) {
		pinfo.refine_level = j["refine_level"];
	}
	pinfo.rank     = j["rank"];
	pinfo.starts   = j["starts"].get<std::array<double, D>>();
	pinfo.spacings = j["lengths"].get<std::array<double, D>>();
	pinfo.ns.fill(1);
	for (const auto &nbr_j : j["nbrs"]) {
		Side<D> s = nbr_j["side"].get<Side<D>>();
		switch (nbr_j["type"].get<NbrType>()) {
			case NbrType::Normal:
				pinfo.setNbrInfo(s, new NormalNbrInfo<D - 1>());
				pinfo.getNormalNbrInfo(s) = nbr_j.get<NormalNbrInfo<D - 1>>();
				break;
			case NbrType::Coarse:
				pinfo.setNbrInfo(s, new CoarseNbrInfo<D - 1>());
				pinfo.getCoarseNbrInfo(s) = nbr_j.get<CoarseNbrInfo<D - 1>>();
				break;
			case NbrType::Fine:
				pinfo.setNbrInfo(s, new FineNbrInfo<D - 1>());
				pinfo.getFineNbrInfo(s) = nbr_j.get<FineNbrInfo<D - 1>>();
				break;
			default:
				throw RuntimeError("Unsupported NbrType");
		}
	}
	if constexpr (D == 3) {
		if (j.contains("edge_nbrs")) {
			for (const auto &nbr_j : j["edge_nbrs"]) {
				Edge e = nbr_j["edge"].get<Edge>();
				switch (nbr_j["type"].get<NbrType>()) {
					case NbrType::Normal:
						pinfo.setNbrInfo(e, new NormalNbrInfo<1>());
						pinfo.getNormalNbrInfo(e) = nbr_j.get<NormalNbrInfo<1>>();
						break;
					case NbrType::Coarse:
						pinfo.setNbrInfo(e, new CoarseNbrInfo<1>());
						pinfo.getCoarseNbrInfo(e) = nbr_j.get<CoarseNbrInfo<1>>();
						break;
					case NbrType::Fine:
						pinfo.setNbrInfo(e, new FineNbrInfo<1>());
						pinfo.getFineNbrInfo(e) = nbr_j.get<FineNbrInfo<1>>();
						break;
					default:
						throw RuntimeError("Unsupported NbrType");
				}
			}
		}
	}
	if constexpr (D >= 2) {
		if (j.contains("corner_nbrs")) {
			for (const auto &nbr_j : j["corner_nbrs"]) {
				Corner<D> c = nbr_j["corner"].get<Corner<D>>();
				switch (nbr_j["type"].get<NbrType>()) {
					case NbrType::Normal:
						pinfo.setNbrInfo(c, new NormalNbrInfo<0>());
						pinfo.getNormalNbrInfo(c) = nbr_j.get<NormalNbrInfo<0>>();
						break;
					case NbrType::Coarse:
						pinfo.setNbrInfo(c, new CoarseNbrInfo<0>());
						pinfo.getCoarseNbrInfo(c) = nbr_j.get<CoarseNbrInfo<0>>();
						break;
					case NbrType::Fine:
						pinfo.setNbrInfo(c, new FineNbrInfo<0>());
						pinfo.getFineNbrInfo(c) = nbr_j.get<FineNbrInfo<0>>();
						break;
					default:
						throw RuntimeError("Unsupported NbrType");
				}
			}
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