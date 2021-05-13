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
#include <ThunderEgg/Corner.h>
#include <ThunderEgg/Edge.h>
#include <ThunderEgg/FineNbrInfo.h>
#include <ThunderEgg/NormalNbrInfo.h>
#include <ThunderEgg/Orthant.h>
#include <ThunderEgg/RuntimeError.h>
#include <ThunderEgg/Serializable.h>
#include <ThunderEgg/TypeDefs.h>
#include <bitset>
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
	std::array<std::unique_ptr<NbrInfo<D>>, Side<D>::num_sides> nbr_info;
	/**
	 * @brief Nbr info objects for each edge
	 * If there is no neighbor, it should be set to nullptr.
	 */
	std::array<std::unique_ptr<NbrInfo<2>>, Edge<D>::num_edges> edge_nbr_info;
	/**
	 * @brief Nbr info objects for each corner
	 * If there is no neighbor, it should be set to nullptr.
	 */
	std::array<std::unique_ptr<NbrInfo<1>>, Corner<D>::num_corners> corner_nbr_info;

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
		for (Side<D> s : Side<D>::getValues()) {
			if (other_pinfo.hasNbr(s)) {
				switch (other_pinfo.getNbrType(s)) {
					case NbrType::Normal:
						nbr_info[s.getIndex()] = std::make_unique<NormalNbrInfo<D>>(other_pinfo.getNormalNbrInfo(s));
						break;
					case NbrType::Fine:
						nbr_info[s.getIndex()] = std::make_unique<FineNbrInfo<D>>(other_pinfo.getFineNbrInfo(s));
						break;
					case NbrType::Coarse:
						nbr_info[s.getIndex()] = std::make_unique<CoarseNbrInfo<D>>(other_pinfo.getCoarseNbrInfo(s));
						break;
					default:
						throw RuntimeError("Unknown NbrType");
						break;
				}
			}
		}
		for (Edge<D> e : Edge<D>::getValues()) {
			if (other_pinfo.hasNbr(e)) {
				switch (other_pinfo.getNbrType(e)) {
					case NbrType::Normal:
						edge_nbr_info[e.getIndex()] = std::make_unique<NormalNbrInfo<2>>(other_pinfo.getNormalNbrInfo(e));
						break;
					case NbrType::Fine:
						edge_nbr_info[e.getIndex()] = std::make_unique<FineNbrInfo<2>>(other_pinfo.getFineNbrInfo(e));
						break;
					case NbrType::Coarse:
						edge_nbr_info[e.getIndex()] = std::make_unique<CoarseNbrInfo<2>>(other_pinfo.getCoarseNbrInfo(e));
						break;
					default:
						throw RuntimeError("Unknown NbrType");
						break;
				}
			}
		}
		for (Corner<D> c : Corner<D>::getValues()) {
			if (other_pinfo.hasNbr(c)) {
				switch (other_pinfo.getNbrType(c)) {
					case NbrType::Normal:
						corner_nbr_info[c.getIndex()] = std::make_unique<NormalNbrInfo<1>>(other_pinfo.getNormalNbrInfo(c));
						break;
					case NbrType::Fine:
						corner_nbr_info[c.getIndex()] = std::make_unique<FineNbrInfo<1>>(other_pinfo.getFineNbrInfo(c));
						break;
					case NbrType::Coarse:
						corner_nbr_info[c.getIndex()] = std::make_unique<CoarseNbrInfo<1>>(other_pinfo.getCoarseNbrInfo(c));
						break;
					default:
						throw RuntimeError("Unknown NbrType");
						break;
				}
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
		for (Side<D> s : Side<D>::getValues()) {
			if (other_pinfo.hasNbr(s)) {
				switch (other_pinfo.getNbrType(s)) {
					case NbrType::Normal:
						nbr_info[s.getIndex()] = std::make_unique<NormalNbrInfo<D>>(other_pinfo.getNormalNbrInfo(s));
						break;
					case NbrType::Fine:
						nbr_info[s.getIndex()] = std::make_unique<FineNbrInfo<D>>(other_pinfo.getFineNbrInfo(s));
						break;
					case NbrType::Coarse:
						nbr_info[s.getIndex()] = std::make_unique<CoarseNbrInfo<D>>(other_pinfo.getCoarseNbrInfo(s));
						break;
					default:
						throw RuntimeError("Unknown NbrType");
						break;
				}
			} else {
				nbr_info[s.getIndex()] = nullptr;
			}
		}
		for (Edge<D> e : Edge<D>::getValues()) {
			if (other_pinfo.hasNbr(e)) {
				switch (other_pinfo.getNbrType(e)) {
					case NbrType::Normal:
						edge_nbr_info[e.getIndex()] = std::make_unique<NormalNbrInfo<2>>(other_pinfo.getNormalNbrInfo(e));
						break;
					case NbrType::Fine:
						edge_nbr_info[e.getIndex()] = std::make_unique<FineNbrInfo<2>>(other_pinfo.getFineNbrInfo(e));
						break;
					case NbrType::Coarse:
						edge_nbr_info[e.getIndex()] = std::make_unique<CoarseNbrInfo<2>>(other_pinfo.getCoarseNbrInfo(e));
						break;
					default:
						throw RuntimeError("Unknown NbrType");
						break;
				}
			}
		}
		for (Corner<D> c : Corner<D>::getValues()) {
			if (other_pinfo.hasNbr(c)) {
				switch (other_pinfo.getNbrType(c)) {
					case NbrType::Normal:
						corner_nbr_info[c.getIndex()] = std::make_unique<NormalNbrInfo<1>>(other_pinfo.getNormalNbrInfo(c));
						break;
					case NbrType::Fine:
						corner_nbr_info[c.getIndex()] = std::make_unique<FineNbrInfo<1>>(other_pinfo.getFineNbrInfo(c));
						break;
					case NbrType::Coarse:
						corner_nbr_info[c.getIndex()] = std::make_unique<CoarseNbrInfo<1>>(other_pinfo.getCoarseNbrInfo(c));
						break;
					default:
						throw RuntimeError("Unknown NbrType");
						break;
				}
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
		return *dynamic_cast<NormalNbrInfo<D> *>(nbr_info[s.getIndex()].get());
	}
	/**
	 * @brief Get the CoarseNbrInfo object
	 *

	 * @param s the side
	 * @return CoarseNbrInfo<D>& the object
	*/
	CoarseNbrInfo<D> &getCoarseNbrInfo(Side<D> s) const
	{
		return *dynamic_cast<CoarseNbrInfo<D> *>(nbr_info[s.getIndex()].get());
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
		return *dynamic_cast<FineNbrInfo<D> *>(nbr_info[s.getIndex()].get());
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
	 * @brief Get the NbrType for a corner
	 *
	 * @param c the corner
	 * @return The NbrType
	 */
	NbrType getNbrType(Corner<D> c) const
	{
		return corner_nbr_info[c.getIndex()]->getNbrType();
	}
	/**
	 * @brief Get the NormalNbrInfo object for a corner
	 *
	 * Neighbor must be of Normal type, otherwise behavior is undefined.
	 *
	 * @param c the corner
	 * @return NormalNbrInfo<1>& the object
	 */
	NormalNbrInfo<1> &getNormalNbrInfo(Corner<D> c) const
	{
		return *dynamic_cast<NormalNbrInfo<1> *>(corner_nbr_info[c.getIndex()].get());
	}
	/**
	 * @brief Get the CoarseNbrInfo object
	 *

	 * @param c the corner
	 * @return CoarseNbrInfo<1>& the object
	*/
	CoarseNbrInfo<1> &getCoarseNbrInfo(Corner<D> c) const
	{
		return *dynamic_cast<CoarseNbrInfo<1> *>(corner_nbr_info[c.getIndex()].get());
	}
	/**
	 * @brief Get the FineNbrInfo object
	 *
	 * Neighbor must be of Fine type, otherwise behavior is undefined.
	 *
	 * @param c the corner
	 * @return FineNbrInfo<1>& the object
	 */
	FineNbrInfo<1> &getFineNbrInfo(Corner<D> c) const
	{
		return *dynamic_cast<FineNbrInfo<1> *>(corner_nbr_info[c.getIndex()].get());
	}
	/**
	 * @brief Return whether the patch has a neighbor on it's corner
	 *
	 * @param c the corner
	 * @return true if the is neighbor on the corner
	 * @return false if not
	 */
	inline bool hasNbr(Corner<D> c) const
	{
		return corner_nbr_info[c.getIndex()] != nullptr;
	}
	/**
	 * @brief Get the NbrType for a edge
	 *
	 * @param e the edge
	 * @return The NbrType
	 */
	NbrType getNbrType(Edge<D> e) const
	{
		return edge_nbr_info[e.getIndex()]->getNbrType();
	}
	/**
	 * @brief Get the NormalNbrInfo object for a edge
	 *
	 * Neighbor must be of Normal type, otherwise behavior is undefined.
	 *
	 * @param e the edge
	 * @return NormalNbrInfo<2>& the object
	 */
	NormalNbrInfo<2> &getNormalNbrInfo(Edge<D> e) const
	{
		return *dynamic_cast<NormalNbrInfo<2> *>(edge_nbr_info[e.getIndex()].get());
	}
	/**
	 * @brief Get the CoarseNbrInfo object
	 *
	 * @param e the edge
	 * @return CoarseNbrInfo<2>& the object
	 */
	CoarseNbrInfo<2> &getCoarseNbrInfo(Edge<D> e) const
	{
		return *dynamic_cast<CoarseNbrInfo<2> *>(edge_nbr_info[e.getIndex()].get());
	}
	/**
	 * @brief Get the FineNbrInfo object
	 *
	 * Neighbor must be of Fine type, otherwise behavior is undefined.
	 *
	 * @param e the edge
	 * @return FineNbrInfo<2>& the object
	 */
	FineNbrInfo<2> &getFineNbrInfo(Edge<D> e) const
	{
		return *dynamic_cast<FineNbrInfo<2> *>(edge_nbr_info[e.getIndex()].get());
	}
	/**
	 * @brief Return whether the patch has a neighbor on it's edge
	 *
	 * @param e the edge
	 * @return true if the is neighbor on the edge
	 * @return false if not
	 */
	inline bool hasNbr(Edge<D> e) const
	{
		return edge_nbr_info[e.getIndex()] != nullptr;
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
		for (Side<D> s : Side<D>::getValues()) {
			if (hasNbr(s)) {
				nbr_info[s.getIndex()]->setLocalIndexes(id_to_local_index_map);
			}
		}
		for (Edge<D> e : Edge<D>::getValues()) {
			if (hasNbr(e)) {
				edge_nbr_info[e.getIndex()]->setLocalIndexes(id_to_local_index_map);
			}
		}
		for (Corner<D> c : Corner<D>::getValues()) {
			if (hasNbr(c)) {
				corner_nbr_info[c.getIndex()]->setLocalIndexes(id_to_local_index_map);
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
		for (Side<D> s : Side<D>::getValues()) {
			if (hasNbr(s)) {
				nbr_info[s.getIndex()]->setGlobalIndexes(id_to_global_index_map);
			}
		}
		for (Edge<D> e : Edge<D>::getValues()) {
			if (hasNbr(e)) {
				edge_nbr_info[e.getIndex()]->setGlobalIndexes(id_to_global_index_map);
			}
		}
		for (Corner<D> c : Corner<D>::getValues()) {
			if (hasNbr(c)) {
				corner_nbr_info[c.getIndex()]->setGlobalIndexes(id_to_global_index_map);
			}
		}
	}
	/**
	 * @brief return a vector of neighbor ids
	 */
	std::deque<int> getNbrIds() const
	{
		std::deque<int> retval;
		for (Side<D> s : Side<D>::getValues()) {
			if (hasNbr(s)) {
				nbr_info[s.getIndex()]->getNbrIds(retval);
			}
		}
		for (Edge<D> e : Edge<D>::getValues()) {
			if (hasNbr(e)) {
				edge_nbr_info[e.getIndex()]->getNbrIds(retval);
			}
		}
		for (Corner<D> c : Corner<D>::getValues()) {
			if (hasNbr(c)) {
				corner_nbr_info[c.getIndex()]->getNbrIds(retval);
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
		for (Side<D> s : Side<D>::getValues()) {
			if (hasNbr(s)) {
				nbr_info[s.getIndex()]->getNbrRanks(retval);
			}
		}
		for (Edge<D> e : Edge<D>::getValues()) {
			if (hasNbr(e)) {
				edge_nbr_info[e.getIndex()]->getNbrRanks(retval);
			}
		}
		for (Corner<D> c : Corner<D>::getValues()) {
			if (hasNbr(c)) {
				corner_nbr_info[c.getIndex()]->getNbrRanks(retval);
			}
		}
		return retval;
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
					default:
						throw RuntimeError("Unsupported NbrType");
				}
			}
		}
		std::bitset<Edge<D>::num_edges> has_edge_nbr;
		for (size_t i = 0; i < Edge<D>::num_edges; i++) {
			has_edge_nbr[i] = edge_nbr_info[i] != nullptr;
		}
		writer << has_edge_nbr;
		for (Edge<D> e : Edge<D>::getValues()) {
			if (hasNbr(e)) {
				NbrType type = getNbrType(e);
				writer << type;
				switch (type) {
					case NbrType::Normal: {
						NormalNbrInfo<2> tmp = getNormalNbrInfo(e);
						writer << tmp;
					} break;
					case NbrType::Fine: {
						FineNbrInfo<2> tmp = getFineNbrInfo(e);
						writer << tmp;
					} break;
					case NbrType::Coarse: {
						CoarseNbrInfo<2> tmp = getCoarseNbrInfo(e);
						writer << tmp;
					} break;
					default:
						throw RuntimeError("Unsupported NbrType");
				}
			}
		}
		std::bitset<Corner<D>::num_corners> has_corner_nbr;
		for (size_t i = 0; i < Corner<D>::num_corners; i++) {
			has_corner_nbr[i] = corner_nbr_info[i] != nullptr;
		}
		writer << has_corner_nbr;
		for (Corner<D> c : Corner<D>::getValues()) {
			if (hasNbr(c)) {
				NbrType type = getNbrType(c);
				writer << type;
				switch (type) {
					case NbrType::Normal: {
						NormalNbrInfo<1> tmp = getNormalNbrInfo(c);
						writer << tmp;
					} break;
					case NbrType::Fine: {
						FineNbrInfo<1> tmp = getFineNbrInfo(c);
						writer << tmp;
					} break;
					case NbrType::Coarse: {
						CoarseNbrInfo<1> tmp = getCoarseNbrInfo(c);
						writer << tmp;
					} break;
					default:
						throw RuntimeError("Unsupported NbrType");
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
		reader >> starts;
		reader >> spacings;
		std::bitset<Side<D>::num_sides> has_nbr;
		reader >> has_nbr;
		for (size_t i = 0; i < Side<D>::num_sides; i++) {
			if (has_nbr[i]) {
				NbrType type;
				reader >> type;
				NbrInfo<D> *info;
				switch (type) {
					case NbrType::Normal:
						info = new NormalNbrInfo<D>();
						reader >> *static_cast<NormalNbrInfo<D> *>(info);
						break;
					case NbrType::Fine:
						info = new FineNbrInfo<D>();
						reader >> *static_cast<FineNbrInfo<D> *>(info);
						break;
					case NbrType::Coarse:
						info = new CoarseNbrInfo<D>();
						reader >> *static_cast<CoarseNbrInfo<D> *>(info);
						break;
					default:
						throw RuntimeError("Unsupported NbrType");
				}
				nbr_info[i].reset(info);
			}
		}
		std::bitset<Edge<D>::num_edges> has_edge_nbr;
		reader >> has_edge_nbr;
		for (size_t i = 0; i < Edge<D>::num_edges; i++) {
			if (has_edge_nbr[i]) {
				NbrType type;
				reader >> type;
				NbrInfo<2> *info;
				switch (type) {
					case NbrType::Normal:
						info = new NormalNbrInfo<2>();
						reader >> *static_cast<NormalNbrInfo<2> *>(info);
						break;
					case NbrType::Fine:
						info = new FineNbrInfo<2>();
						reader >> *static_cast<FineNbrInfo<2> *>(info);
						break;
					case NbrType::Coarse:
						info = new CoarseNbrInfo<2>();
						reader >> *static_cast<CoarseNbrInfo<2> *>(info);
						break;
					default:
						throw RuntimeError("Unsupported NbrType");
				}
				edge_nbr_info[i].reset(info);
			}
		}
		std::bitset<Corner<D>::num_corners> has_corner_nbr;
		reader >> has_corner_nbr;
		for (size_t i = 0; i < Corner<D>::num_corners; i++) {
			if (has_corner_nbr[i]) {
				NbrType type;
				reader >> type;
				NbrInfo<1> *info;
				switch (type) {
					case NbrType::Normal:
						info = new NormalNbrInfo<1>();
						reader >> *static_cast<NormalNbrInfo<1> *>(info);
						break;
					case NbrType::Fine:
						info = new FineNbrInfo<1>();
						reader >> *static_cast<FineNbrInfo<1> *>(info);
						break;
					case NbrType::Coarse:
						info = new CoarseNbrInfo<1>();
						reader >> *static_cast<CoarseNbrInfo<1> *>(info);
						break;
					default:
						throw RuntimeError("Unsupported NbrType");
				}
				corner_nbr_info[i].reset(info);
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
	j["edge_nbrs"] = nlohmann::json::array();
	for (Edge<D> e : Edge<D>::getValues()) {
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
				pinfo.nbr_info[s.getIndex()] = std::make_unique<NormalNbrInfo<D>>();
				pinfo.getNormalNbrInfo(s)    = nbr_j.get<NormalNbrInfo<D>>();
				break;
			case NbrType::Coarse:
				pinfo.nbr_info[s.getIndex()] = std::make_unique<CoarseNbrInfo<D>>();
				pinfo.getCoarseNbrInfo(s)    = nbr_j.get<CoarseNbrInfo<D>>();
				break;
			case NbrType::Fine:
				pinfo.nbr_info[s.getIndex()] = std::make_unique<FineNbrInfo<D>>();
				pinfo.getFineNbrInfo(s)      = nbr_j.get<FineNbrInfo<D>>();
				break;
			default:
				throw RuntimeError("Unsupported NbrType");
		}
	}
	if (j.contains("edge_nbrs")) {
		for (const auto &nbr_j : j["edge_nbrs"]) {
			Edge<D> e = nbr_j["edge"].get<Edge<D>>();
			switch (nbr_j["type"].get<NbrType>()) {
				case NbrType::Normal:
					pinfo.edge_nbr_info[e.getIndex()] = std::make_unique<NormalNbrInfo<2>>();
					pinfo.getNormalNbrInfo(e)         = nbr_j.get<NormalNbrInfo<2>>();
					break;
				case NbrType::Coarse:
					pinfo.edge_nbr_info[e.getIndex()] = std::make_unique<CoarseNbrInfo<2>>();
					pinfo.getCoarseNbrInfo(e)         = nbr_j.get<CoarseNbrInfo<2>>();
					break;
				case NbrType::Fine:
					pinfo.edge_nbr_info[e.getIndex()] = std::make_unique<FineNbrInfo<2>>();
					pinfo.getFineNbrInfo(e)           = nbr_j.get<FineNbrInfo<2>>();
					break;
				default:
					throw RuntimeError("Unsupported NbrType");
			}
		}
	}
	if (j.contains("corner_nbrs")) {
		for (const auto &nbr_j : j["corner_nbrs"]) {
			Corner<D> c = nbr_j["corner"].get<Corner<D>>();
			switch (nbr_j["type"].get<NbrType>()) {
				case NbrType::Normal:
					pinfo.corner_nbr_info[c.getIndex()] = std::make_unique<NormalNbrInfo<1>>();
					pinfo.getNormalNbrInfo(c)           = nbr_j.get<NormalNbrInfo<1>>();
					break;
				case NbrType::Coarse:
					pinfo.corner_nbr_info[c.getIndex()] = std::make_unique<CoarseNbrInfo<1>>();
					pinfo.getCoarseNbrInfo(c)           = nbr_j.get<CoarseNbrInfo<1>>();
					break;
				case NbrType::Fine:
					pinfo.corner_nbr_info[c.getIndex()] = std::make_unique<FineNbrInfo<1>>();
					pinfo.getFineNbrInfo(c)             = nbr_j.get<FineNbrInfo<1>>();
					break;
				default:
					throw RuntimeError("Unsupported NbrType");
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