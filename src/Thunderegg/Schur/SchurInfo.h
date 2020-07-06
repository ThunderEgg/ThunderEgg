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

#ifndef THUNDEREGG_SCHUR_SCHURINFO_H
#define THUNDEREGG_SCHUR_SCHURINFO_H
#include <Thunderegg/PatchInfo.h>
#include <Thunderegg/Schur/IfaceType.h>
#include <deque>
#include <set>
namespace Thunderegg
{
namespace Schur
{
/**
 * @brief The IfaceInfo class represents the information for an interface on a given side of the
 * patch.
 *
 * The information contained  will be the globally unique ID and the local and global index(es) in
 * the interface vector.
 *
 * @tparam D the number of Cartesian dimensions in the patches.
 */
template <size_t D> class IfaceInfo
{
	public:
	/**
	 * @brief The rank that the interface resides on.
	 */
	int rank;
	/**
	 * @brief The globally unique ID of the interface.
	 */
	int id;
	/**
	 * @brief the local index in the interface vector.
	 */
	int local_index;
	/**
	 * @brief the global index in the interface vector.
	 */
	int global_index;
	/**
	 * @brief add to a deque of globally unique ids
	 *
	 * @param ids adds to the deque of ids
	 */
	virtual void getIds(std::deque<int> &ids) = 0;
	/**
	 * @brief add to a deque of local interface indexes
	 *
	 * @param idx adds to the deque of local interface indexes
	 */
	virtual void getLocalIndexes(std::deque<int> &idx) = 0;
	/**
	 * @brief add to a deque of global interface indexes
	 *
	 * @param idx adds to the deque of global interface indexes
	 */
	virtual void getGlobalIndexes(std::deque<int> &idx) = 0;
	/**
	 * @brief add to a deque of IfaceTypes
	 *
	 * @param types adds to the deque of IfaceTypes
	 */
	virtual void getIfaceTypes(std::deque<IfaceType<D>> &types) = 0;
	/**
	 * @brief add to a deque of interface ranks
	 *
	 * @param ranks adds to the deque of interface ranks
	 */
	virtual void getRanks(std::deque<int> &ranks) = 0;
	/**
	 * @brief Set the local indexes in the IfaceInfo objects
	 *
	 * @param rev_map map from id to local_index
	 */
	virtual void setLocalIndexes(const std::map<int, int> &rev_map) = 0;
	/**
	 * @brief Set the global indexes in the IfaceInfo objects
	 *
	 * @param rev_map map form local_index to global_index
	 */
	virtual void setGlobalIndexes(const std::map<int, int> &rev_map) = 0;
};
/**
 * @brief This represents an interface where the neighbor is at the same refinement level
 *
 * @tparam D the number of Cartesian dimensions in a patch
 */
template <size_t D> class NormalIfaceInfo : public IfaceInfo<D>
{
	public:
	/**
	 * @brief convenience pointer to associated NbrInfo object
	 */
	std::shared_ptr<NormalNbrInfo<D>> nbr_info;
	/**
	 * @brief Construct a new NormalIfaceInfo object all values are set to zero
	 */
	NormalIfaceInfo()
	{
		this->id           = 0;
		this->local_index  = 0;
		this->global_index = 0;
	}
	/**
	 * @brief Construct a new NormalIfaceInfo object
	 *
	 * @param pinfo the associated PatchInfo object
	 * @param s the side of the patch that the interface is on
	 */
	NormalIfaceInfo(std::shared_ptr<PatchInfo<D>> pinfo, Side<D> s)
	{
		nbr_info = pinfo->getNormalNbrInfoPtr(s);
		if (s.isLowerOnAxis()) {
			this->id = pinfo->id * Side<D>::num_sides + s.getIndex();
			// lower axis interface belongs to neighboring rank
			this->rank = nbr_info->rank;
		} else {
			this->id = nbr_info->id * Side<D>::num_sides + s.opposite().getIndex();
			// higher axis interafce belongs to this patch's rank
			this->rank = pinfo->rank;
		}
	}
	void getIds(std::deque<int> &ids)
	{
		ids.push_back(this->id);
	}
	void getLocalIndexes(std::deque<int> &idx)
	{
		idx.push_back(this->local_index);
	}
	void getGlobalIndexes(std::deque<int> &idx)
	{
		idx.push_back(this->global_index);
	}
	void getIfaceTypes(std::deque<IfaceType<D>> &types)
	{
		types.push_back(IfaceType<D>::normal);
	}
	void getRanks(std::deque<int> &ranks)
	{
		ranks.push_back(nbr_info->rank);
	}
	void setLocalIndexes(const std::map<int, int> &rev_map)
	{
		this->local_index = rev_map.at(this->id);
	}
	void setGlobalIndexes(const std::map<int, int> &rev_map)
	{
		this->global_index = rev_map.at(this->local_index);
	}
};
/**
 * @brief Represents the interfaces where the neighbor is at a coarser refinement level.
 *
 * There will be two interfaces associated with this object. The interface that lines up with this
 * patch, and interface that lines up with the coarser patch.
 *
 * @tparam D the number of Cartesian dimensions in a patch
 */
template <size_t D> class CoarseIfaceInfo : public IfaceInfo<D>
{
	public:
	/**
	 * @brief convenience pointer to associated NbrInfo object
	 */
	std::shared_ptr<CoarseNbrInfo<D>> nbr_info;
	/**
	 * @brief The orthant that this patch in relation to the coarser patch's interface.
	 */
	Orthant<D - 1> orth_on_coarse;
	/**
	 * @brief Rank of the coarse interface
	 */
	int coarse_rank;
	/**
	 * @brief The id of the coarser patch's interface
	 */
	int coarse_id;
	/**
	 * @brief The local index of the coarser patch's inteface
	 */
	int coarse_local_index;
	/**
	 * @brief The global index of the coarser patch's interface
	 */
	int coarse_global_index;
	/**
	 * @brief Construct a new CoarseIfaceInfo object
	 *
	 * @param pinfo the cooresponding PatchInfo object
	 * @param s the side that the interface is on
	 */
	CoarseIfaceInfo(std::shared_ptr<PatchInfo<D>> pinfo, Side<D> s)
	{
		nbr_info       = pinfo->getCoarseNbrInfoPtr(s);
		this->id       = pinfo->id * Side<D>::num_sides + s.getIndex();
		orth_on_coarse = nbr_info->orth_on_coarse;
		coarse_id      = nbr_info->id * Side<D>::num_sides + s.opposite().getIndex();
		// fine and coarse interfaces always belong to their patches
		this->rank        = pinfo->rank;
		this->coarse_rank = nbr_info->rank;
	}
	void getIds(std::deque<int> &ids)
	{
		ids.push_back(this->id);
		ids.push_back(coarse_id);
	}
	void getLocalIndexes(std::deque<int> &idx)
	{
		idx.push_back(this->local_index);
		idx.push_back(coarse_local_index);
	}
	void getGlobalIndexes(std::deque<int> &idx)
	{
		idx.push_back(this->global_index);
		idx.push_back(coarse_global_index);
	}
	void getIfaceTypes(std::deque<IfaceType<D>> &types)
	{
		IfaceType<D> fine_type(IfaceType<D>::fine_to_fine, orth_on_coarse);
		IfaceType<D> coarse_type(IfaceType<D>::fine_to_coarse, orth_on_coarse);
		types.push_back(fine_type);
		types.push_back(coarse_type);
	}
	void getRanks(std::deque<int> &ranks)
	{
		ranks.push_back(nbr_info->rank);
		ranks.push_back(nbr_info->rank);
	}
	void setLocalIndexes(const std::map<int, int> &rev_map)
	{
		this->local_index = rev_map.at(this->id);
		auto it           = rev_map.find(coarse_id);
		if (it != rev_map.end())
			coarse_local_index = it->second;
	}
	void setGlobalIndexes(const std::map<int, int> &rev_map)
	{
		this->global_index  = rev_map.at(this->local_index);
		coarse_global_index = rev_map.at(coarse_local_index);
	}
};
/**
 * @brief Represents the interfaces where the neighbors are at a finer refinement level.
 *
 * There will be 2^(D-1)+1 interfaces associated with this object. The interface that lines up with
 * this patch, and interface that lines up with the finer patches.
 *
 * @tparam D the number of Cartesian dimensions in a patch
 */
template <size_t D> class FineIfaceInfo : public IfaceInfo<D>
{
	public:
	/**
	 * @brief convenience pointer to associated NbrInfo object
	 */
	std::shared_ptr<FineNbrInfo<D>> nbr_info;
	/**
	 * @brief the ranks of the fine patches' interfaces
	 */
	std::array<int, Orthant<D - 1>::num_orthants> fine_ranks;
	/**
	 * @brief the ids of the fine patches' interfaces
	 */
	std::array<int, Orthant<D - 1>::num_orthants> fine_ids;
	/**
	 * @brief the local indexes of the fine patches' interfaces
	 */
	std::array<int, Orthant<D - 1>::num_orthants> fine_local_indexes;
	/**
	 * @brief the global indexes of the fine patches' interfaces
	 */
	std::array<int, Orthant<D - 1>::num_orthants> fine_global_indexes;
	/**
	 * @brief Construct a new FineIfaceInfo object
	 *
	 * @param pinfo the associated PatchInfo object
	 * @param s the side of the patch that the interface is on
	 */
	FineIfaceInfo(std::shared_ptr<PatchInfo<D>> pinfo, Side<D> s)
	{
		nbr_info   = pinfo->getFineNbrInfoPtr(s);
		this->id   = pinfo->id * Side<D>::num_sides + s.getIndex();
		this->rank = pinfo->rank;
		for (size_t i = 0; i < fine_ids.size(); i++) {
			fine_ids[i]   = nbr_info->ids[i] * Side<D>::num_sides + s.opposite().getIndex();
			fine_ranks[i] = nbr_info->ranks[i];
		}
	}
	void getIdxAndTypes(std::deque<int> &idx, std::deque<IfaceType<D>> &types)
	{
		idx.push_back(this->local_index);
		types.push_back(IfaceType<D>::coarse_to_coarse);
		for (size_t i = 0; i < fine_local_indexes.size(); i++) {
			idx.push_back(fine_local_indexes[i]);
			IfaceType<D> type(IfaceType<D>::coarse_to_fine, i);
			types.push_back(type);
		}
	}
	void getIds(std::deque<int> &ids)
	{
		ids.push_back(this->id);
		for (size_t i = 0; i < fine_ids.size(); i++) {
			ids.push_back(fine_ids[i]);
		}
	}
	void getLocalIndexes(std::deque<int> &idx)
	{
		idx.push_back(this->local_index);
		for (size_t i = 0; i < fine_ids.size(); i++) {
			idx.push_back(fine_local_indexes[i]);
		}
	}
	void getGlobalIndexes(std::deque<int> &idx)
	{
		idx.push_back(this->global_index);
		for (size_t i = 0; i < fine_ids.size(); i++) {
			idx.push_back(fine_global_indexes[i]);
		}
	}
	void getIfaceTypes(std::deque<IfaceType<D>> &types)
	{
		types.push_back(IfaceType<D>::coarse_to_coarse);
		for (size_t i = 0; i < fine_ids.size(); i++) {
			IfaceType<D> type(IfaceType<D>::coarse_to_fine, i);
			types.push_back(type);
		}
	}
	void getRanks(std::deque<int> &ranks)
	{
		ranks.push_back(-1);
		for (size_t i = 0; i < fine_ids.size(); i++) {
			ranks.push_back(nbr_info->ranks[i]);
		}
	}
	void setLocalIndexes(const std::map<int, int> &rev_map)
	{
		this->local_index = rev_map.at(this->id);
		for (size_t i = 0; i < fine_ids.size(); i++) {
			auto it = rev_map.find(this->fine_ids[i]);
			if (it != rev_map.end())
				fine_local_indexes[i] = it->second;
		}
	}
	void setGlobalIndexes(const std::map<int, int> &rev_map)
	{
		this->global_index = rev_map.at(this->local_index);
		for (size_t i = 0; i < fine_local_indexes.size(); i++) {
			fine_global_indexes[i] = rev_map.at(fine_local_indexes[i]);
		}
	}
};
/**
 * @brief
 *
 * @tparam D the number of Cartesian dimensions in a patch
 */
template <size_t D> struct SchurInfo : public Serializable {
	/**
	 * @brief Pointer to associated PatchInfo object
	 */
	std::shared_ptr<PatchInfo<D>> pinfo;
	/**
	 * @brief Array of IfaceInfo objects
	 */
	std::array<std::shared_ptr<IfaceInfo<D>>, Side<D>::num_sides> iface_info;
	/**
	 * @brief Construct a new empty SchurInfo object
	 *
	 */
	SchurInfo()
	{
		iface_info.fill(nullptr);
	}
	/**
	 * @brief Construct a new SchurInfo object.
	 *
	 * Fills in information from the given PatchInfo object.
	 */
	SchurInfo(std::shared_ptr<PatchInfo<D>> &pinfo)
	{
		this->pinfo = pinfo;
		iface_info.fill(nullptr);

		// create iface objects
		for (Side<D> s : Side<D>::getValues()) {
			if (pinfo->hasNbr(s)) {
				switch (pinfo->getNbrType(s)) {
					case NbrType::Normal:
						getIfaceInfoPtr(s).reset(new NormalIfaceInfo<D>(pinfo, s));
						break;
					case NbrType::Fine:
						getIfaceInfoPtr(s).reset(new FineIfaceInfo<D>(pinfo, s));
						break;
					case NbrType::Coarse:
						getIfaceInfoPtr(s).reset(new CoarseIfaceInfo<D>(pinfo, s));
						break;
				}
			}
		}
	}
	/**
	 * @brief Get a reference to the IfaceInfo pointer
	 *
	 * @param s the side of the patch that the interface is on
	 */
	std::shared_ptr<IfaceInfo<D>> &getIfaceInfoPtr(Side<D> s)
	{
		return iface_info[s.getIndex()];
	}
	/*
	 * @brief Get a reference to the NormalIfaceInfo object.
	 *
	 * If there is not a NormalIfaceInfo object on this side, the behavior is undefined.
	 *
	 * @param s the side of the patch that the interface is on
	 */
	NormalIfaceInfo<D> &getNormalIfaceInfo(Side<D> s)
	{
		return *std::dynamic_pointer_cast<NormalIfaceInfo<D>>(iface_info[s.getIndex()]);
	}
	/*
	 * @brief Get a reference to the CoarseIfaceInfo object.
	 *
	 * If there is not a CoarseIfaceInfo object on this side, the behavior is undefined.
	 *
	 * @param s the side of the patch that the interface is on
	 */
	CoarseIfaceInfo<D> &getCoarseIfaceInfo(Side<D> s)
	{
		return *std::dynamic_pointer_cast<CoarseIfaceInfo<D>>(iface_info[s.getIndex()]);
	}
	/*
	 * @brief Get a reference to the FineIfaceInfo object.
	 *
	 * If there is not a FineIfaceInfo object on this side, the behavior is undefined.
	 *
	 * @param s the side of the patch that the interface is on
	 */
	FineIfaceInfo<D> &getFineIfaceInfo(Side<D> s)
	{
		return *std::dynamic_pointer_cast<FineIfaceInfo<D>>(iface_info[s.getIndex()]);
	}
	/*
	 * @brief Get a reference to the NormalIfaceInfo object.
	 *
	 * If there is not a NormalIfaceInfo object on this side, the behavior is undefined.
	 *
	 * @param s the side of the patch that the interface is on
	 */
	const NormalIfaceInfo<D> &getNormalIfaceInfo(Side<D> s) const
	{
		return *std::dynamic_pointer_cast<NormalIfaceInfo<D>>(iface_info[s.getIndex()]);
	}
	/*
	 * @brief Get a reference to the CoarseIfaceInfo object.
	 *
	 * If there is not a CoarseIfaceInfo object on this side, the behavior is undefined.
	 *
	 * @param s the side of the patch that the interface is on
	 */
	const CoarseIfaceInfo<D> &getCoarseIfaceInfo(Side<D> s) const
	{
		return *std::dynamic_pointer_cast<CoarseIfaceInfo<D>>(iface_info[s.getIndex()]);
	}
	/*
	 * @brief Get a reference to the FineIfaceInfo object.
	 *
	 * If there is not a FineIfaceInfo object on this side, the behavior is undefined.
	 *
	 * @param s the side of the patch that the interface is on
	 */
	const FineIfaceInfo<D> &getFineIfaceInfo(Side<D> s) const
	{
		return *std::dynamic_pointer_cast<FineIfaceInfo<D>>(iface_info[s.getIndex()]);
	}
	/**
	 * @brief
	 *
	 * @param ifaces
	 * @param off_proc_ifaces
	 * @param incoming_procs
	 */

	std::deque<int> getIds()
	{
		std::deque<int> retval;
		for (Side<D> s : Side<D>::getValues()) {
			if (pinfo->hasNbr(s)) {
				getIfaceInfoPtr(s)->getIds(retval);
			}
		}
		return retval;
	}
	void setLocalIndexes(const std::map<int, int> &rev_map)
	{
		for (Side<D> s : Side<D>::getValues()) {
			if (pinfo->hasNbr(s)) {
				getIfaceInfoPtr(s)->setLocalIndexes(rev_map);
			}
		}
	}
	void setGlobalIndexes(const std::map<int, int> &rev_map)
	{
		for (Side<D> s : Side<D>::getValues()) {
			if (pinfo->hasNbr(s)) {
				getIfaceInfoPtr(s)->setGlobalIndexes(rev_map);
			}
		}
	}
	int getIfaceLocalIndex(Side<D> s) const
	{
		return iface_info[s.getIndex()]->local_index;
	}
	int serialize(char *buffer) const
	{
		BufferWriter writer(buffer);
		writer << *pinfo;
		return writer.getPos();
	}
	int deserialize(char *buffer)
	{
		BufferReader reader(buffer);
		pinfo.reset(new PatchInfo<D>());
		reader >> *pinfo;
		// create iface objects
		for (Side<D> s : Side<D>::getValues()) {
			if (pinfo->hasNbr(s)) {
				switch (pinfo->getNbrType(s)) {
					case NbrType::Normal:
						getIfaceInfoPtr(s).reset(new NormalIfaceInfo<D>(pinfo, s));
						break;
					case NbrType::Fine:
						getIfaceInfoPtr(s).reset(new FineIfaceInfo<D>(pinfo, s));
						break;
					case NbrType::Coarse:
						getIfaceInfoPtr(s).reset(new CoarseIfaceInfo<D>(pinfo, s));
						break;
				}
			}
		}
		return reader.getPos();
	}
};
extern template struct SchurInfo<2>;
extern template struct SchurInfo<3>;
} // namespace Schur
} // namespace Thunderegg
#endif
