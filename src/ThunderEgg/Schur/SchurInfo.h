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

#ifndef THUNDEREGG_SCHUR_SCHURINFO_H
#define THUNDEREGG_SCHUR_SCHURINFO_H
#include <ThunderEgg/Schur/CoarseIfaceInfo.h>
#include <ThunderEgg/Schur/FineIfaceInfo.h>
#include <ThunderEgg/Schur/NormalIfaceInfo.h>
namespace ThunderEgg
{
namespace Schur
{
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
} // namespace ThunderEgg
#endif
