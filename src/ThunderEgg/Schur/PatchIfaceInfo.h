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

#ifndef THUNDEREGG_SCHUR_PATCHIFACEINFO_H
#define THUNDEREGG_SCHUR_PATCHIFACEINFO_H
#include <ThunderEgg/RuntimeError.h>
#include <ThunderEgg/Schur/CoarseIfaceInfo.h>
#include <ThunderEgg/Schur/FineIfaceInfo.h>
#include <ThunderEgg/Schur/NormalIfaceInfo.h>
namespace ThunderEgg
{
namespace Schur
{
/**
 * @brief This decorates a PatchInfo object with a IfaceInfo object for each side of the patch.
 *
 * @tparam D the number of Cartesian dimensions in a patch
 */
template <size_t D> class PatchIfaceInfo : public Serializable
{
	public:
	/**
	 * @brief Pointer to associated PatchInfo object
	 */
	std::shared_ptr<const PatchInfo<D>> pinfo;
	/**
	 * @brief Array of IfaceInfo objects
	 */
	std::array<std::shared_ptr<IfaceInfo<D>>, Side<D>::num_sides> iface_info;
	/**
	 * @brief Construct a new PatchIfaceInfo object
	 */
	PatchIfaceInfo() = default;
	/**
	 * @brief Construct a new PatchIfaceInfo object.
	 *
	 * Fills in information from the given PatchInfo object.
	 */
	explicit PatchIfaceInfo(std::shared_ptr<const PatchInfo<D>> pinfo) : pinfo(pinfo)
	{
		// create iface objects
		for (Side<D> s : Side<D>::getValues()) {
			if (pinfo->hasNbr(s)) {
				switch (pinfo->getNbrType(s)) {
					case NbrType::Normal:
						setIfaceInfo(s, std::make_shared<NormalIfaceInfo<D>>(pinfo, s));
						break;
					case NbrType::Fine:
						setIfaceInfo(s, std::make_shared<FineIfaceInfo<D>>(pinfo, s));
						break;
					case NbrType::Coarse:
						setIfaceInfo(s, std::make_shared<CoarseIfaceInfo<D>>(pinfo, s));
						break;
					default:
						throw RuntimeError("Unsupported NbrType");
				}
			}
		}
	}
	/**
	 * @brief Set the IfaceInfo object on a given side of the patch
	 *
	 * @param s the side of the patch
	 * @param info the IfaceInfo object
	 */
	void setIfaceInfo(Side<D> s, std::shared_ptr<IfaceInfo<D>> info)
	{
		iface_info[s.getIndex()] = info;
	}
	/**
	 * @brief Get the IfaceInfo object on a given side of the patch
	 *
	 * @param s the side of the patch
	 * @return std::shared_ptr<IfaceInfo<D>> the IfaceInfo object
	 */
	std::shared_ptr<IfaceInfo<D>> getIfaceInfo(Side<D> s)
	{
		return iface_info[s.getIndex()];
	}
	/**
	 * @brief Get the IfaceInfo object on a given side of the patch
	 *
	 * @param s the side of the patch
	 * @return std::shared_ptr<const IfaceInfo<D>> the IfaceInfo object
	 */
	std::shared_ptr<const IfaceInfo<D>> getIfaceInfo(Side<D> s) const
	{
		return iface_info[s.getIndex()];
	}
	/**
	 * @brief Get the NormalIfaceInfo object on a given side of the patch
	 *
	 * @param s the side of the object
	 * @return std::shared_ptr<const NormalIfaceInfo<D>> the NormalIfaceInfo object, nullptr if
	 * there is not a NormalIfaceInfo object on the given side
	 */
	std::shared_ptr<const NormalIfaceInfo<D>> getNormalIfaceInfo(Side<D> s) const
	{
		return std::dynamic_pointer_cast<const NormalIfaceInfo<D>>(iface_info[s.getIndex()]);
	}
	/**
	 * @brief Get the CoarseIfaceInfo object on a given side of the patch
	 *
	 * @param s the side of the object
	 * @return std::shared_ptr<const CoarseIfaceInfo<D>> the CoarseIfaceInfo object, nullptr if
	 * there is not a CoarseIfaceInfo object on the given side
	 */
	std::shared_ptr<const CoarseIfaceInfo<D>> getCoarseIfaceInfo(Side<D> s) const
	{
		return std::dynamic_pointer_cast<const CoarseIfaceInfo<D>>(iface_info[s.getIndex()]);
	}
	/**
	 * @brief Get the FineIfaceInfo object on a given side of the patch
	 *
	 * @param s the side of the object
	 * @return std::shared_ptr<const FineIfaceInfo<D>> the FineIfaceInfo object, nullptr if
	 * there is not a FineIfaceInfo object on the given side
	 */
	std::shared_ptr<const FineIfaceInfo<D>> getFineIfaceInfo(Side<D> s)
	{
		return std::dynamic_pointer_cast<const FineIfaceInfo<D>>(iface_info[s.getIndex()]);
	}
	/**
	 * @brief Set the local indexes from interface ids
	 *
	 * @param rev_map the map from id to local index
	 */
	void setLocalIndexesFromId(const std::map<int, int> &rev_map)
	{
		for (Side<D> s : Side<D>::getValues()) {
			if (getIfaceInfo(s) != nullptr) {
				getIfaceInfo(s)->setLocalIndexesFromId(rev_map);
			}
		}
	}
	/**
	 * @brief Set the global indexes from the local indexes of the interafce
	 *
	 * @param rev_map map from local index to global index
	 */
	void setGlobalIndexesFromLocalIndex(const std::map<int, int> &rev_map)
	{
		for (Side<D> s : Side<D>::getValues()) {
			if (getIfaceInfo(s) != nullptr) {
				getIfaceInfo(s)->setGlobalIndexesFromLocalIndex(rev_map);
			}
		}
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
		auto         pinfo_in = std::make_shared<PatchInfo<D>>();
		reader >> *pinfo_in;
		pinfo = pinfo_in;
		// create iface objects
		for (Side<D> s : Side<D>::getValues()) {
			if (pinfo->hasNbr(s)) {
				switch (pinfo->getNbrType(s)) {
					case NbrType::Normal:
						setIfaceInfo(s, std::make_shared<NormalIfaceInfo<D>>(pinfo, s));
						break;
					case NbrType::Fine:
						setIfaceInfo(s, std::make_shared<FineIfaceInfo<D>>(pinfo, s));
						break;
					case NbrType::Coarse:
						setIfaceInfo(s, std::make_shared<CoarseIfaceInfo<D>>(pinfo, s));
						break;
					default:
						throw RuntimeError("Unsupported NbrType");
				}
			}
		}
		return reader.getPos();
	}
};
extern template class PatchIfaceInfo<2>;
extern template class PatchIfaceInfo<3>;
} // namespace Schur
} // namespace ThunderEgg
#endif
