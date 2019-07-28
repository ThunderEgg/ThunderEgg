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

#ifndef IFACE_H
#define IFACE_H
#include <Thunderegg/BufferIO.h>
#include <Thunderegg/IfaceType.h>
#include <Thunderegg/Side.h>
#include <bitset>
#include <map>
#include <set>
#include <vector>
namespace Thunderegg
{
/**
 * @brief  Contains information for an a specific interface of a patch. (Does not include
 * neighboring information.)
 *
 * @tparam D The number of cartesian dimensions in a patch.
 */
template <size_t D> struct Iface {
	/**
	 * @brief The type of interface.
	 */
	IfaceType<D> type;
	/**
	 * @brief The side of the patch that the interface is on.
	 */
	Side<D> s;
	/**
	 * @brief The id of each interface on each side of the patch.
	 */
	std::array<int, Side<D>::num_sides> ids;
	/**
	 * @brief The local index of the interfaces on each side of the patch.
	 */
	std::array<int, Side<D>::num_sides> local_id;
	/**
	 * @brief The global index of the interfaces on each side of the patch.
	 */
	std::array<int, Side<D>::num_sides> global_id;
	/**
	 * @brief The boundary conditions on the patch.
	 */
	std::bitset<Side<D>::num_sides> neumann;
	/**
	 * @brief Construct a new Iface object
	 *
	 * Everything is set to -1
	 */
	Iface()
	{
		ids.fill(-1);
		local_id.fill(-1);
		global_id.fill(-1);
	}
	/**
	 * @brief Construct a new Iface object
	 *
	 * @param ids the ids of the interfaces on each side of the patch
	 * @param type the type of iface
	 * @param s the side of the patch that the iface is on
	 * @param neumann the boundary conditions of the patch
	 */
	Iface(std::array<int, Side<D>::num_sides> ids, IfaceType<D> type, Side<D> s,
	      std::bitset<Side<D>::num_sides> neumann)
	{
		this->ids = ids;
		local_id.fill(-1);
		global_id.fill(-1);
		this->type    = type;
		this->s       = s;
		this->neumann = neumann;
	}
};
/**
 * @brief This is a set of Iface objects for an interface.
 *
 * This will contain information from each patch that the interface is on.
 *
 * @tparam D the number of cartesian dimensions on a patch.
 */
template <size_t D> struct IfaceSet : public Serializable {
	/**
	 * @brief The id of the interface
	 */
	int id = -1;
	/**
	 * @brief The local index of the interface
	 */
	int id_local = -1;
	/**
	 * @brief The global index of the interface
	 */
	int id_global = -1;
	/**
	 * @brief the set of Iface objects
	 */
	std::vector<Iface<D>> ifaces;
	/**
	 * @brief Get a set of neighboring interfaces. (The interfaces that are on patches connected to
	 * this interface)
	 *
	 * @return std::set<int>
	 */
	std::set<int> getNbrs() const
	{
		std::set<int> retval;
		for (const Iface<D> &iface : ifaces) {
			for (const int i : iface.ids) {
				if (i != -1 && i != id) { retval.insert(i); }
			}
		}
		return retval;
	}
	/**
	 * @brief Add an Iface to this set.
	 */
	void insert(Iface<D> i)
	{
		ifaces.push_back(i);
	}
	/**
	 * @brief Add an Ifaces in an IfaceSet to this set.
	 */
	void insert(IfaceSet<D> ifs)
	{
		id = ifs.id;
		for (Iface<D> &i : ifs.ifaces) {
			ifaces.push_back(i);
		}
	}
	/**
	 * @brief Set the local indexes in the Iface objects
	 *
	 * @param rev_map map from id to local_index
	 */
	void setLocalIndexes(const std::map<int, int> &rev_map)
	{
		id_local = rev_map.at(id);
		for (Iface<D> &iface : ifaces) {
			for (int i = 0; i < Side<D>::num_sides; i++) {
				if (iface.ids[i] != -1) { iface.local_id[i] = rev_map.at(iface.ids[i]); }
			}
		}
	}
	/**
	 * @brief Set the global indexes in the NbrInfo objects
	 *
	 * @param rev_map map form local_index to global_index
	 */
	void setGlobalIndexes(const std::map<int, int> &rev_map)
	{
		id_global = rev_map.at(id_local);
		for (Iface<D> &iface : ifaces) {
			for (int i = 0; i < Side<D>::num_sides; i++) {
				if (iface.local_id[i] != -1) { iface.global_id[i] = rev_map.at(iface.local_id[i]); }
			}
		}
	}
	int serialize(char *buffer) const
	{
		BufferWriter writer(buffer);
		writer << id;
		writer << id_global;
		int size = ifaces.size();
		writer << size;
		for (const Iface<D> &i : ifaces) {
			writer << i;
		}
		return writer.getPos();
	}
	int deserialize(char *buffer)
	{
		BufferReader reader(buffer);
		reader >> id;
		reader >> id_global;
		int size = 0;
		reader >> size;
		for (int i = 0; i < size; i++) {
			Iface<D> iface;
			reader >> iface;
			ifaces.push_back(iface);
		}
		return reader.getPos();
	}
};
} // namespace Thunderegg
#endif
