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
#include <Thunderegg/SchurInfo.h>
#include <bitset>
#include <map>
#include <mpi.h>
#include <set>
#include <vector>
namespace Thunderegg
{
/**
 * @brief This is a set of SchurInfo objects for an interface.
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
	 * @brief the side of the patchinfo object that this interface is on.
	 */
	std::vector<Side<D>> sides;
	/**
	 * @brief the type of interface that this interface is on the patchinfo object
	 */
	std::vector<IfaceType<D>> types;
	/**
	 * @brief the set of PatchInfo objects associated with this interface
	 */
	std::vector<std::shared_ptr<SchurInfo<D>>> sinfos;
	/**
	 * @brief Create IfaceSet objects from a given collection of SchurInfo objects.
	 *
	 * This method will communicate any necessary information over mpi
	 *
	 * @tparam Iter the iterator type
	 * @param begin the begining of the collection of std::shared_ptr<SchurInfo<D>> objects
	 * @param end the end of the collection of std::shared_ptr<SchurInfo<D>> objects
	 * @return std::map<int, IfaceSet<D>> a map of interface id to IfaceSet oject
	 */
	template <class Iter> static std::map<int, IfaceSet<D>> EnumerateIfaces(Iter begin, Iter end);

	/**
	 * @brief Get a set of neighboring interfaces. (The interfaces that are on patches connected to
	 * this interface)
	 *
	 * @return std::set<int>
	 */
	std::set<int> getNbrs() const
	{
		std::set<int> retval;
		for (auto sinfo : sinfos) {
			for (auto info : sinfo->iface_info) {
				if (info != nullptr && info->id != -1 && info->id != id) {
					retval.insert(info->id);
				}
			}
		}
		return retval;
	}
	/**
	 * @brief Add an Iface to this set.
	 */
	void insert(int id, Side<D> s, std::shared_ptr<SchurInfo<D>> sinfo)
	{
		this->id = id;
		sides.push_back(s);
		sinfos.push_back(sinfo);
		switch (sinfo->pinfo->getNbrType(s)) {
			case NbrType::Normal:
				types.push_back(IfaceType<D>::normal);
				break;
			case NbrType::Fine: {
				const FineIfaceInfo<D> &info = sinfo->getFineIfaceInfo(s);
				if (info.id == id) {
					types.push_back(IfaceType<D>::coarse_to_coarse);
				} else {
					for (int i = 0; i < Orthant<D - 1>::num_orthants; i++) {
						if (info.fine_ids[i] == id) {
							types.push_back(IfaceType<D>::coarse_to_fine);
							types.back().setOrthant(i);
							break;
						}
					}
				}
			} break;
			case NbrType::Coarse: {
				const CoarseIfaceInfo<D> &info = sinfo->getCoarseIfaceInfo(s);
				if (info.id == id) {
					types.push_back(IfaceType<D>::fine_to_fine);
				} else {
					types.push_back(IfaceType<D>::fine_to_coarse);
					types.back().setOrthant(info.orth_on_coarse);
				}
			} break;
		}
	}
	/**
	 * @brief Add an Ifaces in an IfaceSet to this set.
	 */
	void merge(IfaceSet<D> ifs)
	{
		id = ifs.id;
		for (size_t i = 0; i < ifs.sinfos.size(); i++) {
			sinfos.push_back(ifs.sinfos[i]);
			sides.push_back(ifs.sides[i]);
			types.push_back(ifs.types[i]);
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
		for (auto &sinfo : sinfos) {
			sinfo->setLocalIndexes(rev_map);
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
		for (auto &sinfo : sinfos) {
			sinfo->setGlobalIndexes(rev_map);
		}
	}
	int serialize(char *buffer) const
	{
		BufferWriter writer(buffer);
		writer << id;
		writer << id_global;
		int size = sinfos.size();
		writer << size;
		for (size_t i = 0; i < sinfos.size(); i++) {
			writer << sides[i];
			writer << types[i];
			writer << *sinfos[i];
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
			Side<D> s;
			reader >> s;
			sides.push_back(s);
			IfaceType<D> type;
			reader >> type;
			types.push_back(type);
			std::shared_ptr<SchurInfo<D>> sinfo(new SchurInfo<D>());
			reader >> *sinfo;
			sinfos.push_back(sinfo);
		}
		return reader.getPos();
	}
};
template <size_t D>
template <class Iter>
std::map<int, IfaceSet<D>> IfaceSet<D>::EnumerateIfaces(Iter begin, Iter end)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	std::map<int, IfaceSet<D>>                ifaces;
	std::set<int>                             incoming_procs;
	std::map<int, std::map<int, IfaceSet<D>>> off_proc_ifaces;
	for (; begin != end; begin++) {
		std::shared_ptr<SchurInfo<D>> sinfo = *begin;
		for (Side<D> s : Side<D>::getValues()) {
			if (sinfo->pinfo->hasNbr(s)) {
				switch (sinfo->pinfo->getNbrType(s)) {
					case NbrType::Normal: {
						const NormalIfaceInfo<D> &info = sinfo->getNormalIfaceInfo(s);
						if (info.rank == rank) {
							ifaces[info.id].insert(info.id, s, sinfo);
							if (info.nbr_info->rank != rank) {
								incoming_procs.insert(info.nbr_info->rank);
							}
						} else {
							off_proc_ifaces[info.rank][info.id].insert(info.id, s, sinfo);
						}
					} break;
					case NbrType::Fine: {
						const FineIfaceInfo<D> &info = sinfo->getFineIfaceInfo(s);
						ifaces[info.id].insert(info.id, s, sinfo);
						for (int i = 0; i < Orthant<D - 1>::num_orthants; i++) {
							if (info.fine_ranks[i] == rank) {
								ifaces[info.fine_ids[i]].insert(info.fine_ids[i], s, sinfo);
							} else {
								off_proc_ifaces[info.fine_ranks[i]][info.fine_ids[i]].insert(
								info.fine_ids[i], s, sinfo);
								incoming_procs.insert(info.fine_ranks[i]);
							}
						}
					} break;
					case NbrType::Coarse: {
						const CoarseIfaceInfo<D> &info = sinfo->getCoarseIfaceInfo(s);
						ifaces[info.id].insert(info.id, s, sinfo);
						if (info.coarse_rank == rank) {
							ifaces[info.coarse_id].insert(info.coarse_id, s, sinfo);
						} else {
							off_proc_ifaces[info.coarse_rank][info.coarse_id].insert(info.coarse_id,
							                                                         s, sinfo);
							incoming_procs.insert(info.coarse_rank);
						}
					} break;
				}
			}
		}
	}
	// send info
	std::deque<char *>       buffers;
	std::vector<MPI_Request> send_requests;
	for (auto &p : off_proc_ifaces) {
		int dest = p.first;
		int size = 0;
		for (auto q : p.second) {
			IfaceSet<D> &iface = q.second;
			size += iface.serialize(nullptr);
		}
		char *buffer = new char[size];
		buffers.push_back(buffer);
		int pos = 0;
		for (auto q : p.second) {
			IfaceSet<D> &iface = q.second;
			pos += iface.serialize(buffer + pos);
		}
		MPI_Request request;
		MPI_Isend(buffer, size, MPI_BYTE, dest, 0, MPI_COMM_WORLD, &request);
		send_requests.push_back(request);
	}
	// recv info
	for (int src : incoming_procs) {
		MPI_Status status;
		MPI_Probe(src, 0, MPI_COMM_WORLD, &status);
		int size;
		MPI_Get_count(&status, MPI_BYTE, &size);
		char *buffer = new char[size];

		MPI_Recv(buffer, size, MPI_BYTE, src, 0, MPI_COMM_WORLD, &status);

		BufferReader reader(buffer);
		while (reader.getPos() < size) {
			IfaceSet<D> ifs;
			reader >> ifs;
			ifaces[ifs.id].merge(ifs);
		}

		delete[] buffer;
	}
	// wait for all
	MPI_Waitall(send_requests.size(), &send_requests[0], MPI_STATUSES_IGNORE);
	// delete send buffers
	for (char *buffer : buffers) {
		delete[] buffer;
	}
	return ifaces;
}
} // namespace Thunderegg
#endif
