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

#ifndef THUNDEREGG_SCHUR_INTERFACE_H
#define THUNDEREGG_SCHUR_INTERFACE_H
#include <ThunderEgg/BufferReader.h>
#include <ThunderEgg/BufferWriter.h>
#include <ThunderEgg/Schur/IfaceType.h>
#include <ThunderEgg/Schur/PatchIfaceInfo.h>
#include <bitset>
#include <map>
#include <mpi.h>
#include <set>
#include <vector>
namespace ThunderEgg
{
namespace Schur
{
/**
 * @brief This is a set of PatchIfaceInfo objects for an interface.
 *
 * This will contain information from each patch that the interface is on.
 *
 * @tparam D the number of cartesian dimensions on a patch.
 */
template <size_t D> class Interface : public Serializable
{
	public:
	/**
	 * @brief The id of the interface
	 */
	int id = -1;
	/**
	 * @brief The local index of the interface
	 */
	int local_index = -1;
	/**
	 * @brief The global index of the interface
	 */
	int global_index = -1;
	/**
	 * @brief A struct for each patch associated with this interface
	 *
	 * The struct contains the side of the patch that this interface lies on, the IfaceType of this
	 * interface on that patch, and the patch.
	 */
	struct SideTypePiinfo {
		/**
		 * @brief the side of the PatchIfaceInfo object that this interface is on.
		 */
		Side<D> side;
		/**
		 * @brief the IfaceType of this interface on the cooresponding PatchIfaceInfo object.
		 */
		IfaceType<D> type;
		/**
		 * @brief the PatchIfaceInfo object
		 */
		std::shared_ptr<const PatchIfaceInfo<D>> piinfo;
		/**
		 * @brief Construct a new Side Type Piinfo object
		 *
		 * @param side the side of the PatchIfaceInfo object that this interface is on.
		 * @param type the IfaceType of this interface on the cooresponding PatchIfaceInfo object.
		 * @param piinfo the PatchIfaceInfo object
		 */
		SideTypePiinfo(Side<D> side, IfaceType<D> type,
		               std::shared_ptr<const PatchIfaceInfo<D>> piinfo)
		: side(side), type(type), piinfo(piinfo)
		{
		}
	};
	/**
	 * @brief SideTypePiinfo structs associated with this interface.
	 *
	 * The struct contains the side of the patch that this interface lies on, the IfaceType of this
	 * interface on that patch, and the patch.
	 */
	std::vector<SideTypePiinfo> patches;

	/**
	 * @brief Construct a new Interface object
	 */
	Interface() = default;

	/**
	 * @brief Construct a new Interface object
	 *
	 * @param id the id of the interface
	 */
	explicit Interface(int id) : id(id) {}

	/**
	 * @brief Add an associated patch to this interface
	 *
	 * @param s the side of the patch that this interface is on
	 * @param piinfo the PatchIfaceInfo object for the patch
	 */
	void insert(Side<D> s, std::shared_ptr<const PatchIfaceInfo<D>> piinfo)
	{
		NbrType nbr_type = piinfo->pinfo->getNbrType(s);
		if (nbr_type == NbrType::Normal) {
			patches.emplace_back(s, IfaceType<D>::Normal(), piinfo);
		} else if (nbr_type == NbrType::Fine) {
			auto info = piinfo->getFineIfaceInfo(s);
			if (info->id == id) {
				patches.emplace_back(s, IfaceType<D>::CoarseToCoarse(), piinfo);
			} else {
				for (Orthant<D - 1> o : Orthant<D - 1>::getValues()) {
					if (info->fine_ids[o.getIndex()] == id) {
						patches.emplace_back(s, IfaceType<D>::CoarseToFine(o), piinfo);
						break;
					}
				}
			}
		} else if (nbr_type == NbrType::Coarse) {
			auto info = piinfo->getCoarseIfaceInfo(s);
			if (info->id == id) {
				patches.emplace_back(s, IfaceType<D>::FineToFine(info->orth_on_coarse), piinfo);
			} else {
				patches.emplace_back(s, IfaceType<D>::FineToCoarse(info->orth_on_coarse), piinfo);
			}
		} else {
			throw RuntimeError("Unsupported NbrType");
		}
	}
	/**
	 * @brief Add the patches from the given Interface object to this object
	 *
	 * @param ifs the Interface object from which to add the patches
	 */
	void merge(Interface<D> ifs)
	{
		for (auto patch : ifs.patches) {
			patches.push_back(patch);
		}
	}
	int serialize(char *buffer) const
	{
		BufferWriter writer(buffer);
		writer << id;
		writer << global_index;
		int size = (int) patches.size();
		writer << size;
		for (auto patch : patches) {
			writer << patch.side;
			writer << patch.type;
			writer << *patch.piinfo;
		}
		return writer.getPos();
	}
	int deserialize(char *buffer)
	{
		BufferReader reader(buffer);
		reader >> id;
		reader >> global_index;
		int size = 0;
		reader >> size;
		for (int i = 0; i < size; i++) {
			Side<D>      s;
			IfaceType<D> type;
			auto         piinfo = std::make_shared<PatchIfaceInfo<D>>();
			reader >> s;
			reader >> type;
			reader >> *piinfo;
			patches.emplace_back(s, type, piinfo);
		}
		return reader.getPos();
	}

	private:
	/**
	 * @brief Insert the given PatchIfaceInfo object into the Interface map
	 *
	 * @param rank_id_iface_map the map from rank to the interface's id to Interface
	 * @param rank the rank of the interface
	 * @param id the id of the interface
	 * @param s the side of the patch that the interface is on
	 * @param piinfo the PatchIfaceInfo object
	 */
	static void InsertPatchToInterface(
	std::map<int, std::map<int, std::shared_ptr<Interface<D>>>> &rank_id_iface_map, int rank,
	int id, Side<D> s, std::shared_ptr<const PatchIfaceInfo<D>> piinfo)
	{
		std::shared_ptr<Interface<D>> &iface_ptr = rank_id_iface_map[rank][id];
		if (iface_ptr == nullptr) {
			iface_ptr.reset(new Interface<D>(id));
		}
		iface_ptr->insert(s, piinfo);
	}
	/**
	 * @brief Will insert the interface shared with a normal neighbor into rank_id_iface_map, and
	 * will add any ranks that are sending messages to us into incoming_procs
	 *
	 * @param rank_id_iface_map the map from rank to the interface's id to Interface
	 * @param incoming_procs the set of ranks that are sending messages this rank
	 * @param piinfo the PatchIfaceInfo object
	 * @param s the side of the patch that the interface is on
	 */
	static void InsertInterfaceWithNormalNbr(
	std::map<int, std::map<int, std::shared_ptr<Interface<D>>>> &rank_id_iface_map,
	std::set<int> &incoming_procs, std::shared_ptr<const PatchIfaceInfo<D>> piinfo, Side<D> s)
	{
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		auto info = piinfo->getNormalIfaceInfo(s);
		InsertPatchToInterface(rank_id_iface_map, info->rank, info->id, s, piinfo);
		if (info->rank == piinfo->pinfo->rank && info->nbr_info->rank != piinfo->pinfo->rank) {
			incoming_procs.insert(info->nbr_info->rank);
		}
	}
	/**
	 * @brief Will insert the interfaces shared with a finer neighbor into rank_id_iface_map, and
	 * will add any ranks that are sending messages to us into incoming_procs
	 *
	 * @param rank_id_iface_map the map from rank to the interface's id to Interface
	 * @param incoming_procs the set of ranks that are sending messages this rank
	 * @param piinfo the PatchIfaceInfo object
	 * @param s the side of the patch that the interface is on
	 */
	static void InsertInterfacesWithFineNbr(
	std::map<int, std::map<int, std::shared_ptr<Interface<D>>>> &rank_id_iface_map,
	std::set<int> &incoming_procs, std::shared_ptr<const PatchIfaceInfo<D>> piinfo, Side<D> s)
	{
		auto info = piinfo->getFineIfaceInfo(s);

		InsertPatchToInterface(rank_id_iface_map, info->rank, info->id, s, piinfo);

		for (size_t i = 0; i < Orthant<D - 1>::num_orthants; i++) {
			InsertPatchToInterface(rank_id_iface_map, info->fine_ranks[i], info->fine_ids[i], s,
			                       piinfo);

			if (info->fine_ranks[i] != piinfo->pinfo->rank) {
				incoming_procs.insert(info->fine_ranks[i]);
			}
		}
	}
	/**
	 * @brief Will insert the interfaces shared with a coarser neighbor into rank_id_iface_map, and
	 * will add any ranks that are sending messages to us into incoming_procs
	 *
	 * @param rank_id_iface_map the map from rank to the interface's id to Interface
	 * @param incoming_procs the set of ranks that are sending messages this rank
	 * @param piinfo the PatchIfaceInfo object
	 * @param s the side of the patch that the interface is on
	 */
	static void InsertInterfacesWithCoarseNbr(
	std::map<int, std::map<int, std::shared_ptr<Interface<D>>>> &rank_id_iface_map,
	std::set<int> &incoming_procs, std::shared_ptr<const PatchIfaceInfo<D>> piinfo, Side<D> s)
	{
		auto info = piinfo->getCoarseIfaceInfo(s);

		InsertPatchToInterface(rank_id_iface_map, info->rank, info->id, s, piinfo);
		InsertPatchToInterface(rank_id_iface_map, info->coarse_rank, info->coarse_id, s, piinfo);

		if (info->coarse_rank != piinfo->pinfo->rank) {
			incoming_procs.insert(info->coarse_rank);
		}
	}

	public:
	/**
	 * @brief Will enumerate a map from interface id to this rank's interfaces, will also do any
	 * neccesary communication to get additional information. This is collective on all processors.
	 *
	 * @param piinfos vector of this ranks piinfo objects
	 * @return std::map<int, std::shared_ptr<Interface<D>>> the map from interface id to interface
	 */
	static void EnumerateIfacesFromPiinfoVector(
	std::vector<std::shared_ptr<const PatchIfaceInfo<D>>> piinfos,
	std::map<int, std::shared_ptr<Interface<D>>> &        id_to_iface_map,
	std::vector<std::shared_ptr<PatchIfaceInfo<D>>> &     off_proc_piinfos)
	{
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		std::map<int, std::map<int, std::shared_ptr<Interface<D>>>> rank_id_iface_map;
		std::set<int>                                               incoming_procs;
		for (auto piinfo : piinfos) {
			for (Side<D> s : Side<D>::getValues()) {
				if (piinfo->pinfo->hasNbr(s)) {
					switch (piinfo->pinfo->getNbrType(s)) {
						case NbrType::Normal:
							InsertInterfaceWithNormalNbr(rank_id_iface_map, incoming_procs, piinfo,
							                             s);
							break;
						case NbrType::Fine:
							InsertInterfacesWithFineNbr(rank_id_iface_map, incoming_procs, piinfo,
							                            s);
							break;
						case NbrType::Coarse:
							InsertInterfacesWithCoarseNbr(rank_id_iface_map, incoming_procs, piinfo,
							                              s);
							break;
						default:
							throw RuntimeError("Unsupported NbrType value");
					}
				}
			}
		}
		// send info
		std::deque<std::vector<char>> buffers;
		std::vector<MPI_Request>      send_requests;
		for (auto &p : rank_id_iface_map) {
			int dest = p.first;
			if (dest != rank) {
				int size = 0;
				for (auto q : p.second) {
					auto &iface = q.second;
					size += iface->serialize(nullptr);
				}
				buffers.emplace_back(size);
				BufferWriter writer(buffers.back().data());
				for (auto q : p.second) {
					auto &iface = q.second;
					writer << *iface;
				}
				MPI_Request request;
				MPI_Isend(buffers.back().data(), size, MPI_BYTE, dest, 0, MPI_COMM_WORLD, &request);
				send_requests.push_back(request);
			}
		}
		// recv info
		off_proc_piinfos.clear();
		size_t num_incoming = incoming_procs.size();
		for (size_t i = 0; i < num_incoming; i++) {
			MPI_Status status;
			MPI_Probe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
			int size;
			MPI_Get_count(&status, MPI_BYTE, &size);
			std::vector<char> buffer(size);

			MPI_Recv(buffer.data(), size, MPI_BYTE, status.MPI_SOURCE, 0, MPI_COMM_WORLD, &status);

			BufferReader reader(buffer.data());
			while (reader.getPos() < size) {
				Interface<D> ifs;
				reader >> ifs;
				rank_id_iface_map.at(rank).at(ifs.id)->merge(ifs);
				for (auto patch : ifs.patches) {
					// need to cast to remove const modifier
					off_proc_piinfos.push_back(
					std::const_pointer_cast<PatchIfaceInfo<D>>(patch.piinfo));
				}
			}
		}
		// wait for all
		MPI_Waitall((int) send_requests.size(), &send_requests[0], MPI_STATUSES_IGNORE);
		id_to_iface_map = rank_id_iface_map[rank];
	}
};
extern template class Interface<2>;
extern template class Interface<3>;
} // namespace Schur
} // namespace ThunderEgg
#endif
