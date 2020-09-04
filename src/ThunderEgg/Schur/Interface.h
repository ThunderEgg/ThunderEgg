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
	private:
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
	Interface(int id) : id(id) {}

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
		int size = patches.size();
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
	/**
	 * @brief Create Interface objects from a given collection of PatchIfaceInfo objects.
	 *
	 * This method will communicate any necessary information over mpi
	 *
	 * @tparam Iter the iterator type
	 * @param begin the begining of the collection of std::shared_ptr<PatchIfaceInfo<D>> objects
	 * @param end the end of the collection of std::shared_ptr<PatchIfaceInfo<D>> objects
	 * @return std::map<int, Interface<D>> a map of interface id to Interface oject
	 */
	static std::map<int, std::shared_ptr<Interface<D>>>
	EnumerateIfacesFromPiinfoVector(std::vector<std::shared_ptr<const PatchIfaceInfo<D>>> piinfos)
	{
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		std::map<int, std::shared_ptr<Interface<D>>> ifaces;
		auto                                         insert_patch
		= [&](int id, Side<2> s, std::shared_ptr<const PatchIfaceInfo<D>> piinfo) {
			  std::shared_ptr<Interface<D>> &iface_ptr = ifaces[id];
			  if (iface_ptr == nullptr) {
				  iface_ptr.reset(new Interface<D>(id));
			  }
			  iface_ptr->insert(s, piinfo);
		  };
		std::set<int>                              incoming_procs;
		std::map<int, std::map<int, Interface<D>>> off_proc_ifaces;
		auto                                       insert_off_proc_patch
		= [&](int rank, int id, Side<2> s, std::shared_ptr<const PatchIfaceInfo<D>> piinfo) {
			  off_proc_ifaces[rank].emplace(id, id);
			  off_proc_ifaces[rank].at(id).insert(s, piinfo);
		  };
		for (auto piinfo : piinfos) {
			for (Side<D> s : Side<D>::getValues()) {
				if (piinfo->pinfo->hasNbr(s)) {
					switch (piinfo->pinfo->getNbrType(s)) {
						case NbrType::Normal: {
							auto info = piinfo->getNormalIfaceInfo(s);
							if (info->rank == rank) {
								insert_patch(info->id, s, piinfo);
								if (info->nbr_info->rank != rank) {
									incoming_procs.insert(info->nbr_info->rank);
								}
							} else {
								insert_off_proc_patch(info->rank, info->id, s, piinfo);
							}
						} break;
						case NbrType::Fine: {
							auto info = piinfo->getFineIfaceInfo(s);
							insert_patch(info->id, s, piinfo);
							for (size_t i = 0; i < Orthant<D - 1>::num_orthants; i++) {
								if (info->fine_ranks[i] == rank) {
									insert_patch(info->fine_ids[i], s, piinfo);
								} else {
									insert_off_proc_patch(info->fine_ranks[i], info->fine_ids[i], s,
									                      piinfo);
									incoming_procs.insert(info->fine_ranks[i]);
								}
							}
						} break;
						case NbrType::Coarse: {
							auto info = piinfo->getCoarseIfaceInfo(s);
							insert_patch(info->id, s, piinfo);
							if (info->coarse_rank == rank) {
								insert_patch(info->coarse_id, s, piinfo);
							} else {
								insert_off_proc_patch(info->coarse_rank, info->coarse_id, s,
								                      piinfo);
								incoming_procs.insert(info->coarse_rank);
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
				Interface<D> &iface = q.second;
				size += iface.serialize(nullptr);
			}
			char *buffer = new char[size];
			buffers.push_back(buffer);
			int pos = 0;
			for (auto q : p.second) {
				Interface<D> &iface = q.second;
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
				Interface<D> ifs;
				reader >> ifs;
				ifaces.at(ifs.id)->merge(ifs);
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
};
} // namespace Schur
} // namespace ThunderEgg
#endif
