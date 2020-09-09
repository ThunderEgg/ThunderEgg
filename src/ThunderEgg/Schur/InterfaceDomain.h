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

#ifndef THUNDEREGG_SCHUR_InterfaceDomain_H
#define THUNDEREGG_SCHUR_InterfaceDomain_H
#include <ThunderEgg/Domain.h>
#include <ThunderEgg/Schur/Interface.h>
#include <ThunderEgg/Schur/PatchIfaceInfo.h>
#include <ThunderEgg/ValVector.h>
#include <ThunderEgg/VectorGenerator.h>
#include <deque>
#include <memory>
#include <valarray>
namespace ThunderEgg
{
namespace Schur
{
/**
 * @brief Represents the Schur compliment domain of the problem.
 *
 * This class mainly manages a set of interface that makes up the Schur compliment system. It is
 * responsible for setting up the indexing of the interfaces, which is used in the rest of the
 * ThunderEgg library.
 *
 * @tparam D the number of Cartesian dimensions
 */
template <size_t D> class InterfaceDomain
{
	private:
	std::shared_ptr<Domain<D>> domain;

	/**
	 * @brief Vector of PatchIfaceInfo pointers where index in the vector corresponds to the patch's
	 * local index
	 */
	std::vector<std::shared_ptr<PatchIfaceInfo<D>>>   sinfo_vector;
	std::map<int, std::shared_ptr<PatchIfaceInfo<D>>> id_sinfo_map;
	std::map<int, Interface<D>>                       ifaces;

	std::vector<int>                       id_map_vec;
	std::vector<int>                       global_map_vec;
	int                                    ghost_start;
	int                                    matrix_extra_ghost_start;
	int                                    rank;
	std::vector<std::tuple<int, int, int>> matrix_out_id_local_rank_vec;
	std::vector<std::tuple<int, int, int>> matrix_in_id_local_rank_vec;
	std::vector<std::tuple<int, int, int>> patch_out_id_local_rank_vec;
	std::vector<std::tuple<int, int, int>> patch_in_id_local_rank_vec;
	void                                   indexIfacesLocal()
	{
		using namespace std;
		int curr_i = 0;
		id_map_vec.reserve(ifaces.size());
		id_map_vec.resize(0);
		map<int, int>       rev_map;
		set<int>            enqueued;
		set<pair<int, int>> out_matrix_offs;
		set<pair<int, int>> in_matrix_offs;
		// set local indexes in schur compliment matrix
		if (!ifaces.empty()) {
			set<int> todo;
			for (auto &p : ifaces) {
				todo.insert(p.first);
			}
			while (!todo.empty()) {
				deque<int> queue;
				queue.push_back(*todo.begin());
				enqueued.insert(*todo.begin());
				while (!queue.empty()) {
					int id = queue.front();
					todo.erase(id);
					queue.pop_front();

					// set local index for interface
					id_map_vec.push_back(id);
					Interface<D> &ifs = ifaces.at(id);
					rev_map[id]       = curr_i;

					curr_i++;
					for (auto patch : ifs.patches) {
						if (patch.piinfo->getIfaceInfo(patch.side)->id == id) {
							// outgoing affects all ifaces
							for (Side<D> s : Side<D>::getValues()) {
								if (patch.piinfo->pinfo->hasNbr(s)) {
									switch (patch.piinfo->pinfo->getNbrType(s)) {
										case NbrType::Normal: {
											auto info = patch.piinfo->getNormalIfaceInfo(s);
											if (info->rank != rank) {
												in_matrix_offs.emplace(info->rank, info->id);
												out_matrix_offs.emplace(info->rank, id);
											}
										} break;
										case NbrType::Fine: {
											auto info = patch.piinfo->getFineIfaceInfo(s);
											if (info->rank != rank) {
												in_matrix_offs.emplace(info->rank, info->id);
												out_matrix_offs.emplace(info->rank, id);
											}
											for (size_t i = 0; i < Orthant<D - 1>::num_orthants;
											     i++) {
												if (info->fine_ranks[i] != rank) {
													out_matrix_offs.emplace(info->fine_ranks[i],
													                        id);
												}
											}
										} break;
										case NbrType::Coarse: {
											auto info = patch.piinfo->getCoarseIfaceInfo(s);
											if (info->rank != rank) {
												in_matrix_offs.emplace(info->rank, info->id);
												out_matrix_offs.emplace(info->rank, id);
											}
											if (info->coarse_rank != rank) {
												out_matrix_offs.emplace(info->coarse_rank, id);
											}
										} break;
									}
								}
							}
						} else {
							// outgoing affects only ifaces on iface_side
							// iface side has to handled specially
							switch (patch.piinfo->pinfo->getNbrType(patch.side)) {
								case NbrType::Normal: {
									auto info = patch.piinfo->getNormalIfaceInfo(patch.side);
									if (info->rank != rank) {
										in_matrix_offs.emplace(info->rank, info->id);
										out_matrix_offs.emplace(info->rank, id);
									}
								} break;
								case NbrType::Fine: {
									auto info = patch.piinfo->getFineIfaceInfo(patch.side);
									if (info->rank != rank) {
										in_matrix_offs.emplace(info->rank, info->id);
										out_matrix_offs.emplace(info->rank, id);
									}
									for (size_t i = 0; i < Orthant<D - 1>::num_orthants; i++) {
										if (info->fine_ranks[i] != rank) {
											out_matrix_offs.emplace(info->fine_ranks[i], id);
										}
									}
								} break;
								case NbrType::Coarse: {
									auto info = patch.piinfo->getCoarseIfaceInfo(patch.side);
									if (info->rank != rank) {
										in_matrix_offs.emplace(info->rank, info->id);
										out_matrix_offs.emplace(info->rank, id);
									}
									if (info->coarse_rank != rank) {
										out_matrix_offs.emplace(info->coarse_rank, id);
									}
								} break;
							}
							// handle the rest of the cases
							for (Side<D> s : Side<D>::getValues()) {
								if (s != patch.side && patch.piinfo->pinfo->hasNbr(s)) {
									switch (patch.piinfo->pinfo->getNbrType(s)) {
										case NbrType::Normal: {
											auto info = patch.piinfo->getNormalIfaceInfo(s);
											if (info->rank != rank) {
												in_matrix_offs.emplace(info->rank, info->id);
											}
										} break;
										case NbrType::Fine: {
											auto info = patch.piinfo->getFineIfaceInfo(s);
											if (info->rank != rank) {
												in_matrix_offs.emplace(info->rank, info->id);
											}
										} break;
										case NbrType::Coarse: {
											auto info = patch.piinfo->getCoarseIfaceInfo(s);
											if (info->rank != rank) {
												in_matrix_offs.emplace(info->rank, info->id);
											}
										} break;
									}
								}
							}
						}
					}
					/*
					for (int nbr_id : ifs.getNbrs()) {
					    if (ifaces.count(nbr_id) && !enqueued.count(nbr_id)) {
					        enqueued.insert(nbr_id);
					        queue.push_back(nbr_id);
					    }
					}
					*/
				}
			}
		}
		set<pair<int, int>> out_patch_offs;
		set<pair<int, int>> in_patch_offs;
		ghost_start = curr_i;
		for (shared_ptr<PatchIfaceInfo<D>> &sinfo : sinfo_vector) {
			for (Side<D> s : Side<D>::getValues()) {
				if (sinfo->pinfo->hasNbr(s)) {
					int id = sinfo->getIfaceInfo(s)->id;
					if (!rev_map.count(id)) {
						id_map_vec.push_back(id);
						rev_map[id] = curr_i;
						curr_i++;
					}
					switch (sinfo->pinfo->getNbrType(s)) {
						case NbrType::Normal: {
							auto info = sinfo->getNormalIfaceInfo(s);
							if (info->nbr_info->rank != rank) {
								in_patch_offs.insert(make_pair(info->nbr_info->rank, id));
								out_patch_offs.insert(make_pair(info->nbr_info->rank, id));
							}
						} break;
						case NbrType::Coarse: {
							auto info = sinfo->getCoarseIfaceInfo(s);
							if (info->nbr_info->rank != rank) {
								in_patch_offs.insert(make_pair(info->nbr_info->rank, id));
								out_patch_offs.insert(
								make_pair(info->nbr_info->rank, info->coarse_id));
								if (!rev_map.count(info->coarse_id)) {
									id_map_vec.push_back(info->coarse_id);
									rev_map[info->coarse_id] = curr_i;
									curr_i++;
								}
							}
						} break;
						case NbrType::Fine: {
							auto info = sinfo->getFineIfaceInfo(s);
							for (Orthant<D - 1> o : Orthant<D - 1>::getValues()) {
								if (info->nbr_info->ranks[o.getIndex()] != rank) {
									in_patch_offs.insert(
									make_pair(info->nbr_info->ranks[o.getIndex()], id));
									out_patch_offs.insert(
									make_pair(info->nbr_info->ranks[o.getIndex()],
									          info->fine_ids[o.getIndex()]));
									if (!rev_map.count(info->fine_ids[o.getIndex()])) {
										id_map_vec.push_back(info->fine_ids[o.getIndex()]);
										rev_map[info->fine_ids[o.getIndex()]] = curr_i;
										curr_i++;
									}
								}
							}
						} break;
					}
				}
			}
		}
		matrix_extra_ghost_start = curr_i;
		// get off proc data
		for (auto pair : out_patch_offs) {
			int rank = pair.first;
			int id   = pair.second;
			patch_out_id_local_rank_vec.emplace_back(id, rev_map[id], rank);
		}
		for (auto pair : in_patch_offs) {
			int rank = pair.first;
			int id   = pair.second;
			patch_in_id_local_rank_vec.emplace_back(id, rev_map[id], rank);
		}
		for (auto pair : out_matrix_offs) {
			int rank = pair.first;
			int id   = pair.second;
			matrix_out_id_local_rank_vec.emplace_back(id, rev_map[id], rank);
		}
		for (auto pair : in_matrix_offs) {
			int rank = pair.first;
			int id   = pair.second;
			if (!rev_map.count(id)) {
				id_map_vec.push_back(id);
				rev_map[id] = curr_i;
				curr_i++;
			}
			matrix_in_id_local_rank_vec.emplace_back(id, rev_map[id], rank);
		}
		// TODO rest of matrix
		for (auto &p : ifaces) {
			// p.second.setLocalIndexes(rev_map);
		}
		for (auto &sinfo : sinfo_vector) {
			sinfo->setLocalIndexesFromId(rev_map);
		}
		indexIfacesGlobal();
	}
	void indexDomainIfacesLocal() {}
	void indexIfacesGlobal()
	{
		using namespace std;
		// global indices are going to be sequentially increasing with rank
		int local_size = ifaces.size();
		int start_i;
		MPI_Scan(&local_size, &start_i, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		start_i -= local_size;
		vector<int> new_global(id_map_vec.size());
		iota(new_global.begin(), new_global.begin() + ghost_start, start_i);

		// send outgoing messages
		vector<MPI_Request> requests;
		requests.reserve(matrix_in_id_local_rank_vec.size() + matrix_out_id_local_rank_vec.size());

		// recv info
		for (const auto &tuple : matrix_in_id_local_rank_vec) {
			int         id          = std::get<0>(tuple);
			int         local_index = std::get<1>(tuple);
			int         source      = std::get<2>(tuple);
			MPI_Request request;
			MPI_Irecv(&new_global[local_index], 1, MPI_INT, source, id, MPI_COMM_WORLD, &request);
			requests.push_back(request);
		}
		// send info
		for (const auto &tuple : matrix_out_id_local_rank_vec) {
			int         id          = std::get<0>(tuple);
			int         local_index = std::get<1>(tuple);
			int         dest        = std::get<2>(tuple);
			MPI_Request request;
			MPI_Isend(&new_global[local_index], 1, MPI_INT, dest, id, MPI_COMM_WORLD, &request);
			requests.push_back(request);
		}
		// wait for all
		MPI_Waitall(requests.size(), &requests[0], MPI_STATUSES_IGNORE);
		global_map_vec = new_global;
	}

	int                    num_global_ifaces = 0;
	int                    iface_stride;
	std::array<int, D - 1> lengths;

	public:
	InterfaceDomain() = default;
	/**
	 * @brief Create a InterfaceDomain from a given DomainCollection
	 *
	 * @param domain the DomainCollection
	 * @param comm the teuchos communicator
	 */
	explicit InterfaceDomain(std::shared_ptr<Domain<D>> domain)
	{
		this->domain = domain;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		iface_stride = 1;
		for (size_t i = 0; i < D - 1; i++) {
			iface_stride *= domain->getNs()[i];
			lengths[i] = domain->getNs()[i];
		}
		sinfo_vector.reserve(domain->getNumLocalPatches());
		for (auto &pinfo : domain->getPatchInfoVector()) {
			sinfo_vector.emplace_back(new PatchIfaceInfo<D>(pinfo));
			id_sinfo_map[pinfo->id] = sinfo_vector.back();
		}
		// ifaces = Interface<D>::EnumerateIfaces(sinfo_vector.begin(), sinfo_vector.end());
		indexDomainIfacesLocal();
		indexIfacesLocal();
		int num_ifaces = ifaces.size();
		MPI_Allreduce(&num_ifaces, &num_global_ifaces, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	}

	void updateInterfaceDist(std::shared_ptr<Vector<D - 1>> gamma)
	{
		// send outgoing messages
		std::vector<MPI_Request> requests;
		requests.reserve(patch_in_id_local_rank_vec.size() + patch_out_id_local_rank_vec.size());

		// send info
		for (const auto &tuple : patch_out_id_local_rank_vec) {
			int         id          = std::get<0>(tuple);
			int         local_index = std::get<1>(tuple);
			int         dest        = std::get<2>(tuple);
			MPI_Request request;
			MPI_Isend(gamma->getLocalData(local_index).getPtr(), iface_stride, MPI_DOUBLE, dest, id,
			          MPI_COMM_WORLD, &request);
			requests.push_back(request);
		}
		// recv info
		for (const auto &tuple : patch_in_id_local_rank_vec) {
			int id          = std::get<0>(tuple);
			int local_index = std::get<1>(tuple);
			int source      = std::get<2>(tuple);

			double buffer[iface_stride];
			MPI_Recv(buffer, iface_stride, MPI_DOUBLE, source, id, MPI_COMM_WORLD,
			         MPI_STATUS_IGNORE);
			double *data = gamma->getLocalData(local_index).getPtr();
			for (int i = 0; i < iface_stride; i++) {
				data[i] += buffer[i];
			}
		}

		// wait for all
		MPI_Waitall(requests.size(), &requests[0], MPI_STATUSES_IGNORE);
	}
	void scatterInterfaceDist(std::shared_ptr<const Vector<D - 1>> gamma)
	{
		// send outgoing messages
		std::vector<MPI_Request> requests;
		requests.reserve(patch_in_id_local_rank_vec.size() + patch_out_id_local_rank_vec.size());

		// send info
		for (const auto &tuple : patch_out_id_local_rank_vec) {
			int id          = std::get<0>(tuple);
			int local_index = std::get<1>(tuple);
			int dest        = std::get<2>(tuple);
			if (local_index < ghost_start) {
				MPI_Request request;
				MPI_Isend(gamma->getLocalData(local_index).getPtr(), iface_stride, MPI_DOUBLE, dest,
				          id, MPI_COMM_WORLD, &request);
				requests.push_back(request);
			}
		}
		// recv info
		for (const auto &tuple : patch_in_id_local_rank_vec) {
			int id          = std::get<0>(tuple);
			int local_index = std::get<1>(tuple);
			int source      = std::get<2>(tuple);

			if (local_index >= ghost_start) {
				MPI_Recv(gamma->getLocalData(local_index).getPtr(), iface_stride, MPI_DOUBLE,
				         source, id, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
		}

		// wait for all
		MPI_Waitall(requests.size(), &requests[0], MPI_STATUSES_IGNORE);
	}
	std::shared_ptr<ValVector<D - 1>> getNewGlobalInterfaceVector()
	{
		return std::shared_ptr<ValVector<D - 1>>(
		new ValVector<D - 1>(MPI_COMM_WORLD, lengths, 0, matrix_extra_ghost_start));
	}
	/**
	 * @brief Get the number of local Interfaces on this rank
	 *
	 * @return int the number of local Interfaces
	 */
	int getNumLocalInterfaces() const
	{
		return 0;
	}
	/**
	 * @brief Get the number of Interfaces on all ranks
	 *
	 * @return int the number of Interfaces on all ranks
	 */
	int getNumGlobalInterfaces() const
	{
		return 0;
	}
	/**
	 * @brief Get the vector Interfaces objects for this rank
	 *
	 * The location of each Interface in the vector will coorespond to the Interafce's local index
	 *
	 * @return const std::vector<std::shared_ptr<const Interface<D>>> the vector of Interface
	 * objects
	 */
	const std::vector<std::shared_ptr<const Interface<D>>> getInterfaces() const
	{
		return std::vector<std::shared_ptr<const Interface<D>>>();
	}
	/**
	 * @brief Get the vector PatchIfaceInfo objects for this rank
	 *
	 * The location of each PatchIfaceInfo in the vector will coorespond to the patch's local index
	 *
	 * @return const std::vector<std::shared_ptr<const PatchIfaceInfo<D>>> the vector of
	 * PatchIfaceInfo objects
	 */
	const std::vector<std::shared_ptr<const PatchIfaceInfo<D>>> getPatchIfaceInfos()
	{
		return std::vector<std::shared_ptr<const PatchIfaceInfo<D>>>();
	}
	/**
	 * @brief Get the Domain object that cooresponds to this InterfaceDomain
	 *
	 * @return std::shared_ptr<Domain<D>> the Domain object
	 */
	std::shared_ptr<Domain<D>> getDomain()
	{
		return domain;
	}
};
template <size_t D> class InterfaceDomainVG : public VectorGenerator<D>
{
	private:
	std::shared_ptr<InterfaceDomain<D + 1>> sh;

	public:
	InterfaceDomainVG(std::shared_ptr<InterfaceDomain<D + 1>> sh)
	{
		this->sh = sh;
	}
	std::shared_ptr<Vector<D>> getNewVector()
	{
		return sh->getNewSchurVec();
	}
};
extern template class InterfaceDomain<2>;
extern template class InterfaceDomain<3>;
} // namespace Schur
} // namespace ThunderEgg
#endif
