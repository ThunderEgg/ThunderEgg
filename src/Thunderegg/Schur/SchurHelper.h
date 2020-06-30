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

#ifndef THUNDEREGG_SCHUR_SCHURHELPER_H
#define THUNDEREGG_SCHUR_SCHURHELPER_H
#include <Thunderegg/Domain.h>
#include <Thunderegg/Schur/IfaceSet.h>
#include <Thunderegg/Schur/SchurInfo.h>
#include <Thunderegg/ValVector.h>
#include <deque>
#include <memory>
#include <valarray>
namespace Thunderegg
{
namespace Schur
{
/**
 * @brief Represents the Schur compliment domain of the problem.
 *
 * This class mainly manages a set of interface that makes up the Schur compliment system. It is
 * responsible for setting up the indexing of the interfaces, which is used in the rest of the
 * Thunderegg library.
 *
 * @tparam D the number of Cartesian dimensions
 */
template <size_t D> class SchurHelper
{
	private:
	std::shared_ptr<Domain<D>> domain;

	/**
	 * @brief Vector of SchurInfo pointers where index in the vector corresponds to the patch's
	 * local index
	 */
	std::vector<std::shared_ptr<SchurInfo<D>>>   sinfo_vector;
	std::map<int, std::shared_ptr<SchurInfo<D>>> id_sinfo_map;
	std::map<int, IfaceSet<D>>                   ifaces;

	std::vector<int>                       id_map_vec;
	std::vector<int>                       global_map_vec;
	int                                    ghost_start;
	int                                    matrix_extra_ghost_start;
	int                                    rank;
	std::vector<std::tuple<int, int, int>> matrix_out_id_local_rank_vec;
	std::vector<std::tuple<int, int, int>> matrix_in_id_local_rank_vec;
	std::vector<std::tuple<int, int, int>> patch_out_id_local_rank_vec;
	std::vector<std::tuple<int, int, int>> patch_in_id_local_rank_vec;
	void                                   indexIfacesLocal();
	void                                   indexDomainIfacesLocal();
	void                                   indexIfacesGlobal();

	int                    num_global_ifaces = 0;
	int                    iface_stride;
	std::array<int, D - 1> lengths;

	public:
	SchurHelper() = default;
	/**
	 * @brief Create a SchurHelper from a given DomainCollection
	 *
	 * @param domain the DomainCollection
	 * @param comm the teuchos communicator
	 */
	SchurHelper(std::shared_ptr<Domain<D>> domain);

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
	std::shared_ptr<ValVector<D - 1>> getNewSchurVec()
	{
		return std::shared_ptr<ValVector<D - 1>>(
		new ValVector<D - 1>(MPI_COMM_WORLD, lengths, 0, matrix_extra_ghost_start));
	}
	std::shared_ptr<ValVector<D - 1>> getNewSchurDistVec()
	{
		return std::shared_ptr<ValVector<D - 1>>(
		new ValVector<D - 1>(MPI_COMM_SELF, lengths, 0, matrix_extra_ghost_start));
	}

	int getSchurVecLocalSize() const
	{
		return ghost_start * iface_stride;
	}
	int getSchurVecGlobalSize() const
	{
		return num_global_ifaces * iface_stride;
	}
	// getters
	const std::map<int, IfaceSet<D>> getIfaces() const
	{
		return ifaces;
	}
	const std::array<int, D - 1> getLengths() const
	{
		return lengths;
	}
	const std::vector<std::shared_ptr<SchurInfo<D>>> getSchurInfoVector()
	{
		return sinfo_vector;
	}
	std::shared_ptr<Domain<D>> getDomain()
	{
		return domain;
	}
};
template <size_t D> inline SchurHelper<D>::SchurHelper(std::shared_ptr<Domain<D>> domain)
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
		sinfo_vector.emplace_back(new SchurInfo<D>(pinfo));
		id_sinfo_map[pinfo->id] = sinfo_vector.back();
	}
	ifaces = IfaceSet<D>::EnumerateIfaces(sinfo_vector.begin(), sinfo_vector.end());
	indexDomainIfacesLocal();
	indexIfacesLocal();
	int num_ifaces = ifaces.size();
	MPI_Allreduce(&num_ifaces, &num_global_ifaces, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
}
template <size_t D> inline void SchurHelper<D>::indexDomainIfacesLocal() {}
template <size_t D> inline void SchurHelper<D>::indexIfacesLocal()
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
				IfaceSet<D> &ifs = ifaces.at(id);
				rev_map[id]      = curr_i;

				curr_i++;
				for (size_t idx = 0; idx < ifs.sinfos.size(); idx++) {
					Side<D> iface_side = ifs.sides[idx];
					auto    sinfo      = ifs.sinfos[idx];
					if (sinfo->getIfaceInfoPtr(iface_side)->id == id) {
						// outgoing affects all ifaces
						for (Side<D> s : Side<D>::getValues()) {
							if (sinfo->pinfo->hasNbr(s)) {
								switch (sinfo->pinfo->getNbrType(s)) {
									case NbrType::Normal: {
										const NormalIfaceInfo<D> &info
										= sinfo->getNormalIfaceInfo(s);
										if (info.rank != rank) {
											in_matrix_offs.emplace(info.rank, info.id);
											out_matrix_offs.emplace(info.rank, id);
										}
									} break;
									case NbrType::Fine: {
										const FineIfaceInfo<D> &info = sinfo->getFineIfaceInfo(s);
										if (info.rank != rank) {
											in_matrix_offs.emplace(info.rank, info.id);
											out_matrix_offs.emplace(info.rank, id);
										}
										for (int i = 0; i < Orthant<D - 1>::num_orthants; i++) {
											if (info.fine_ranks[i] != rank) {
												out_matrix_offs.emplace(info.fine_ranks[i], id);
											}
										}
									} break;
									case NbrType::Coarse: {
										const CoarseIfaceInfo<D> &info
										= sinfo->getCoarseIfaceInfo(s);
										if (info.rank != rank) {
											in_matrix_offs.emplace(info.rank, info.id);
											out_matrix_offs.emplace(info.rank, id);
										}
										if (info.coarse_rank != rank) {
											out_matrix_offs.emplace(info.coarse_rank, id);
										}
									} break;
								}
							}
						}
					} else {
						// outgoing affects only ifaces on iface_side
						// iface side has to handled specially
						switch (sinfo->pinfo->getNbrType(iface_side)) {
							case NbrType::Normal: {
								const NormalIfaceInfo<D> &info
								= sinfo->getNormalIfaceInfo(iface_side);
								if (info.rank != rank) {
									in_matrix_offs.emplace(info.rank, info.id);
									out_matrix_offs.emplace(info.rank, id);
								}
							} break;
							case NbrType::Fine: {
								const FineIfaceInfo<D> &info = sinfo->getFineIfaceInfo(iface_side);
								if (info.rank != rank) {
									in_matrix_offs.emplace(info.rank, info.id);
									out_matrix_offs.emplace(info.rank, id);
								}
								for (int i = 0; i < Orthant<D - 1>::num_orthants; i++) {
									if (info.fine_ranks[i] != rank) {
										out_matrix_offs.emplace(info.fine_ranks[i], id);
									}
								}
							} break;
							case NbrType::Coarse: {
								const CoarseIfaceInfo<D> &info
								= sinfo->getCoarseIfaceInfo(iface_side);
								if (info.rank != rank) {
									in_matrix_offs.emplace(info.rank, info.id);
									out_matrix_offs.emplace(info.rank, id);
								}
								if (info.coarse_rank != rank) {
									out_matrix_offs.emplace(info.coarse_rank, id);
								}
							} break;
						}
						// handle the rest of the cases
						for (Side<D> s : Side<D>::getValues()) {
							if (s != iface_side && sinfo->pinfo->hasNbr(s)) {
								switch (sinfo->pinfo->getNbrType(s)) {
									case NbrType::Normal: {
										const NormalIfaceInfo<D> &info
										= sinfo->getNormalIfaceInfo(s);
										if (info.rank != rank) {
											in_matrix_offs.emplace(info.rank, info.id);
										}
									} break;
									case NbrType::Fine: {
										const FineIfaceInfo<D> &info = sinfo->getFineIfaceInfo(s);
										if (info.rank != rank) {
											in_matrix_offs.emplace(info.rank, info.id);
										}
									} break;
									case NbrType::Coarse: {
										const CoarseIfaceInfo<D> &info
										= sinfo->getCoarseIfaceInfo(s);
										if (info.rank != rank) {
											in_matrix_offs.emplace(info.rank, info.id);
										}
									} break;
								}
							}
						}
					}
				}
				for (int nbr_id : ifs.getNbrs()) {
					if (ifaces.count(nbr_id) && !enqueued.count(nbr_id)) {
						enqueued.insert(nbr_id);
						queue.push_back(nbr_id);
					}
				}
			}
		}
	}
	set<pair<int, int>> out_patch_offs;
	set<pair<int, int>> in_patch_offs;
	ghost_start = curr_i;
	for (shared_ptr<SchurInfo<D>> &sinfo : sinfo_vector) {
		for (Side<D> s : Side<D>::getValues()) {
			if (sinfo->pinfo->hasNbr(s)) {
				int id = sinfo->getIfaceInfoPtr(s)->id;
				if (!rev_map.count(id)) {
					id_map_vec.push_back(id);
					rev_map[id] = curr_i;
					curr_i++;
				}
				switch (sinfo->pinfo->getNbrType(s)) {
					case NbrType::Normal: {
						NormalIfaceInfo<D> &info = sinfo->getNormalIfaceInfo(s);
						if (info.nbr_info->ptr == nullptr) {
							in_patch_offs.insert(make_pair(info.nbr_info->rank, id));
							out_patch_offs.insert(make_pair(info.nbr_info->rank, id));
						}
					} break;
					case NbrType::Coarse: {
						CoarseIfaceInfo<D> &info = sinfo->getCoarseIfaceInfo(s);
						if (info.nbr_info->ptr == nullptr) {
							in_patch_offs.insert(make_pair(info.nbr_info->rank, id));
							out_patch_offs.insert(make_pair(info.nbr_info->rank, info.coarse_id));
							if (!rev_map.count(info.coarse_id)) {
								id_map_vec.push_back(info.coarse_id);
								rev_map[info.coarse_id] = curr_i;
								curr_i++;
							}
						}
					} break;
					case NbrType::Fine: {
						FineIfaceInfo<D> &info = sinfo->getFineIfaceInfo(s);
						for (Orthant<D - 1> o : Orthant<D - 1>::getValues()) {
							if (info.nbr_info->ptrs[o.toInt()] == nullptr) {
								in_patch_offs.insert(
								make_pair(info.nbr_info->ranks[o.toInt()], id));
								out_patch_offs.insert(make_pair(info.nbr_info->ranks[o.toInt()],
								                                info.fine_ids[o.toInt()]));
								if (!rev_map.count(info.fine_ids[o.toInt()])) {
									id_map_vec.push_back(info.fine_ids[o.toInt()]);
									rev_map[info.fine_ids[o.toInt()]] = curr_i;
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
		p.second.setLocalIndexes(rev_map);
	}
	for (auto &sinfo : sinfo_vector) {
		sinfo->setLocalIndexes(rev_map);
	}
	indexIfacesGlobal();
} // namespace Thunderegg
template <size_t D> inline void SchurHelper<D>::indexIfacesGlobal()
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
template <size_t D> class SchurHelperVG : public VectorGenerator<D>
{
	private:
	std::shared_ptr<SchurHelper<D + 1>> sh;

	public:
	SchurHelperVG(std::shared_ptr<SchurHelper<D + 1>> sh)
	{
		this->sh = sh;
	}
	std::shared_ptr<Vector<D>> getNewVector()
	{
		return sh->getNewSchurVec();
	}
};
extern template class SchurHelper<2>;
extern template class SchurHelper<3>;
} // namespace Schur
} // namespace Thunderegg
#endif
