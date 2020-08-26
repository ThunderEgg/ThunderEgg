/***************************************************************************
 *  ThunderEgg, a library for solving Poisson's equation on adaptively
 *  refined block-structured Cartesian grids
 *
 *  Copyright (C) 2019-2020 ThunderEgg Developers. See AUTHORS.md file at the
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
#include "../tpl/json.hpp"
#include <ThunderEgg/Domain.h>
#include <fstream>
#include <string>
template <size_t D> class DomainReader
{
	private:
	bool                                      neumann;
	std::array<int, D>                        ns;
	int                                       num_ghost;
	std::shared_ptr<ThunderEgg::Domain<D>>    coarser_domain;
	std::shared_ptr<ThunderEgg::Domain<D>>    finer_domain;
	std::shared_ptr<ThunderEgg::PatchInfo<D>> parsePatch(nlohmann::json &patch_j)
	{
		auto pinfo             = std::make_shared<ThunderEgg::PatchInfo<D>>();
		pinfo->num_ghost_cells = num_ghost;
		patch_j.at("id").get_to(pinfo->id);
		patch_j.at("rank").get_to(pinfo->rank);
		patch_j.at("rank").get_to(pinfo->rank);
		if (patch_j.find("child_ids") != patch_j.end()) {
			patch_j.at("child_ids").get_to(pinfo->child_ids);
			patch_j.at("child_ranks").get_to(pinfo->child_ranks);
		}
		if (patch_j.find("parent_id") != patch_j.end()) {
			patch_j.at("parent_id").get_to(pinfo->parent_id);
			patch_j.at("parent_rank").get_to(pinfo->parent_rank);
			if (patch_j.contains("orth_on_parent")) {
				pinfo->orth_on_parent = GetOrthant(patch_j.at("orth_on_parent"));
			}
		}
		patch_j.at("starts").get_to(pinfo->starts);
		patch_j.at("lengths").get_to(pinfo->spacings);
		for (size_t d = 0; d < D; d++) {
			pinfo->spacings[d] /= ns[d];
		}

		pinfo->ns              = ns;
		nlohmann::json &nbrs_j = patch_j.at("nbrs");
		for (nlohmann::json &nbr_j : nbrs_j) {
			std::string         nbr_type_str = nbr_j["type"];
			ThunderEgg::NbrType nbr_type;
			if (nbr_type_str == "NORMAL") {
				nbr_type = ThunderEgg::NbrType::Normal;
			} else if (nbr_type_str == "COARSE") {
				nbr_type = ThunderEgg::NbrType::Coarse;
			} else if (nbr_type_str == "FINE") {
				nbr_type = ThunderEgg::NbrType::Fine;
			} else {
				throw "parsing error";
			}
			std::string         side_str = nbr_j["side"];
			ThunderEgg::Side<D> side;
			if (side_str == "WEST") {
				side = ThunderEgg::Side<D>::west();
			} else if (side_str == "EAST") {
				side = ThunderEgg::Side<D>::east();
			} else if (side_str == "SOUTH") {
				side = ThunderEgg::Side<D>::south();
			} else if (side_str == "NORTH") {
				side = ThunderEgg::Side<D>::north();
			} else if (side_str == "BOTTOM") {
				side = ThunderEgg::Side<D>::bottom();
			} else if (side_str == "TOP") {
				side = ThunderEgg::Side<D>::top();
			} else {
				throw "parsing error";
			}
			using array  = std::array<int, ThunderEgg::Orthant<D - 1>::num_orthants>;
			using array1 = std::array<int, 1>;
			switch (nbr_type) {
				case ThunderEgg::NbrType::Normal: {
					pinfo->nbr_info[side.getIndex()].reset(
					new ThunderEgg::NormalNbrInfo<D>(nbr_j.at("ids").get<array1>()[0]));
					pinfo->getNormalNbrInfo(side).rank = nbr_j.at("ranks").get<array1>()[0];
				} break;
				case ThunderEgg::NbrType::Fine: {
					pinfo->nbr_info[side.getIndex()].reset(
					new ThunderEgg::FineNbrInfo<D>(nbr_j.at("ids").get<array>()));
					pinfo->getFineNbrInfo(side).ranks = nbr_j.at("ranks").get<array>();
				} break;
				case ThunderEgg::NbrType::Coarse: {
					pinfo->nbr_info[side.getIndex()].reset(new ThunderEgg::CoarseNbrInfo<D>(
					nbr_j.at("ids").get<array1>()[0],
					GetOrthant(nbr_j.at("orth_on_coarse")).collapseOnAxis(side.getAxisIndex())));
					pinfo->getCoarseNbrInfo(side).rank = nbr_j.at("ranks").get<array1>()[0];
				} break;
				default:
					throw "parsing error";
			}
		}
		for (ThunderEgg::Side<D> s : ThunderEgg::Side<D>::getValues()) {
			if (!pinfo->hasNbr(s)) {
				pinfo->neumann[s.getIndex()] = neumann;
			}
		}
		return pinfo;
	}

	static ThunderEgg::Orthant<D> GetOrthant(nlohmann::json &orth_j);

	public:
	DomainReader(std::string file_name, std::array<int, D> ns_in, int num_ghost_in,
	             bool neumann_in = false)
	: neumann(neumann_in), ns(ns_in), num_ghost(num_ghost_in)
	{
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		nlohmann::json j;
		std::ifstream  input_stream(file_name);
		if (!input_stream.good()) {
			throw "could not open file";
		}
		input_stream >> j;
		input_stream.close();
		std::map<int, std::shared_ptr<ThunderEgg::PatchInfo<D>>> finer_map;
		for (nlohmann::json &patch_j : j.at("levels")[0]) {
			auto patch = parsePatch(patch_j);
			if (patch->rank == rank)
				finer_map[patch->id] = patch;
		}
		finer_domain = std::make_shared<ThunderEgg::Domain<D>>(finer_map, ns, num_ghost);
		std::map<int, std::shared_ptr<ThunderEgg::PatchInfo<D>>> coarser_map;
		for (nlohmann::json &patch_j : j.at("levels")[1]) {
			auto patch = parsePatch(patch_j);
			if (patch->rank == rank)
				coarser_map[patch->id] = patch;
		}
		coarser_domain = std::make_shared<ThunderEgg::Domain<D>>(coarser_map, ns, num_ghost);
	}
	std::shared_ptr<ThunderEgg::Domain<D>> getCoarserDomain()
	{
		return coarser_domain;
	}
	std::shared_ptr<ThunderEgg::Domain<D>> getFinerDomain()
	{
		return finer_domain;
	}
};
template <> ThunderEgg::Orthant<2> DomainReader<2>::GetOrthant(nlohmann::json &orth_j);
template <> ThunderEgg::Orthant<3> DomainReader<3>::GetOrthant(nlohmann::json &orth_j);
extern template class DomainReader<2>;