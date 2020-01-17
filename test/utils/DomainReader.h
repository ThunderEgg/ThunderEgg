/***************************************************************************
 *  Thunderegg, a library for solving Poisson's equation on adaptively
 *  refined block-structured Cartesian grids
 *
 *  Copyright (C) 2019-2020 Thunderegg Developers. See AUTHORS.md file at the
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
#include <Thunderegg/Domain.h>
#include <fstream>
#include <string>
template <size_t D> class DomainReader
{
	private:
	std::array<int, 2>                        ns;
	int                                       num_ghost;
	std::shared_ptr<Thunderegg::Domain<2>>    coarser_domain;
	std::shared_ptr<Thunderegg::Domain<2>>    finer_domain;
	std::shared_ptr<Thunderegg::PatchInfo<2>> parsePatch(nlohmann::json &patch_j)
	{
		auto pinfo             = std::make_shared<Thunderegg::PatchInfo<2>>();
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
			Thunderegg::Orthant<2> orth(patch_j.at("orth_on_parent"));
			pinfo->orth_on_parent = orth;
		}
		patch_j.at("starts").get_to(pinfo->starts);
		patch_j.at("lengths").get_to(pinfo->spacings);
		for (size_t d = 0; d < 2; d++) {
			pinfo->spacings[d] /= ns[d];
		}

		pinfo->ns              = ns;
		nlohmann::json &nbrs_j = patch_j.at("nbrs");
		for (nlohmann::json &nbr_j : nbrs_j) {
			std::string         nbr_type_str = nbr_j["type"];
			Thunderegg::NbrType nbr_type;
			if (nbr_type_str == "normal") {
				nbr_type = Thunderegg::NbrType::Normal;
			} else if (nbr_type_str == "coarse") {
				nbr_type = Thunderegg::NbrType::Coarse;
			} else if (nbr_type_str == "fine") {
				nbr_type = Thunderegg::NbrType::Fine;
			} else {
				throw "parsing error";
			}
			std::string         side_str = nbr_j["side"];
			Thunderegg::Side<2> side;
			if (side_str == "west") {
				side = Thunderegg::Side<2>::west;
			} else if (side_str == "east") {
				side = Thunderegg::Side<2>::east;
			} else if (side_str == "south") {
				side = Thunderegg::Side<2>::south;
			} else if (side_str == "north") {
				side = Thunderegg::Side<2>::north;
			} else {
				throw "parsing error";
			}
			switch (nbr_type) {
				case Thunderegg::NbrType::Normal: {
					pinfo->nbr_info[side.toInt()].reset(
					new Thunderegg::NormalNbrInfo<2>(nbr_j.at("id")));
					pinfo->getNormalNbrInfo(side).updateRank(nbr_j.at("rank"));
				} break;
				default:
					throw "parsing error";
			}
		}
		return pinfo;
	}

	public:
	DomainReader(std::string file_name, std::array<int, D> ns_in, int num_ghost_in)
	: ns(ns_in), num_ghost(num_ghost_in)
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
		std::map<int, std::shared_ptr<Thunderegg::PatchInfo<2>>> finer_map;
		for (nlohmann::json &patch_j : j.at("finer")) {
			auto patch = parsePatch(patch_j);
			if (patch->rank == rank)
				finer_map[patch->id] = patch;
		}
		finer_domain = std::make_shared<Thunderegg::Domain<2>>(finer_map);
		std::map<int, std::shared_ptr<Thunderegg::PatchInfo<2>>> coarser_map;
		for (nlohmann::json &patch_j : j.at("coarser")) {
			auto patch = parsePatch(patch_j);
			if (patch->rank == rank)
				coarser_map[patch->id] = patch;
		}
		coarser_domain = std::make_shared<Thunderegg::Domain<2>>(coarser_map);
	}
	std::shared_ptr<Thunderegg::Domain<2>> getCoarserDomain()
	{
		return coarser_domain;
	}
	std::shared_ptr<Thunderegg::Domain<2>> getFinerDomain()
	{
		return finer_domain;
	}
};
extern template class DomainReader<2>;