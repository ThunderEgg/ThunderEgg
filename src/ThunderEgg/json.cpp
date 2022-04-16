/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured
 *  Cartesian grids.
 *
 *  Copyright (c) 2020-2021 Scott Aiton
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

#include <ThunderEgg/Domain.h>
#include <ThunderEgg/tpl/json.hpp>

namespace ThunderEgg {

void
to_json(tpl::nlohmann::json& j, const NbrType& o)
{
  switch (o) {
    case NbrType::Normal:
      j = "NORMAL";
      break;
    case NbrType::Coarse:
      j = "COARSE";
      break;
    case NbrType::Fine:
      j = "FINE";
      break;
    default:
      j = nullptr;
      break;
  }
}
void
from_json(const tpl::nlohmann::json& j, NbrType& o)
{
  if (j.is_string()) {
    if (j.get<std::string>() == "NORMAL") {
      o = NbrType::Normal;
    } else if (j.get<std::string>() == "COARSE") {
      o = NbrType::Coarse;
    } else if (j.get<std::string>() == "FINE") {
      o = NbrType::Fine;
    } else {
      throw RuntimeError("unsuppored NbrType");
    }
  } else {
    throw RuntimeError("unsuppored NbrType");
  }
}

template<int D>
void
to_json(tpl::nlohmann::json& j, const NormalNbrInfo<D>& n)
{
  j["type"] = NbrType::Normal;
  j["ids"] = { n.id };
  j["ranks"] = { n.rank };
}

template<int D>
void
from_json(const tpl::nlohmann::json& j, NormalNbrInfo<D>& n)
{
  n.id = j["ids"][0];
  n.rank = j["ranks"][0];
}

template void
to_json(tpl::nlohmann::json& j, const NormalNbrInfo<0>& n);
template void
to_json(tpl::nlohmann::json& j, const NormalNbrInfo<1>& n);
template void
to_json(tpl::nlohmann::json& j, const NormalNbrInfo<2>& n);

template void
from_json(const tpl::nlohmann::json& j, NormalNbrInfo<0>& n);
template void
from_json(const tpl::nlohmann::json& j, NormalNbrInfo<1>& n);
template void
from_json(const tpl::nlohmann::json& j, NormalNbrInfo<2>& n);

template<int D>
void
to_json(tpl::nlohmann::json& j, const CoarseNbrInfo<D>& n)
{
  j["type"] = NbrType::Coarse;
  j["ids"] = { n.id };
  j["ranks"] = { n.rank };
  j["orth_on_coarse"] = n.orth_on_coarse;
}

template<int D>
void
from_json(const tpl::nlohmann::json& j, CoarseNbrInfo<D>& n)
{
  n.id = j["ids"][0];
  n.rank = j["ranks"][0];
  if (j.contains("orth_on_coarse")) {
    n.orth_on_coarse = j["orth_on_coarse"].get<Orthant<D>>();
  } else {
    n.orth_on_coarse = Orthant<D>::null();
  }
}

template void
to_json(tpl::nlohmann::json& j, const CoarseNbrInfo<0>& n);
template void
to_json(tpl::nlohmann::json& j, const CoarseNbrInfo<1>& n);
template void
to_json(tpl::nlohmann::json& j, const CoarseNbrInfo<2>& n);

template void
from_json(const tpl::nlohmann::json& j, CoarseNbrInfo<0>& n);
template void
from_json(const tpl::nlohmann::json& j, CoarseNbrInfo<1>& n);
template void
from_json(const tpl::nlohmann::json& j, CoarseNbrInfo<2>& n);

template<int D>
void
to_json(tpl::nlohmann::json& j, const FineNbrInfo<D>& n)
{
  j["type"] = NbrType::Fine;
  j["ids"] = n.ids;
  j["ranks"] = n.ranks;
}

template<int D>
void
from_json(const tpl::nlohmann::json& j, FineNbrInfo<D>& n)
{
  n.ids = j["ids"].get<std::array<int, Orthant<D>::num_orthants>>();
  n.ranks = j["ranks"].get<std::array<int, Orthant<D>::num_orthants>>();
}

template void
to_json(tpl::nlohmann::json& j, const FineNbrInfo<0>& n);
template void
to_json(tpl::nlohmann::json& j, const FineNbrInfo<1>& n);
template void
to_json(tpl::nlohmann::json& j, const FineNbrInfo<2>& n);

template void
from_json(const tpl::nlohmann::json& j, FineNbrInfo<0>& n);
template void
from_json(const tpl::nlohmann::json& j, FineNbrInfo<1>& n);
template void
from_json(const tpl::nlohmann::json& j, FineNbrInfo<2>& n);

template<int D>
void
to_json(tpl::nlohmann::json& j, const PatchInfo<D>& pinfo)
{
  j["id"] = pinfo.id;
  j["parent_id"] = pinfo.parent_id;
  j["parent_rank"] = pinfo.parent_rank;
  j["orth_on_parent"] = pinfo.orth_on_parent;
  j["rank"] = pinfo.rank;
  j["starts"] = pinfo.starts;
  j["lengths"] = pinfo.spacings;
  j["refine_level"] = pinfo.refine_level;
  for (int i = 0; i < D; i++) {
    j["lengths"][i] = pinfo.spacings[i] * pinfo.ns[i];
  }
  j["nbrs"] = tpl::nlohmann::json::array();
  for (Side<D> s : Side<D>::getValues()) {
    if (pinfo.hasNbr(s)) {
      switch (pinfo.getNbrType(s)) {
        case NbrType::Normal:
          j["nbrs"].push_back(pinfo.getNormalNbrInfo(s));
          break;
        case NbrType::Coarse:
          j["nbrs"].push_back(pinfo.getCoarseNbrInfo(s));
          break;
        case NbrType::Fine:
          j["nbrs"].push_back(pinfo.getFineNbrInfo(s));
          break;
        default:
          throw RuntimeError("Unsupported NbrType");
      }
      j["nbrs"].back()["type"] = pinfo.getNbrType(s);
      j["nbrs"].back()["side"] = s;
    }
  }
  if constexpr (D == 3) {
    j["edge_nbrs"] = tpl::nlohmann::json::array();
    for (Edge e : Edge::getValues()) {
      if (pinfo.hasNbr(e)) {
        switch (pinfo.getNbrType(e)) {
          case NbrType::Normal:
            j["edge_nbrs"].push_back(pinfo.getNormalNbrInfo(e));
            break;
          case NbrType::Coarse:
            j["edge_nbrs"].push_back(pinfo.getCoarseNbrInfo(e));
            break;
          case NbrType::Fine:
            j["edge_nbrs"].push_back(pinfo.getFineNbrInfo(e));
            break;
          default:
            throw RuntimeError("Unsupported NbrType");
        }
        j["edge_nbrs"].back()["type"] = pinfo.getNbrType(e);
        j["edge_nbrs"].back()["edge"] = e;
      }
    }
  }
  if constexpr (D >= 2) {
    j["corner_nbrs"] = tpl::nlohmann::json::array();
    for (Corner<D> c : Corner<D>::getValues()) {
      if (pinfo.hasNbr(c)) {
        switch (pinfo.getNbrType(c)) {
          case NbrType::Normal:
            j["corner_nbrs"].push_back(pinfo.getNormalNbrInfo(c));
            break;
          case NbrType::Coarse:
            j["corner_nbrs"].push_back(pinfo.getCoarseNbrInfo(c));
            break;
          case NbrType::Fine:
            j["corner_nbrs"].push_back(pinfo.getFineNbrInfo(c));
            break;
          default:
            throw RuntimeError("Unsupported NbrType");
        }
        j["corner_nbrs"].back()["type"] = pinfo.getNbrType(c);
        j["corner_nbrs"].back()["corner"] = c;
      }
    }
  }
  if (pinfo.child_ids[0] != -1) {
    j["child_ids"] = pinfo.child_ids;
    j["child_ranks"] = pinfo.child_ranks;
  }
}
template<int D>
void
from_json(const tpl::nlohmann::json& j, PatchInfo<D>& pinfo)
{
  pinfo.id = j["id"];
  pinfo.parent_id = j["parent_id"];
  pinfo.parent_rank = j["parent_rank"];
  if (j.contains("orth_on_parent")) {
    j["orth_on_parent"].get_to(pinfo.orth_on_parent);
  }
  if (j.contains("refine_level")) {
    pinfo.refine_level = j["refine_level"];
  }
  pinfo.rank = j["rank"];
  pinfo.starts = j["starts"].get<std::array<double, D>>();
  pinfo.spacings = j["lengths"].get<std::array<double, D>>();
  pinfo.ns.fill(1);
  for (const auto& nbr_j : j["nbrs"]) {
    Side<D> s = nbr_j["side"].get<Side<D>>();
    switch (nbr_j["type"].get<NbrType>()) {
      case NbrType::Normal:
        pinfo.setNbrInfo(s, new NormalNbrInfo<D - 1>());
        pinfo.getNormalNbrInfo(s) = nbr_j.get<NormalNbrInfo<D - 1>>();
        break;
      case NbrType::Coarse:
        pinfo.setNbrInfo(s, new CoarseNbrInfo<D - 1>());
        pinfo.getCoarseNbrInfo(s) = nbr_j.get<CoarseNbrInfo<D - 1>>();
        break;
      case NbrType::Fine:
        pinfo.setNbrInfo(s, new FineNbrInfo<D - 1>());
        pinfo.getFineNbrInfo(s) = nbr_j.get<FineNbrInfo<D - 1>>();
        break;
      default:
        throw RuntimeError("Unsupported NbrType");
    }
  }
  if constexpr (D == 3) {
    if (j.contains("edge_nbrs")) {
      for (const auto& nbr_j : j["edge_nbrs"]) {
        Edge e = nbr_j["edge"].get<Edge>();
        switch (nbr_j["type"].get<NbrType>()) {
          case NbrType::Normal:
            pinfo.setNbrInfo(e, new NormalNbrInfo<1>());
            pinfo.getNormalNbrInfo(e) = nbr_j.get<NormalNbrInfo<1>>();
            break;
          case NbrType::Coarse:
            pinfo.setNbrInfo(e, new CoarseNbrInfo<1>());
            pinfo.getCoarseNbrInfo(e) = nbr_j.get<CoarseNbrInfo<1>>();
            break;
          case NbrType::Fine:
            pinfo.setNbrInfo(e, new FineNbrInfo<1>());
            pinfo.getFineNbrInfo(e) = nbr_j.get<FineNbrInfo<1>>();
            break;
          default:
            throw RuntimeError("Unsupported NbrType");
        }
      }
    }
  }
  if constexpr (D >= 2) {
    if (j.contains("corner_nbrs")) {
      for (const auto& nbr_j : j["corner_nbrs"]) {
        Corner<D> c = nbr_j["corner"].get<Corner<D>>();
        switch (nbr_j["type"].get<NbrType>()) {
          case NbrType::Normal:
            pinfo.setNbrInfo(c, new NormalNbrInfo<0>());
            pinfo.getNormalNbrInfo(c) = nbr_j.get<NormalNbrInfo<0>>();
            break;
          case NbrType::Coarse:
            pinfo.setNbrInfo(c, new CoarseNbrInfo<0>());
            pinfo.getCoarseNbrInfo(c) = nbr_j.get<CoarseNbrInfo<0>>();
            break;
          case NbrType::Fine:
            pinfo.setNbrInfo(c, new FineNbrInfo<0>());
            pinfo.getFineNbrInfo(c) = nbr_j.get<FineNbrInfo<0>>();
            break;
          default:
            throw RuntimeError("Unsupported NbrType");
        }
      }
    }
  }
  if (j.contains("child_ids")) {
    pinfo.child_ids = j["child_ids"].get<std::array<int, Orthant<D>::num_orthants>>();
    pinfo.child_ranks = j["child_ranks"].get<std::array<int, Orthant<D>::num_orthants>>();
  }
}

template void
to_json(tpl::nlohmann::json& j, const PatchInfo<2>& pinfo);
template void
to_json(tpl::nlohmann::json& j, const PatchInfo<3>& pinfo);

template void
from_json(const tpl::nlohmann::json& j, PatchInfo<2>& pinfo);
template void
from_json(const tpl::nlohmann::json& j, PatchInfo<3>& pinfo);

template<int D>
void
to_json(tpl::nlohmann::json& j, const Domain<D>& domain)
{
  for (auto pinfo : domain.getPatchInfoVector()) {
    j.push_back(pinfo);
  }
}

template void
ThunderEgg::to_json<2>(ThunderEgg::tpl::nlohmann::json& j, const Domain<2>& domain);
template void
ThunderEgg::to_json<3>(ThunderEgg::tpl::nlohmann::json& j, const Domain<3>& domain);

template<int D>
void
Domain<D>::setTimer(std::shared_ptr<Timer> timer) const
{
  this->timer = timer;
  timer->addDomain(id, *this);
}

template void
Domain<2>::setTimer(std::shared_ptr<Timer> timer) const;
template void
Domain<3>::setTimer(std::shared_ptr<Timer> timer) const;

} // namespace ThunderEgg