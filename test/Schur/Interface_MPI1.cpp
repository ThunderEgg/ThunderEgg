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

#include "../utils/DomainReader.h"
#include <ThunderEgg/Schur/Interface.h>

#include <algorithm>

#include <doctest.h>

using namespace std;
using namespace ThunderEgg;

TEST_CASE("Schur::Interface id constructor")
{
  Schur::Interface<2> iface(2);
  CHECK(iface.id == 2);
  CHECK(iface.global_index == -1);
  CHECK(iface.local_index == -1);
  CHECK(iface.patches.empty());
}
TEST_CASE("Schur::Interface insert Normal interface")
{
  DomainReader<2> domain_reader("mesh_inputs/2d_uniform_1x2_mpi1.json", { 10, 10 }, 0);
  auto domain = domain_reader.getFinerDomain();
  vector<shared_ptr<Schur::PatchIfaceInfo<2>>> piinfos;
  for (auto patch : domain.getPatchInfoVector()) {
    piinfos.push_back(make_shared<Schur::PatchIfaceInfo<2>>(patch));
  }
  int id = piinfos[0]->getNormalIfaceInfo(Side<2>::east())->id;
  Schur::Interface<2> iface(id);

  iface.insert(Side<2>::east(), piinfos[0]);
  REQUIRE(iface.patches.size() == 1);
  CHECK(iface.patches[0].side == Side<2>::east());
  CHECK(iface.patches[0].type.isNormal());
  CHECK(iface.patches[0].piinfo == piinfos[0]);

  iface.insert(Side<2>::west(), piinfos[1]);
  REQUIRE(iface.patches.size() == 2);
  CHECK(iface.patches[1].side == Side<2>::west());
  CHECK(iface.patches[1].type.isNormal());
  CHECK(iface.patches[0].piinfo == piinfos[0]);
  CHECK(iface.patches[1].piinfo == piinfos[1]);
}
TEST_CASE("Schur::Interface insert Coarse interface")
{
  DomainReader<2> domain_reader("mesh_inputs/2d_refined_east_1x2_mpi1.json", { 10, 10 }, 0);
  auto domain = domain_reader.getFinerDomain();
  vector<shared_ptr<Schur::PatchIfaceInfo<2>>> piinfos;
  for (auto patch : domain.getPatchInfoVector()) {
    piinfos.push_back(make_shared<Schur::PatchIfaceInfo<2>>(patch));
  }
  int id = piinfos[0]->getFineIfaceInfo(Side<2>::east())->id;
  Schur::Interface<2> iface(id);

  iface.insert(Side<2>::east(), piinfos[0]);
  REQUIRE(iface.patches.size() == 1);
  CHECK(iface.patches[0].side == Side<2>::east());
  CHECK(iface.patches[0].type.isCoarseToCoarse());
  CHECK(iface.patches[0].piinfo == piinfos[0]);

  iface.insert(Side<2>::west(), piinfos[piinfos[0]->pinfo.getFineNbrInfo(Side<2>::east()).local_indexes[0]]);
  REQUIRE(iface.patches.size() == 2);
  CHECK(iface.patches[1].side == Side<2>::west());
  CHECK(iface.patches[1].type.isFineToCoarse());
  CHECK(iface.patches[1].type.getOrthant() == Orthant<1>::lower());
  CHECK(iface.patches[0].piinfo == piinfos[0]);
  CHECK(iface.patches[1].piinfo == piinfos[piinfos[0]->pinfo.getFineNbrInfo(Side<2>::east()).local_indexes[0]]);

  iface.insert(Side<2>::west(), piinfos[piinfos[0]->pinfo.getFineNbrInfo(Side<2>::east()).local_indexes[1]]);
  REQUIRE(iface.patches.size() == 3);
  CHECK(iface.patches[2].side == Side<2>::west());
  CHECK(iface.patches[2].type.isFineToCoarse());
  CHECK(iface.patches[2].type.getOrthant() == Orthant<1>::upper());
  CHECK(iface.patches[0].piinfo == piinfos[0]);
  CHECK(iface.patches[1].piinfo == piinfos[piinfos[0]->pinfo.getFineNbrInfo(Side<2>::east()).local_indexes[0]]);
  CHECK(iface.patches[2].piinfo == piinfos[piinfos[0]->pinfo.getFineNbrInfo(Side<2>::east()).local_indexes[1]]);
}
TEST_CASE("Schur::Interface insert Fine interface")
{
  DomainReader<2> domain_reader("mesh_inputs/2d_refined_east_1x2_mpi1.json", { 10, 10 }, 0);
  auto domain = domain_reader.getFinerDomain();
  vector<shared_ptr<Schur::PatchIfaceInfo<2>>> piinfos;
  for (auto patch : domain.getPatchInfoVector()) {
    piinfos.push_back(make_shared<Schur::PatchIfaceInfo<2>>(patch));
  }
  int id = piinfos[0]->getFineIfaceInfo(Side<2>::east())->fine_ids[0];
  Schur::Interface<2> iface(id);

  iface.insert(Side<2>::east(), piinfos[0]);
  REQUIRE(iface.patches.size() == 1);
  CHECK(iface.patches[0].side == Side<2>::east());
  CHECK(iface.patches[0].type.isCoarseToFine());
  CHECK(iface.patches[0].type.getOrthant() == Orthant<1>::lower());
  CHECK(iface.patches[0].piinfo == piinfos[0]);

  iface.insert(Side<2>::west(), piinfos[piinfos[0]->pinfo.getFineNbrInfo(Side<2>::east()).local_indexes[0]]);
  REQUIRE(iface.patches.size() == 2);
  CHECK(iface.patches[1].side == Side<2>::west());
  CHECK(iface.patches[1].type.isFineToFine());
  CHECK(iface.patches[1].type.getOrthant() == Orthant<1>::lower());
  CHECK(iface.patches[0].piinfo == piinfos[0]);
  CHECK(iface.patches[1].piinfo == piinfos[piinfos[0]->pinfo.getFineNbrInfo(Side<2>::east()).local_indexes[0]]);
}
TEST_CASE("Schur::Interface merge Fine interface")
{
  DomainReader<2> domain_reader("mesh_inputs/2d_refined_east_1x2_mpi1.json", { 10, 10 }, 0);
  auto domain = domain_reader.getFinerDomain();
  vector<shared_ptr<Schur::PatchIfaceInfo<2>>> piinfos;
  for (auto patch : domain.getPatchInfoVector()) {
    piinfos.push_back(make_shared<Schur::PatchIfaceInfo<2>>(patch));
  }
  int id = piinfos[0]->getFineIfaceInfo(Side<2>::east())->fine_ids[0];

  Schur::Interface<2> iface(id);
  iface.insert(Side<2>::east(), piinfos[0]);

  Schur::Interface<2> iface2(id);
  iface2.insert(Side<2>::west(), piinfos[piinfos[0]->pinfo.getFineNbrInfo(Side<2>::east()).local_indexes[0]]);

  iface.merge(iface2);
  REQUIRE(iface.patches.size() == 2);
  CHECK(iface.patches[0].side == Side<2>::east());
  CHECK(iface.patches[0].type.isCoarseToFine());
  CHECK(iface.patches[0].type.getOrthant() == Orthant<1>::lower());
  CHECK(iface.patches[0].piinfo == piinfos[0]);
  CHECK(iface.patches[1].side == Side<2>::west());
  CHECK(iface.patches[1].type.isFineToFine());
  CHECK(iface.patches[1].type.getOrthant() == Orthant<1>::lower());
  CHECK(iface.patches[1].piinfo == piinfos[piinfos[0]->pinfo.getFineNbrInfo(Side<2>::east()).local_indexes[0]]);
}
TEST_CASE("Schur::Interface enumerateIfacesFromPiinfoVector")
{
  DomainReader<2> domain_reader("mesh_inputs/2d_refined_east_1x2_mpi1.json", { 10, 10 }, 0);
  auto domain = domain_reader.getFinerDomain();
  vector<shared_ptr<const Schur::PatchIfaceInfo<2>>> piinfos;
  for (auto patch : domain.getPatchInfoVector()) {
    piinfos.push_back(make_shared<Schur::PatchIfaceInfo<2>>(patch));
  }

  map<int, map<int, std::shared_ptr<Schur::Interface<2>>>> ifaces;
  vector<std::shared_ptr<Schur::PatchIfaceInfo<2>>> off_proc_piinfos;
  Schur::Interface<2>::EnumerateIfacesFromPiinfoVector(piinfos, ifaces, off_proc_piinfos);
  CHECK(ifaces.size() == 1);
  CHECK(ifaces[0].size() == 7);
  CHECK(off_proc_piinfos.size() == 0);

  auto coarse_piinfo = piinfos[0];
  auto ref_sw_piinfo = piinfos[coarse_piinfo->pinfo.getFineNbrInfo(Side<2>::east()).local_indexes[0]];
  auto ref_se_piinfo = piinfos[ref_sw_piinfo->pinfo.getNormalNbrInfo(Side<2>::east()).local_index];
  auto ref_nw_piinfo = piinfos[coarse_piinfo->pinfo.getFineNbrInfo(Side<2>::east()).local_indexes[1]];
  auto ref_ne_piinfo = piinfos[ref_nw_piinfo->pinfo.getNormalNbrInfo(Side<2>::east()).local_index];

  // check coarse interface
  {
    int id = coarse_piinfo->getFineIfaceInfo(Side<2>::east())->id;
    auto iface = ifaces[0].at(id);
    CHECK(iface->id == id);
    CHECK(iface->patches.size() == 3);

    Schur::Interface<2>::SideTypePiinfo* coarse_patch = nullptr;
    for (auto& patch : iface->patches) {
      if (patch.piinfo->pinfo.id == coarse_piinfo->pinfo.id) {
        coarse_patch = &patch;
        break;
      }
    }

    REQUIRE(coarse_patch != nullptr);
    CHECK(coarse_patch->type.isCoarseToCoarse());
    CHECK(coarse_patch->side == Side<2>::east());

    Schur::Interface<2>::SideTypePiinfo* ref_sw_patch = nullptr;
    for (auto& patch : iface->patches) {
      if (patch.piinfo->pinfo.id == ref_sw_piinfo->pinfo.id) {
        ref_sw_patch = &patch;
        break;
      }
    }

    REQUIRE(ref_sw_patch != nullptr);
    CHECK(ref_sw_patch->type.isFineToCoarse());
    CHECK(ref_sw_patch->side == Side<2>::west());

    Schur::Interface<2>::SideTypePiinfo* ref_nw_patch = nullptr;
    for (auto& patch : iface->patches) {
      if (patch.piinfo->pinfo.id == ref_nw_piinfo->pinfo.id) {
        ref_nw_patch = &patch;
        break;
      }
    }

    REQUIRE(ref_nw_patch != nullptr);
    CHECK(ref_nw_patch->type.isFineToCoarse());
    CHECK(ref_nw_patch->side == Side<2>::west());
  }
  // check sw fine interface
  {
    int id = coarse_piinfo->getFineIfaceInfo(Side<2>::east())->fine_ids[0];
    auto iface = ifaces[0].at(id);
    CHECK(iface->id == id);
    CHECK(iface->patches.size() == 2);

    Schur::Interface<2>::SideTypePiinfo* coarse_patch = nullptr;
    for (auto& patch : iface->patches) {
      if (patch.piinfo->pinfo.id == coarse_piinfo->pinfo.id) {
        coarse_patch = &patch;
        break;
      }
    }

    REQUIRE(coarse_patch != nullptr);
    CHECK(coarse_patch->type.isCoarseToFine());
    CHECK(coarse_patch->side == Side<2>::east());

    Schur::Interface<2>::SideTypePiinfo* ref_sw_patch = nullptr;
    for (auto& patch : iface->patches) {
      if (patch.piinfo->pinfo.id == ref_sw_piinfo->pinfo.id) {
        ref_sw_patch = &patch;
        break;
      }
    }

    REQUIRE(ref_sw_patch != nullptr);
    CHECK(ref_sw_patch->type.isFineToFine());
    CHECK(ref_sw_patch->side == Side<2>::west());
  }
  // check nw fine interface
  {
    int id = coarse_piinfo->getFineIfaceInfo(Side<2>::east())->fine_ids[1];
    auto iface = ifaces[0].at(id);
    CHECK(iface->id == id);
    CHECK(iface->patches.size() == 2);

    Schur::Interface<2>::SideTypePiinfo* coarse_patch = nullptr;
    for (auto& patch : iface->patches) {
      if (patch.piinfo->pinfo.id == coarse_piinfo->pinfo.id) {
        coarse_patch = &patch;
        break;
      }
    }

    REQUIRE(coarse_patch != nullptr);
    CHECK(coarse_patch->type.isCoarseToFine());
    CHECK(coarse_patch->side == Side<2>::east());

    Schur::Interface<2>::SideTypePiinfo* ref_nw_patch = nullptr;
    for (auto& patch : iface->patches) {
      if (patch.piinfo->pinfo.id == ref_nw_piinfo->pinfo.id) {
        ref_nw_patch = &patch;
        break;
      }
    }

    REQUIRE(ref_nw_patch != nullptr);
    CHECK(ref_nw_patch->type.isFineToFine());
    CHECK(ref_nw_patch->side == Side<2>::west());
  }
  // check interface between se and ne
  {
    int id = ref_se_piinfo->getNormalIfaceInfo(Side<2>::north())->id;
    auto iface = ifaces[0].at(id);
    CHECK(iface->id == id);
    CHECK(iface->patches.size() == 2);

    Schur::Interface<2>::SideTypePiinfo* ref_se_patch = nullptr;
    for (auto& patch : iface->patches) {
      if (patch.piinfo->pinfo.id == ref_se_piinfo->pinfo.id) {
        ref_se_patch = &patch;
        break;
      }
    }

    REQUIRE(ref_se_patch != nullptr);
    CHECK(ref_se_patch->type.isNormal());
    CHECK(ref_se_patch->side == Side<2>::north());

    Schur::Interface<2>::SideTypePiinfo* ref_ne_patch = nullptr;
    for (auto& patch : iface->patches) {
      if (patch.piinfo->pinfo.id == ref_ne_piinfo->pinfo.id) {
        ref_ne_patch = &patch;
        break;
      }
    }

    REQUIRE(ref_ne_patch != nullptr);
    CHECK(ref_ne_patch->type.isNormal());
    CHECK(ref_ne_patch->side == Side<2>::south());
  }
}
TEST_CASE("Schur::Interface serialize Normal interface")
{
  DomainReader<2> domain_reader("mesh_inputs/2d_uniform_1x2_mpi1.json", { 10, 10 }, 0);
  auto domain = domain_reader.getFinerDomain();
  vector<shared_ptr<Schur::PatchIfaceInfo<2>>> piinfos;
  for (auto patch : domain.getPatchInfoVector()) {
    piinfos.push_back(make_shared<Schur::PatchIfaceInfo<2>>(patch));
  }
  int id = piinfos[0]->getNormalIfaceInfo(Side<2>::east())->id;
  Schur::Interface<2> iface_in(id);
  iface_in.insert(Side<2>::east(), piinfos[0]);
  iface_in.insert(Side<2>::west(), piinfos[1]);

  char* buff = new char[iface_in.serialize(nullptr)];
  iface_in.serialize(buff);
  Schur::Interface<2> iface;
  iface.deserialize(buff);
  delete[] buff;

  REQUIRE(iface.patches.size() == 2);

  CHECK(iface.patches[0].side == Side<2>::east());
  CHECK(iface.patches[0].type.isNormal());
  CHECK(iface.patches[0].piinfo->pinfo.id == piinfos[0]->pinfo.id);

  CHECK(iface.patches[1].side == Side<2>::west());
  CHECK(iface.patches[1].type.isNormal());
  CHECK(iface.patches[1].piinfo->pinfo.id == piinfos[1]->pinfo.id);
}
TEST_CASE("Schur::Interface serialize Coarse interface")
{
  DomainReader<2> domain_reader("mesh_inputs/2d_refined_east_1x2_mpi1.json", { 10, 10 }, 0);
  auto domain = domain_reader.getFinerDomain();
  vector<shared_ptr<Schur::PatchIfaceInfo<2>>> piinfos;
  for (auto patch : domain.getPatchInfoVector()) {
    piinfos.push_back(make_shared<Schur::PatchIfaceInfo<2>>(patch));
  }
  int id = piinfos[0]->getFineIfaceInfo(Side<2>::east())->id;
  Schur::Interface<2> iface_in(id);

  iface_in.insert(Side<2>::east(), piinfos[0]);
  iface_in.insert(Side<2>::west(), piinfos[piinfos[0]->pinfo.getFineNbrInfo(Side<2>::east()).local_indexes[0]]);
  iface_in.insert(Side<2>::west(), piinfos[piinfos[0]->pinfo.getFineNbrInfo(Side<2>::east()).local_indexes[1]]);

  char* buff = new char[iface_in.serialize(nullptr)];
  iface_in.serialize(buff);
  Schur::Interface<2> iface;
  iface.deserialize(buff);
  delete[] buff;

  REQUIRE(iface.patches.size() == 3);
  CHECK(iface.patches[0].side == Side<2>::east());
  CHECK(iface.patches[0].type.isCoarseToCoarse());
  CHECK(iface.patches[0].piinfo->pinfo.id == piinfos[0]->pinfo.id);

  CHECK(iface.patches[1].side == Side<2>::west());
  CHECK(iface.patches[1].type.isFineToCoarse());
  CHECK(iface.patches[1].type.getOrthant() == Orthant<1>::lower());
  CHECK(iface.patches[1].piinfo->pinfo.id == piinfos[0]->pinfo.getFineNbrInfo(Side<2>::east()).ids[0]);

  CHECK(iface.patches[2].side == Side<2>::west());
  CHECK(iface.patches[2].type.isFineToCoarse());
  CHECK(iface.patches[2].type.getOrthant() == Orthant<1>::upper());
  CHECK(iface.patches[2].piinfo->pinfo.id == piinfos[0]->pinfo.getFineNbrInfo(Side<2>::east()).ids[1]);
}
TEST_CASE("Schur::Interface serialize Fine interface")
{
  DomainReader<2> domain_reader("mesh_inputs/2d_refined_east_1x2_mpi1.json", { 10, 10 }, 0);
  auto domain = domain_reader.getFinerDomain();
  vector<shared_ptr<Schur::PatchIfaceInfo<2>>> piinfos;
  for (auto patch : domain.getPatchInfoVector()) {
    piinfos.push_back(make_shared<Schur::PatchIfaceInfo<2>>(patch));
  }
  int id = piinfos[0]->getFineIfaceInfo(Side<2>::east())->fine_ids[0];
  Schur::Interface<2> iface_in(id);

  iface_in.insert(Side<2>::east(), piinfos[0]);
  iface_in.insert(Side<2>::west(), piinfos[piinfos[0]->pinfo.getFineNbrInfo(Side<2>::east()).local_indexes[0]]);

  char* buff = new char[iface_in.serialize(nullptr)];
  iface_in.serialize(buff);
  Schur::Interface<2> iface;
  iface.deserialize(buff);
  delete[] buff;

  REQUIRE(iface.patches.size() == 2);

  CHECK(iface.patches[0].side == Side<2>::east());
  CHECK(iface.patches[0].type.isCoarseToFine());
  CHECK(iface.patches[0].type.getOrthant() == Orthant<1>::lower());
  CHECK(iface.patches[0].piinfo->pinfo.id == piinfos[0]->pinfo.id);

  CHECK(iface.patches[1].side == Side<2>::west());
  CHECK(iface.patches[1].type.isFineToFine());
  CHECK(iface.patches[1].type.getOrthant() == Orthant<1>::lower());
  CHECK(iface.patches[1].piinfo->pinfo.id == piinfos[0]->pinfo.getFineNbrInfo(Side<2>::east()).ids[0]);
}
