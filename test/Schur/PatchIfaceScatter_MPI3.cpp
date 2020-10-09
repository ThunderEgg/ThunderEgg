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

#include "../utils/DomainReader.h"
#include "catch.hpp"
#include <ThunderEgg/Schur/PatchIfaceScatter.h>
#include <ThunderEgg/Schur/ValVectorGenerator.h>
#include <limits>
using namespace std;
using namespace ThunderEgg;
#define MESHES "mesh_inputs/2d_refined_complicated_mpi3.json"
/*****************
 *
 *
 *  Exception cases
 *
 *
 *****************/
TEST_CASE("Schur::PatchIfaceScatter<2> throws exception for non-square patches",
          "[Schur::PatchIfaceScatter]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH: " << mesh_file);
	auto nx = GENERATE(5, 7);
	INFO("NX: " << nx);
	auto ny = GENERATE(6, 8);
	INFO("NY: " << ny);
	DomainReader<2> domain_reader(mesh_file, {nx, ny}, 0);
	auto            domain       = domain_reader.getFinerDomain();
	auto            iface_domain = make_shared<Schur::InterfaceDomain<2>>(domain);

	CHECK_THROWS_AS(Schur::PatchIfaceScatter<2>(iface_domain), RuntimeError);
}
TEST_CASE("Schur::PatchIfaceScatter<2> scatterStart throws exception when called twice",
          "[Schur::PatchIfaceScatter]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH: " << mesh_file);
	auto n = GENERATE(5, 10);
	INFO("N" << n);

	DomainReader<2> domain_reader(mesh_file, {n, n}, 0);
	auto            domain       = domain_reader.getFinerDomain();
	auto            iface_domain = make_shared<Schur::InterfaceDomain<2>>(domain);

	Schur::PatchIfaceScatter<2>  scatter(iface_domain);
	Schur::ValVectorGenerator<1> vg(iface_domain);

	auto global_vector = vg.getNewVector();
	auto local_vector  = scatter.getNewLocalPatchIfaceVector();

	scatter.scatterStart(global_vector, local_vector);
	CHECK_THROWS_AS(scatter.scatterStart(global_vector, local_vector), RuntimeError);
}
TEST_CASE(
"Schur::PatchIfaceScatter<2> scatterFinish throws exception when called with different global vectors",
"[Schur::PatchIfaceScatter]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH: " << mesh_file);
	auto n = GENERATE(5, 10);
	INFO("N" << n);

	DomainReader<2> domain_reader(mesh_file, {n, n}, 0);
	auto            domain       = domain_reader.getFinerDomain();
	auto            iface_domain = make_shared<Schur::InterfaceDomain<2>>(domain);

	Schur::PatchIfaceScatter<2>  scatter(iface_domain);
	Schur::ValVectorGenerator<1> vg(iface_domain);

	auto global_vector   = vg.getNewVector();
	auto global_vector_2 = vg.getNewVector();
	auto local_vector    = scatter.getNewLocalPatchIfaceVector();

	scatter.scatterStart(global_vector, local_vector);
	CHECK_THROWS_AS(scatter.scatterFinish(global_vector_2, local_vector), RuntimeError);
}
TEST_CASE(
"Schur::PatchIfaceScatter<2> scatterFinish throws exception when called with different local vectors",
"[Schur::PatchIfaceScatter]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH: " << mesh_file);
	auto n = GENERATE(5, 10);
	INFO("N" << n);

	DomainReader<2> domain_reader(mesh_file, {n, n}, 0);
	auto            domain       = domain_reader.getFinerDomain();
	auto            iface_domain = make_shared<Schur::InterfaceDomain<2>>(domain);

	Schur::PatchIfaceScatter<2>  scatter(iface_domain);
	Schur::ValVectorGenerator<1> vg(iface_domain);

	auto global_vector  = vg.getNewVector();
	auto local_vector   = scatter.getNewLocalPatchIfaceVector();
	auto local_vector_2 = scatter.getNewLocalPatchIfaceVector();

	scatter.scatterStart(global_vector, local_vector);
	CHECK_THROWS_AS(scatter.scatterFinish(global_vector, local_vector_2), RuntimeError);
}
TEST_CASE(
"Schur::PatchIfaceScatter<2> scatterFinish throws exception when called with different local and global vectors",
"[Schur::PatchIfaceScatter]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH: " << mesh_file);
	auto n = GENERATE(5, 10);
	INFO("N" << n);

	DomainReader<2> domain_reader(mesh_file, {n, n}, 0);
	auto            domain       = domain_reader.getFinerDomain();
	auto            iface_domain = make_shared<Schur::InterfaceDomain<2>>(domain);

	Schur::PatchIfaceScatter<2>  scatter(iface_domain);
	Schur::ValVectorGenerator<1> vg(iface_domain);

	auto global_vector   = vg.getNewVector();
	auto global_vector_2 = vg.getNewVector();
	auto local_vector    = scatter.getNewLocalPatchIfaceVector();
	auto local_vector_2  = scatter.getNewLocalPatchIfaceVector();

	scatter.scatterStart(global_vector, local_vector);
	CHECK_THROWS_AS(scatter.scatterFinish(global_vector_2, local_vector_2), RuntimeError);
}
/******
 *
 *  getNewLocalPatchIfaceVector
 *
 *
 ********/
TEST_CASE(
"Schur::PatchIfaceScatter<2> getNewLocalPatchIfaceVector returns vector of expected length",
"[Schur::PatchIfaceScatter]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH: " << mesh_file);
	auto n = GENERATE(5, 10);
	INFO("N" << n);

	DomainReader<2> domain_reader(mesh_file, {n, n}, 0);
	auto            domain       = domain_reader.getFinerDomain();
	auto            iface_domain = make_shared<Schur::InterfaceDomain<2>>(domain);

	Schur::PatchIfaceScatter<2> scatter(iface_domain);

	auto     local_vector = scatter.getNewLocalPatchIfaceVector();
	set<int> patch_iface_interfaces;

	for (auto piinfo : iface_domain->getPatchIfaceInfos()) {
		for (Side<2> s : Side<2>::getValues()) {
			if (piinfo->pinfo->hasNbr(s)) {
				auto iface_info = piinfo->getIfaceInfo(s);
				patch_iface_interfaces.insert(iface_info->id);
			}
		}
	}

	CHECK(local_vector->getNumLocalPatches() == (int) patch_iface_interfaces.size());
	CHECK(local_vector->getNumLocalCells() == n * patch_iface_interfaces.size());
}
TEST_CASE(
"Schur::PatchIfaceScatter<2> getNewLocalPatchIfaceVector returns vector with local MPI_Comm",
"[Schur::PatchIfaceScatter]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH: " << mesh_file);
	auto n = GENERATE(5, 10);
	INFO("N" << n);

	DomainReader<2> domain_reader(mesh_file, {n, n}, 0);
	auto            domain       = domain_reader.getFinerDomain();
	auto            iface_domain = make_shared<Schur::InterfaceDomain<2>>(domain);

	Schur::PatchIfaceScatter<2> scatter(iface_domain);

	auto local_vector = scatter.getNewLocalPatchIfaceVector();

	CHECK(local_vector->getMPIComm() == MPI_COMM_SELF);
}
/******
 *
 *  scatter
 *
 *
 ********/
TEST_CASE("Schur::PatchIfaceScatter<2> scatter local interfaces are copied",
          "[Schur::PatchIfaceScatter]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH: " << mesh_file);
	auto n = GENERATE(5, 10);
	INFO("N" << n);

	DomainReader<2> domain_reader(mesh_file, {n, n}, 0);
	auto            domain       = domain_reader.getFinerDomain();
	auto            iface_domain = make_shared<Schur::InterfaceDomain<2>>(domain);

	Schur::PatchIfaceScatter<2>  scatter(iface_domain);
	Schur::ValVectorGenerator<1> vg(iface_domain);

	auto global_vector = vg.getNewVector();
	auto local_vector  = scatter.getNewLocalPatchIfaceVector();

	for (int i = 0; i < global_vector->getNumLocalPatches(); i++) {
		auto iface      = iface_domain->getInterfaces()[i];
		auto local_data = global_vector->getLocalData(0, i);
		nested_loop<1>(local_data.getStart(), local_data.getEnd(),
		               [&](const std::array<int, 1> &coord) {
			               local_data[coord] = iface->global_index + 1 + coord[0];
		               });
	}
	scatter.scatterStart(global_vector, local_vector);
	for (auto iface : iface_domain->getInterfaces()) {
		INFO("IFACE_ID: " << iface->id);
		auto local_data = local_vector->getLocalData(0, iface->local_index);
		nested_loop<1>(local_data.getStart(), local_data.getEnd(),
		               [&](const std::array<int, 1> &coord) {
			               CHECK(local_data[coord] == Approx(iface->global_index + 1 + coord[0]));
		               });
	}
}
TEST_CASE("Schur::PatchIfaceScatter<2> scatter", "[Schur::PatchIfaceScatter]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH: " << mesh_file);
	auto n = GENERATE(5, 10);
	INFO("N" << n);

	DomainReader<2> domain_reader(mesh_file, {n, n}, 0);
	auto            domain       = domain_reader.getFinerDomain();
	auto            iface_domain = make_shared<Schur::InterfaceDomain<2>>(domain);

	Schur::PatchIfaceScatter<2>  scatter(iface_domain);
	Schur::ValVectorGenerator<1> vg(iface_domain);

	auto global_vector = vg.getNewVector();
	auto local_vector  = scatter.getNewLocalPatchIfaceVector();

	for (int i = 0; i < global_vector->getNumLocalPatches(); i++) {
		auto iface      = iface_domain->getInterfaces()[i];
		auto local_data = global_vector->getLocalData(0, i);
		nested_loop<1>(local_data.getStart(), local_data.getEnd(),
		               [&](const std::array<int, 1> &coord) {
			               local_data[coord] = iface->global_index + 1 + coord[0];
		               });
	}
	scatter.scatterStart(global_vector, local_vector);
	scatter.scatterFinish(global_vector, local_vector);
	for (auto piinfo : iface_domain->getPatchIfaceInfos()) {
		INFO("PATCH_ID: " << piinfo->pinfo->id);
		for (Side<2> s : Side<2>::getValues()) {
			if (piinfo->pinfo->hasNbr(s)) {
				INFO("Side: " << s);
				auto iface_info = piinfo->getIfaceInfo(s);
				auto local_data = local_vector->getLocalData(0, iface_info->patch_local_index);
				nested_loop<1>(
				local_data.getStart(), local_data.getEnd(), [&](const std::array<int, 1> &coord) {
					CHECK(local_data[coord] == Approx(iface_info->global_index + 1 + coord[0]));
				});
			}
		}
	}
}
TEST_CASE("Schur::PatchIfaceScatter<2> scatter twice", "[Schur::PatchIfaceScatter]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH: " << mesh_file);
	auto n = GENERATE(5, 10);
	INFO("N" << n);

	DomainReader<2> domain_reader(mesh_file, {n, n}, 0);
	auto            domain       = domain_reader.getFinerDomain();
	auto            iface_domain = make_shared<Schur::InterfaceDomain<2>>(domain);

	Schur::PatchIfaceScatter<2>  scatter(iface_domain);
	Schur::ValVectorGenerator<1> vg(iface_domain);

	auto global_vector = vg.getNewVector();
	auto local_vector  = scatter.getNewLocalPatchIfaceVector();

	for (int i = 0; i < global_vector->getNumLocalPatches(); i++) {
		auto iface      = iface_domain->getInterfaces()[i];
		auto local_data = global_vector->getLocalData(0, i);
		nested_loop<1>(local_data.getStart(), local_data.getEnd(),
		               [&](const std::array<int, 1> &coord) {
			               local_data[coord] = iface->global_index + 1 + coord[0];
		               });
	}
	scatter.scatterStart(global_vector, local_vector);
	scatter.scatterFinish(global_vector, local_vector);
	scatter.scatterStart(global_vector, local_vector);
	scatter.scatterFinish(global_vector, local_vector);
	for (auto piinfo : iface_domain->getPatchIfaceInfos()) {
		INFO("PATCH_ID: " << piinfo->pinfo->id);
		for (Side<2> s : Side<2>::getValues()) {
			if (piinfo->pinfo->hasNbr(s)) {
				INFO("Side: " << s);
				auto iface_info = piinfo->getIfaceInfo(s);
				auto local_data = local_vector->getLocalData(0, iface_info->patch_local_index);
				nested_loop<1>(
				local_data.getStart(), local_data.getEnd(), [&](const std::array<int, 1> &coord) {
					CHECK(local_data[coord] == Approx(iface_info->global_index + 1 + coord[0]));
				});
			}
		}
	}
}
TEST_CASE("Schur::PatchIfaceScatter<2> scatter with local vector already filled",
          "[Schur::PatchIfaceScatter]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH: " << mesh_file);
	auto n = GENERATE(5, 10);
	INFO("N" << n);

	DomainReader<2> domain_reader(mesh_file, {n, n}, 0);
	auto            domain       = domain_reader.getFinerDomain();
	auto            iface_domain = make_shared<Schur::InterfaceDomain<2>>(domain);

	Schur::PatchIfaceScatter<2>  scatter(iface_domain);
	Schur::ValVectorGenerator<1> vg(iface_domain);

	auto global_vector = vg.getNewVector();
	auto local_vector  = scatter.getNewLocalPatchIfaceVector();
	local_vector->setWithGhost(99);

	for (int i = 0; i < global_vector->getNumLocalPatches(); i++) {
		auto iface      = iface_domain->getInterfaces()[i];
		auto local_data = global_vector->getLocalData(0, i);
		nested_loop<1>(local_data.getStart(), local_data.getEnd(),
		               [&](const std::array<int, 1> &coord) {
			               local_data[coord] = iface->global_index + 1 + coord[0];
		               });
	}
	scatter.scatterStart(global_vector, local_vector);
	scatter.scatterFinish(global_vector, local_vector);
	for (auto piinfo : iface_domain->getPatchIfaceInfos()) {
		INFO("PATCH_ID: " << piinfo->pinfo->id);
		for (Side<2> s : Side<2>::getValues()) {
			if (piinfo->pinfo->hasNbr(s)) {
				INFO("Side: " << s);
				auto iface_info = piinfo->getIfaceInfo(s);
				auto local_data = local_vector->getLocalData(0, iface_info->patch_local_index);
				nested_loop<1>(
				local_data.getStart(), local_data.getEnd(), [&](const std::array<int, 1> &coord) {
					CHECK(local_data[coord] == Approx(iface_info->global_index + 1 + coord[0]));
				});
			}
		}
	}
}