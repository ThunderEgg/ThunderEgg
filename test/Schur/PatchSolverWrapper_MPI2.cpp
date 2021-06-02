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
#include "PatchSolverWrapper_MOCKS.h"
#include <ThunderEgg/Schur/PatchSolverWrapper.h>
#include <ThunderEgg/Schur/ValVectorGenerator.h>
#include <ThunderEgg/ValVectorGenerator.h>

#include <limits>

#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators.hpp>

using namespace std;
using namespace ThunderEgg;

#define MESHES                                             \
	"mesh_inputs/2d_refined_east_1x2_east_on_1_mpi2.json", \
	"mesh_inputs/2d_uniform_1x2_east_on_1_mpi2.json"

TEST_CASE("Schur::PatchSolverWrapper<2> throws exception for non-square patches",
          "[Schur::PatchSolverWrapper]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH: " << mesh_file);
	auto nx = GENERATE(5, 7);
	INFO("NX: " << nx);
	auto ny = GENERATE(6, 8);
	INFO("NY: " << ny);
	DomainReader<2> domain_reader(mesh_file, {nx, ny}, 1);
	auto            domain       = domain_reader.getFinerDomain();
	auto            iface_domain = make_shared<Schur::InterfaceDomain<2>>(domain);
	auto            ghost_filler = make_shared<MockGhostFiller<2>>();
	auto            solver       = make_shared<MockPatchSolver<2>>(domain, ghost_filler);

	CHECK_THROWS_AS(Schur::PatchSolverWrapper<2>(iface_domain, solver), RuntimeError);
}
TEST_CASE("Schur::PatchSolverWrapper<2> apply fills ghost in rhs as expected",
          "[Schur::PatchSolverWrapper]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH: " << mesh_file);
	auto n = GENERATE(5, 7);
	INFO("N: " << n);
	auto            schur_fill_value = GENERATE(1, 1.3, 8, 2, -1);
	DomainReader<2> domain_reader(mesh_file, {n, n}, 1);
	auto            domain       = domain_reader.getFinerDomain();
	auto            iface_domain = make_shared<Schur::InterfaceDomain<2>>(domain);
	auto            ghost_filler = make_shared<MockGhostFiller<2>>();
	auto            solver
	= make_shared<RHSGhostCheckingPatchSolver<2>>(domain, ghost_filler, schur_fill_value);

	Schur::ValVectorGenerator<1> vg(iface_domain);

	auto x = vg.getNewVector();
	auto b = vg.getNewVector();

	x->set(schur_fill_value);

	// checking will be done in the solver
	Schur::PatchSolverWrapper<2> psw(iface_domain, solver);
	psw.apply(x, b);
	CHECK(solver->wasCalled());
}
TEST_CASE("Schur::PatchSolverWrapper<2> apply gives expected rhs value for Schur matrix",
          "[Schur::PatchSolverWrapper]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH: " << mesh_file);
	auto n = GENERATE(5, 7);
	INFO("N: " << n);
	auto            schur_fill_value  = GENERATE(1, 1.3, 8, 2, -1);
	auto            domain_fill_value = GENERATE(1, 1.3, 8, 2, -1);
	DomainReader<2> domain_reader(mesh_file, {n, n}, 1);
	auto            domain       = domain_reader.getFinerDomain();
	auto            iface_domain = make_shared<Schur::InterfaceDomain<2>>(domain);
	auto            ghost_filler = make_shared<PatchFillingGhostFiller<2>>(domain_fill_value);
	auto            solver       = make_shared<MockPatchSolver<2>>(domain, ghost_filler);

	Schur::ValVectorGenerator<1> vg(iface_domain);

	auto x = vg.getNewVector();
	auto b = vg.getNewVector();

	x->set(schur_fill_value);

	Schur::PatchSolverWrapper<2> psw(iface_domain, solver);
	psw.apply(x, b);
	CHECK(solver->allPatchesCalled());
	CHECK(ghost_filler->wasCalled());
	for (int i = 0; i < b->getNumLocalPatches(); i++) {
		auto local_data = b->getView(0, i);
		nested_loop<1>(local_data.getStart(), local_data.getEnd(),
		               [&](const std::array<int, 1> &coord) {
			               CHECK(local_data[coord] == Catch::Approx(schur_fill_value - domain_fill_value));
		               });
	}
}
TEST_CASE(
"Schur::PatchSolverWrapper<2> apply gives expected rhs value for Schur matrix with rhs already set",
"[Schur::PatchSolverWrapper]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH: " << mesh_file);
	auto n = GENERATE(5, 7);
	INFO("N: " << n);
	auto            schur_fill_value  = GENERATE(1, 1.3, 8, 2, -1);
	auto            domain_fill_value = GENERATE(1, 1.3, 8, 2, -1);
	DomainReader<2> domain_reader(mesh_file, {n, n}, 1);
	auto            domain       = domain_reader.getFinerDomain();
	auto            iface_domain = make_shared<Schur::InterfaceDomain<2>>(domain);
	auto            ghost_filler = make_shared<PatchFillingGhostFiller<2>>(domain_fill_value);
	auto            solver       = make_shared<MockPatchSolver<2>>(domain, ghost_filler);

	Schur::ValVectorGenerator<1> vg(iface_domain);

	auto x = vg.getNewVector();
	auto b = vg.getNewVector();

	x->set(schur_fill_value);
	b->set(99);

	Schur::PatchSolverWrapper<2> psw(iface_domain, solver);
	psw.apply(x, b);
	CHECK(solver->allPatchesCalled());
	CHECK(ghost_filler->wasCalled());
	for (int i = 0; i < b->getNumLocalPatches(); i++) {
		auto local_data = b->getView(0, i);
		nested_loop<1>(local_data.getStart(), local_data.getEnd(),
		               [&](const std::array<int, 1> &coord) {
			               CHECK(local_data[coord] == Catch::Approx(schur_fill_value - domain_fill_value));
		               });
	}
}
TEST_CASE("Schur::PatchSolverWrapper<2> getSchurRHSFromDomainRHS fills ghost in rhs as expected",
          "[Schur::PatchSolverWrapper]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH: " << mesh_file);
	auto n = GENERATE(5, 7);
	INFO("N: " << n);
	auto            domain_fill_value = GENERATE(1, 1.3, 8, 2, -1);
	DomainReader<2> domain_reader(mesh_file, {n, n}, 1);
	auto            domain       = domain_reader.getFinerDomain();
	auto            iface_domain = make_shared<Schur::InterfaceDomain<2>>(domain);
	auto            ghost_filler = make_shared<MockGhostFiller<2>>();
	auto            solver       = make_shared<RHSGhostCheckingPatchSolver<2>>(domain, ghost_filler, 0);

	Schur::ValVectorGenerator<1> vg(iface_domain);
	ValVectorGenerator<2>        domain_vg(domain, 1);

	auto schur_b  = vg.getNewVector();
	auto domain_b = domain_vg.getNewVector();

	domain_b->set(domain_fill_value);

	// checking will be done in the solver
	Schur::PatchSolverWrapper<2> psw(iface_domain, solver);
	psw.getSchurRHSFromDomainRHS(domain_b, schur_b);
	CHECK(solver->wasCalled());
}
TEST_CASE(
"Schur::PatchSolverWrapper<2> getSchurRHSFromDomainRHS gives expected rhs value for Schur matrix",
"[Schur::PatchSolverWrapper]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH: " << mesh_file);
	auto n = GENERATE(5, 7);
	INFO("N: " << n);
	auto            domain_fill_value = GENERATE(1, 1.3, 8, 2, -1);
	DomainReader<2> domain_reader(mesh_file, {n, n}, 1);
	auto            domain       = domain_reader.getFinerDomain();
	auto            iface_domain = make_shared<Schur::InterfaceDomain<2>>(domain);
	auto            ghost_filler = make_shared<PatchFillingGhostFiller<2>>(domain_fill_value);
	auto            solver       = make_shared<MockPatchSolver<2>>(domain, ghost_filler);

	Schur::ValVectorGenerator<1> vg(iface_domain);
	ValVectorGenerator<2>        domain_vg(domain, 1);

	auto schur_b  = vg.getNewVector();
	auto domain_b = domain_vg.getNewVector();

	Schur::PatchSolverWrapper<2> psw(iface_domain, solver);
	psw.getSchurRHSFromDomainRHS(domain_b, schur_b);
	CHECK(solver->allPatchesCalled());
	CHECK(ghost_filler->wasCalled());
	for (int i = 0; i < schur_b->getNumLocalPatches(); i++) {
		auto local_data = schur_b->getView(0, i);
		nested_loop<1>(local_data.getStart(), local_data.getEnd(),
		               [&](const std::array<int, 1> &coord) {
			               CHECK(local_data[coord] == Catch::Approx(domain_fill_value));
		               });
	}
}
TEST_CASE(
"Schur::PatchSolverWrapper<2> getSchurRHSFromDomainRHS gives expected rhs value for Schur matrix with rhs already set",
"[Schur::PatchSolverWrapper]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH: " << mesh_file);
	auto n = GENERATE(5, 7);
	INFO("N: " << n);
	auto            domain_fill_value = GENERATE(1, 1.3, 8, 2, -1);
	DomainReader<2> domain_reader(mesh_file, {n, n}, 1);
	auto            domain       = domain_reader.getFinerDomain();
	auto            iface_domain = make_shared<Schur::InterfaceDomain<2>>(domain);
	auto            ghost_filler = make_shared<PatchFillingGhostFiller<2>>(domain_fill_value);
	auto            solver       = make_shared<MockPatchSolver<2>>(domain, ghost_filler);

	Schur::ValVectorGenerator<1> vg(iface_domain);
	ValVectorGenerator<2>        domain_vg(domain, 1);

	auto schur_b  = vg.getNewVector();
	auto domain_b = domain_vg.getNewVector();

	schur_b->set(99);

	Schur::PatchSolverWrapper<2> psw(iface_domain, solver);
	psw.getSchurRHSFromDomainRHS(domain_b, schur_b);
	CHECK(solver->allPatchesCalled());
	CHECK(ghost_filler->wasCalled());
	for (int i = 0; i < schur_b->getNumLocalPatches(); i++) {
		auto local_data = schur_b->getView(0, i);
		nested_loop<1>(local_data.getStart(), local_data.getEnd(),
		               [&](const std::array<int, 1> &coord) {
			               CHECK(local_data[coord] == Catch::Approx(domain_fill_value));
		               });
	}
}