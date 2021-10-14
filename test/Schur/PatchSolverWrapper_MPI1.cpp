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
#include "PatchSolverWrapper_MOCKS.h"
#include <ThunderEgg/Schur/PatchSolverWrapper.h>

#include <limits>

#include <catch2/generators/catch_generators.hpp>

using namespace std;
using namespace ThunderEgg;
#define MESHES "mesh_inputs/2d_refined_east_1x2_mpi1.json", "mesh_inputs/2d_uniform_1x2_mpi1.json"
TEST_CASE("Schur::PatchSolverWrapper<2> throws exception for non-square patches",
          "[Schur::PatchSolverWrapper]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH: " << mesh_file);
	auto nx = GENERATE(5, 7);
	INFO("NX: " << nx);
	auto ny = GENERATE(6, 8);
	INFO("NY: " << ny);
	DomainReader<2>           domain_reader(mesh_file, {nx, ny}, 1);
	auto                      domain = domain_reader.getFinerDomain();
	Schur::InterfaceDomain<2> iface_domain(domain);
	MockGhostFiller<2>        ghost_filler;
	MockPatchSolver<2>        solver(domain, ghost_filler);

	CHECK_THROWS_AS(Schur::PatchSolverWrapper<2>(iface_domain, solver), RuntimeError);
}
TEST_CASE("Schur::PatchSolverWrapper<2> apply fills ghost in rhs as expected",
          "[Schur::PatchSolverWrapper]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH: " << mesh_file);
	auto n = GENERATE(5, 7);
	INFO("N: " << n);
	auto                           schur_fill_value = GENERATE(1, 1.3, 8, 2, -1);
	DomainReader<2>                domain_reader(mesh_file, {n, n}, 1);
	auto                           domain = domain_reader.getFinerDomain();
	Schur::InterfaceDomain<2>      iface_domain(domain);
	MockGhostFiller<2>             ghost_filler;
	RHSGhostCheckingPatchSolver<2> solver(domain, ghost_filler, schur_fill_value);

	Vector<1> x = iface_domain.getNewVector();
	Vector<1> b = iface_domain.getNewVector();

	x.set(schur_fill_value);

	// checking will be done in the solver
	Schur::PatchSolverWrapper<2> psw(iface_domain, solver);
	psw.apply(x, b);
	CHECK(solver.wasCalled());
}
TEST_CASE("Schur::PatchSolverWrapper<2> apply gives expected rhs value for Schur matrix",
          "[Schur::PatchSolverWrapper]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH: " << mesh_file);
	auto n = GENERATE(5, 7);
	INFO("N: " << n);
	auto                       schur_fill_value  = GENERATE(1, 1.3, 8, 2, -1);
	auto                       domain_fill_value = GENERATE(1, 1.3, 8, 2, -1);
	DomainReader<2>            domain_reader(mesh_file, {n, n}, 1);
	auto                       domain = domain_reader.getFinerDomain();
	Schur::InterfaceDomain<2>  iface_domain(domain);
	PatchFillingGhostFiller<2> ghost_filler(domain_fill_value);
	MockPatchSolver<2>         solver(domain, ghost_filler);

	Vector<1> x = iface_domain.getNewVector();
	Vector<1> b = iface_domain.getNewVector();

	x.set(schur_fill_value);

	Schur::PatchSolverWrapper<2> psw(iface_domain, solver);
	psw.apply(x, b);
	CHECK(solver.allPatchesCalled());
	CHECK(ghost_filler.wasCalled());
	for (int i = 0; i < b.getNumLocalPatches(); i++) {
		auto local_data = b.getComponentView(0, i);
		Loop::Nested<1>(local_data.getStart(), local_data.getEnd(),
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
	auto                       schur_fill_value  = GENERATE(1, 1.3, 8, 2, -1);
	auto                       domain_fill_value = GENERATE(1, 1.3, 8, 2, -1);
	DomainReader<2>            domain_reader(mesh_file, {n, n}, 1);
	auto                       domain = domain_reader.getFinerDomain();
	Schur::InterfaceDomain<2>  iface_domain(domain);
	PatchFillingGhostFiller<2> ghost_filler(domain_fill_value);
	MockPatchSolver<2>         solver(domain, ghost_filler);

	Vector<1> x = iface_domain.getNewVector();
	Vector<1> b = iface_domain.getNewVector();

	x.set(schur_fill_value);
	b.set(99);

	Schur::PatchSolverWrapper<2> psw(iface_domain, solver);
	psw.apply(x, b);
	CHECK(solver.allPatchesCalled());
	CHECK(ghost_filler.wasCalled());
	for (int i = 0; i < b.getNumLocalPatches(); i++) {
		auto local_data = b.getComponentView(0, i);
		Loop::Nested<1>(local_data.getStart(), local_data.getEnd(),
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
	auto                           domain_fill_value = GENERATE(1, 1.3, 8, 2, -1);
	DomainReader<2>                domain_reader(mesh_file, {n, n}, 1);
	auto                           domain = domain_reader.getFinerDomain();
	Schur::InterfaceDomain<2>      iface_domain(domain);
	MockGhostFiller<2>             ghost_filler;
	RHSGhostCheckingPatchSolver<2> solver(domain, ghost_filler, 0);

	Vector<1> schur_b = iface_domain.getNewVector();
	Vector<2> domain_b(domain, 1);

	domain_b.set(domain_fill_value);

	// checking will be done in the solver
	Schur::PatchSolverWrapper<2> psw(iface_domain, solver);
	psw.getSchurRHSFromDomainRHS(domain_b, schur_b);
	CHECK(solver.wasCalled());
}
TEST_CASE(
"Schur::PatchSolverWrapper<2> getSchurRHSFromDomainRHS gives expected rhs value for Schur matrix",
"[Schur::PatchSolverWrapper]")
{
	auto mesh_file = GENERATE(as<std::string>{}, MESHES);
	INFO("MESH: " << mesh_file);
	auto n = GENERATE(5, 7);
	INFO("N: " << n);
	auto                       domain_fill_value = GENERATE(1, 1.3, 8, 2, -1);
	DomainReader<2>            domain_reader(mesh_file, {n, n}, 1);
	auto                       domain = domain_reader.getFinerDomain();
	Schur::InterfaceDomain<2>  iface_domain(domain);
	PatchFillingGhostFiller<2> ghost_filler(domain_fill_value);
	MockPatchSolver<2>         solver(domain, ghost_filler);

	Vector<1> schur_b = iface_domain.getNewVector();
	Vector<2> domain_b(domain, 1);

	Schur::PatchSolverWrapper<2> psw(iface_domain, solver);
	psw.getSchurRHSFromDomainRHS(domain_b, schur_b);
	CHECK(solver.allPatchesCalled());
	CHECK(ghost_filler.wasCalled());
	for (int i = 0; i < schur_b.getNumLocalPatches(); i++) {
		auto local_data = schur_b.getComponentView(0, i);
		Loop::Nested<1>(local_data.getStart(), local_data.getEnd(),
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
	auto                       domain_fill_value = GENERATE(1, 1.3, 8, 2, -1);
	DomainReader<2>            domain_reader(mesh_file, {n, n}, 1);
	auto                       domain = domain_reader.getFinerDomain();
	Schur::InterfaceDomain<2>  iface_domain(domain);
	PatchFillingGhostFiller<2> ghost_filler(domain_fill_value);
	MockPatchSolver<2>         solver(domain, ghost_filler);

	Vector<1> schur_b = iface_domain.getNewVector();
	Vector<2> domain_b(domain, 1);

	schur_b.set(99);

	Schur::PatchSolverWrapper<2> psw(iface_domain, solver);
	psw.getSchurRHSFromDomainRHS(domain_b, schur_b);
	CHECK(solver.allPatchesCalled());
	CHECK(ghost_filler.wasCalled());
	for (int i = 0; i < schur_b.getNumLocalPatches(); i++) {
		auto local_data = schur_b.getComponentView(0, i);
		Loop::Nested<1>(local_data.getStart(), local_data.getEnd(),
		                [&](const std::array<int, 1> &coord) {
			                CHECK(local_data[coord] == Catch::Approx(domain_fill_value));
		                });
	}
}