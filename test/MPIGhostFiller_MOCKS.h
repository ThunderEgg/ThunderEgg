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

#ifndef MPIGHOSTFILLER_MOCKS_H
#define MPIGHOSTFILLER_MOCKS_H

#include <ThunderEgg/DomainTools.h>
#include <ThunderEgg/MPIGhostFiller.h>
#include <ThunderEgg/ValVector.h>

#include <list>
#include <tuple>

#include <catch2/catch_test_macros.hpp>

namespace ThunderEgg
{
template <int D>
class CallMockMPIGhostFiller : public MPIGhostFiller<D>
{
	private:
	int num_components;

	mutable std::list<std::tuple<const PatchInfo<D> *, Side<D>, const NbrType, const Orthant<D - 1>, int, int>> nbr_calls;

	mutable std::list<std::tuple<const PatchInfo<D> *, Edge, const NbrType, const Orthant<1>, int, int>> edge_nbr_calls;

	mutable std::list<std::tuple<const PatchInfo<D> *, Corner<D>, const NbrType, int, int>> corner_nbr_calls;

	mutable std::list<std::tuple<const PatchInfo<D> *, int>> local_calls;

	public:
	void fillGhostCellsForNbrPatch(const PatchInfo<D> &             pinfo,
	                               const std::vector<LocalData<D>> &local_datas,
	                               std::vector<LocalData<D>> &nbr_datas, Side<D> side,
	                               NbrType nbr_type, Orthant<D - 1> orth_on_coarse) const override
	{
		called = true;
		nbr_calls.emplace_back(&pinfo, side, nbr_type, orth_on_coarse, local_datas.size(), nbr_datas.size());
	}

	void fillGhostCellsForEdgeNbrPatch(const PatchInfo<D> &             pinfo,
	                                   const std::vector<LocalData<D>> &local_datas,
	                                   std::vector<LocalData<D>> &nbr_datas, Edge edge,
	                                   NbrType nbr_type, Orthant<1> orth_on_coarse) const override
	{
		called = true;
		edge_nbr_calls.emplace_back(&pinfo, edge, nbr_type, orth_on_coarse, local_datas.size(), nbr_datas.size());
	}

	void fillGhostCellsForCornerNbrPatch(const PatchInfo<D> &             pinfo,
	                                     const std::vector<LocalData<D>> &local_datas,
	                                     std::vector<LocalData<D>> &nbr_datas, Corner<D> corner,
	                                     NbrType nbr_type) const override
	{
		called = true;
		corner_nbr_calls.emplace_back(&pinfo, corner, nbr_type, local_datas.size(), nbr_datas.size());
	}

	void fillGhostCellsForLocalPatch(const PatchInfo<D> &       pinfo,
	                                 std::vector<LocalData<D>> &local_datas) const override
	{
		called = true;
		local_calls.emplace_back(&pinfo, local_datas.size());
	}

	CallMockMPIGhostFiller(std::shared_ptr<const Domain<D>> domain, int num_components,
	                       GhostFillingType fill_type)
	: MPIGhostFiller<D>(domain, fill_type), num_components(num_components)
	{
	}

	/**
	 * @brief was any of the functions called
	 */
	mutable bool called = false;

	void checkNbrCalls()
	{
		// remove from this collection the calls
		auto remaining_nbr_calls = nbr_calls;

		auto check_for_nbr_call
		= [&](const std::tuple<const PatchInfo<D> *, Side<D>, const NbrType,
		                       const Orthant<D - 1>, int, int> &call) {
			  // check if call was made for neighbor
			  auto found_call
			  = std::find(remaining_nbr_calls.begin(), remaining_nbr_calls.end(), call);

			  CHECK(found_call != remaining_nbr_calls.end());

			  if (found_call != remaining_nbr_calls.end()) {
				  remaining_nbr_calls.erase(found_call);
			  }
		  };
		for (const PatchInfo<D> &patch : this->domain->getPatchInfoVector()) {
			INFO("id: " << patch.id);
			std::string starts = "starts: ";
			for (size_t i = 0; i < D; i++) {
				starts = starts + " " + std::to_string(patch.starts[i]);
			}
			INFO(starts);

			for (auto side : Side<D>::getValues()) {
				INFO("side: " << side);
				if (patch.hasNbr(side)) {
					switch (patch.getNbrType(side)) {
						case NbrType::Normal: {
							INFO("NbrType: Normal");

							auto call = std::make_tuple(&patch, side, NbrType::Normal, Orthant<D - 1>::null(), num_components, num_components);
							check_for_nbr_call(call);

						} break;
						case NbrType::Fine: {
							INFO("NbrType: Fine");

							for (auto orthant : Orthant<D - 1>::getValues()) {
								INFO("Orthant: " << orthant.getIndex());
								auto call = std::make_tuple(&patch, side, NbrType::Fine, orthant, num_components, num_components);
								check_for_nbr_call(call);
							}
						} break;
						case NbrType::Coarse: {
							INFO("NbrType: Coarse");

							auto orthant = patch.getCoarseNbrInfo(side).orth_on_coarse;

							INFO("Orthant: " << orthant.getIndex());

							auto call = std::make_tuple(&patch, side, NbrType::Coarse, orthant, num_components, num_components);
							check_for_nbr_call(call);

						} break;
						default:
							CHECK(false);
					}
				}
			}
		}
		CHECK(remaining_nbr_calls.empty());
		INFO("REMAINING CALL");
		for (auto call : remaining_nbr_calls) {
			auto patch = *std::get<0>(call);
			INFO("id: " << patch.id);
			std::string starts = "starts: ";
			for (size_t i = 0; i < D; i++) {
				starts = starts + " " + std::to_string(patch.starts[i]);
			}
			INFO(starts);
			INFO("side: " << std::get<1>(call));
			std::string nbrtype = "";
			switch (std::get<2>(call)) {
				case NbrType::Normal:
					nbrtype = "normal";
					break;
				case NbrType::Fine:
					nbrtype = "fine";
					break;
				case NbrType::Coarse:
					nbrtype = "coarse";
					break;
				default:
					nbrtype = "???";
			}
			INFO("nbrtype: " << nbrtype);
			INFO("orthant: " << std::get<3>(call).getIndex());
			CHECK(false);
		}
	}
	void checkEdgeNbrCalls()
	{
		if constexpr (D == 3) {
			// remove from this collection the calls
			auto remaining_nbr_calls = edge_nbr_calls;

			auto check_for_nbr_call
			= [&](const std::tuple<const PatchInfo<D> *, Edge, const NbrType,
			                       const Orthant<1>, int, int> &call) {
				  // check if call was made for neighbor
				  auto found_call
				  = std::find(remaining_nbr_calls.begin(), remaining_nbr_calls.end(), call);

				  CHECK(found_call != remaining_nbr_calls.end());

				  if (found_call != remaining_nbr_calls.end()) {
					  remaining_nbr_calls.erase(found_call);
				  }
			  };
			for (const PatchInfo<D> &patch : this->domain->getPatchInfoVector()) {
				INFO("id: " << patch.id);
				std::string starts = "starts: ";
				for (size_t i = 0; i < D; i++) {
					starts = starts + " " + std::to_string(patch.starts[i]);
				}
				INFO(starts);

				for (auto edge : Edge::getValues()) {
					INFO("side: " << edge);
					if (patch.hasNbr(edge)) {
						switch (patch.getNbrType(edge)) {
							case NbrType::Normal: {
								INFO("NbrType: Normal");

								auto call = std::make_tuple(&patch, edge, NbrType::Normal, Orthant<1>::null(), num_components, num_components);
								check_for_nbr_call(call);

							} break;
							case NbrType::Fine: {
								INFO("NbrType: Fine");

								for (auto orthant : Orthant<1>::getValues()) {
									INFO("Orthant: " << orthant.getIndex());
									auto call = std::make_tuple(&patch, edge, NbrType::Fine, orthant, num_components, num_components);
									check_for_nbr_call(call);
								}
							} break;
							case NbrType::Coarse: {
								INFO("NbrType: Coarse");

								auto orthant = patch.getCoarseNbrInfo(edge).orth_on_coarse;

								INFO("Orthant: " << orthant.getIndex());

								auto call = std::make_tuple(&patch, edge, NbrType::Coarse, orthant, num_components, num_components);
								check_for_nbr_call(call);

							} break;
							default:
								CHECK(false);
						}
					}
				}
			}
			CHECK(remaining_nbr_calls.empty());
			INFO("REMAINING CALL");
			for (auto call : remaining_nbr_calls) {
				auto patch = *std::get<0>(call);
				INFO("id: " << patch.id);
				std::string starts = "starts: ";
				for (size_t i = 0; i < D; i++) {
					starts = starts + " " + std::to_string(patch.starts[i]);
				}
				INFO(starts);
				INFO("side: " << std::get<1>(call));
				std::string nbrtype = "";
				switch (std::get<2>(call)) {
					case NbrType::Normal:
						nbrtype = "normal";
						break;
					case NbrType::Fine:
						nbrtype = "fine";
						break;
					case NbrType::Coarse:
						nbrtype = "coarse";
						break;
					default:
						nbrtype = "???";
				}
				INFO("nbrtype: " << nbrtype);
				INFO("orthant: " << std::get<3>(call).getIndex());
				CHECK(false);
			}
		}
	}
	void checkCornerNbrCalls()
	{
		// remove from this collection the calls
		auto remaining_nbr_calls = corner_nbr_calls;

		auto check_for_nbr_call
		= [&](const std::tuple<const PatchInfo<D> *, Corner<D>, const NbrType,
		                       int, int> &call) {
			  // check if call was made for neighbor
			  auto found_call
			  = std::find(remaining_nbr_calls.begin(), remaining_nbr_calls.end(), call);

			  CHECK(found_call != remaining_nbr_calls.end());

			  if (found_call != remaining_nbr_calls.end()) {
				  remaining_nbr_calls.erase(found_call);
			  }
		  };
		for (const PatchInfo<D> &patch : this->domain->getPatchInfoVector()) {
			INFO("id: " << patch.id);
			std::string starts = "starts: ";
			for (size_t i = 0; i < D; i++) {
				starts = starts + " " + std::to_string(patch.starts[i]);
			}
			INFO(starts);

			for (auto corner : Corner<D>::getValues()) {
				INFO("side: " << corner);
				if (patch.hasNbr(corner)) {
					switch (patch.getNbrType(corner)) {
						case NbrType::Normal: {
							INFO("NbrType: Normal");

							auto call = std::make_tuple(&patch, corner, NbrType::Normal, num_components, num_components);
							check_for_nbr_call(call);

						} break;
						case NbrType::Fine: {
							INFO("NbrType: Fine");

							auto call = std::make_tuple(&patch, corner, NbrType::Fine, num_components, num_components);
							check_for_nbr_call(call);
						} break;
						case NbrType::Coarse: {
							INFO("NbrType: Coarse");

							auto call = std::make_tuple(&patch, corner, NbrType::Coarse, num_components, num_components);
							check_for_nbr_call(call);

						} break;
						default:
							CHECK(false);
					}
				}
			}
		}
		CHECK(remaining_nbr_calls.empty());
		INFO("REMAINING CALL");
		for (auto call : remaining_nbr_calls) {
			auto patch = *std::get<0>(call);
			INFO("id: " << patch.id);
			std::string starts = "starts: ";
			for (size_t i = 0; i < D; i++) {
				starts = starts + " " + std::to_string(patch.starts[i]);
			}
			INFO(starts);
			INFO("side: " << std::get<1>(call));
			std::string nbrtype = "";
			switch (std::get<2>(call)) {
				case NbrType::Normal:
					nbrtype = "normal";
					break;
				case NbrType::Fine:
					nbrtype = "fine";
					break;
				case NbrType::Coarse:
					nbrtype = "coarse";
					break;
				default:
					nbrtype = "???";
			}
			INFO("nbrtype: " << nbrtype);
			CHECK(false);
		}
	}
	void checkLocalCalls()
	{
		// remove from this collection the calls
		auto remaining_local_calls = local_calls;

		for (const PatchInfo<D> &patch : this->domain->getPatchInfoVector()) {
			INFO("id: " << patch.id);
			std::string starts = "starts: ";
			for (size_t i = 0; i < D; i++) {
				starts = starts + " " + std::to_string(patch.starts[i]);
			}
			INFO(starts);

			// check for local call
			auto found_call = std::find(remaining_local_calls.begin(), remaining_local_calls.end(),
			                            std::make_tuple(&patch, num_components));

			CHECK(found_call != remaining_local_calls.end());

			if (found_call != remaining_local_calls.end()) {
				remaining_local_calls.erase(found_call);
			}
		}
		CHECK(remaining_local_calls.empty());
	}
	void checkCalls()
	{
		INFO("NumLocalPatches: " << this->domain->getNumLocalPatches());
		INFO("NumGlobalPatches: " << this->domain->getNumGlobalPatches());

		switch (this->fill_type) {
			case GhostFillingType::Corners:
				checkCornerNbrCalls();
			case GhostFillingType::Edges:
				checkEdgeNbrCalls();
			case GhostFillingType::Faces:
				checkNbrCalls();
		}
		checkLocalCalls();
	}
};
template <int D>
class ExchangeMockMPIGhostFiller : public MPIGhostFiller<D>
{
	public:
	void fillGhostCellsForNbrPatch(const PatchInfo<D> &             pinfo,
	                               const std::vector<LocalData<D>> &local_datas,
	                               std::vector<LocalData<D>> &nbr_datas, Side<D> side,
	                               NbrType nbr_type, Orthant<D - 1> orthant) const override
	{
		for (size_t c = 0; c < nbr_datas.size(); c++) {
			int index = 0;
			for (int i = 0; i < pinfo.num_ghost_cells; i++) {
				LocalData<D - 1> slice = nbr_datas[c].getSliceOn(side.opposite(), {-i - 1});
				nested_loop<D - 1>(slice.getStart(), slice.getEnd(),
				                   [&](const std::array<int, D - 1> &coord) {
					                   slice[coord] += pinfo.id + index + c;
					                   index++;
				                   });
			}
		}
	}

	void fillGhostCellsForEdgeNbrPatch(const PatchInfo<D> &             pinfo,
	                                   const std::vector<LocalData<D>> &local_datas,
	                                   std::vector<LocalData<D>> &nbr_datas, Edge edge,
	                                   NbrType nbr_type, Orthant<1> orthant) const override
	{
		if constexpr (D == 3) {
			for (size_t c = 0; c < nbr_datas.size(); c++) {
				int index = 0;
				for (int j = 0; j < pinfo.num_ghost_cells; j++) {
					for (int i = 0; i < pinfo.num_ghost_cells; i++) {
						LocalData<1> slice = nbr_datas[c].getSliceOn(edge.opposite(), {-1 - i, -1 - j});
						nested_loop<1>(slice.getStart(), slice.getEnd(),
						               [&](const std::array<int, 1> &coord) {
							               slice[coord] += pinfo.id + index + c;
							               index++;
						               });
					}
				}
			}
		}
	}

	void fillGhostCellsForCornerNbrPatch(const PatchInfo<D> &             pinfo,
	                                     const std::vector<LocalData<D>> &local_datas,
	                                     std::vector<LocalData<D>> &nbr_datas, Corner<D> corner,
	                                     NbrType nbr_type) const override
	{
		for (size_t c = 0; c < nbr_datas.size(); c++) {
			int                index = 0;
			std::array<int, D> start = nbr_datas[c].getGhostStart();
			std::array<int, D> end;
			end.fill(-1);
			nested_loop<D>(start, end, [&](const std::array<int, D> &offset) {
				nbr_datas[c].getSliceOn(corner.opposite(), offset)[{}] += pinfo.id + index + c;
				index++;
			});
		}
	}

	void fillGhostCellsForLocalPatch(const PatchInfo<D> &       pinfo,
	                                 std::vector<LocalData<D>> &local_data) const override {}

	ExchangeMockMPIGhostFiller(std::shared_ptr<const Domain<D>> domain_in, GhostFillingType fill_type)
	: MPIGhostFiller<D>(domain_in, fill_type) {}

	void checkInterior(std::shared_ptr<const Vector<D>> vec)
	{
		for (auto pinfo : this->domain->getPatchInfoVector()) {
			INFO("id: " << pinfo.id);
			std::string starts = "start: ";
			for (size_t i = 0; i < D; i++) {
				starts += " " + std::to_string(pinfo.starts[i]);
			}
			INFO(starts);
			std::string ns = "ns: ";
			for (size_t i = 0; i < D; i++) {
				ns += " " + std::to_string(pinfo.ns[i]);
			}
			INFO(ns);
			INFO("num_ghost_cells: " << pinfo.num_ghost_cells);
			for (int c = 0; c < vec->getNumComponents(); c++) {
				INFO("c: " << c);
				auto data = vec->getLocalData(c, pinfo.local_index);
				// check that vector was not modified on the interior
				nested_loop<D>(data.getStart(), data.getEnd(),
				               [&](const std::array<int, D> &coord) {
					               std::string coord_str = "coord: ";
					               for (size_t i = 0; i < D; i++) {
						               coord_str += " " + std::to_string(coord[i]);
					               }
					               INFO(coord_str);
					               CHECK(data[coord] == pinfo.id);
				               });
			}
		}
	}
	void checkCorners(std::shared_ptr<const Vector<D>> vec)
	{
		for (auto pinfo : this->domain->getPatchInfoVector()) {
			INFO("id: " << pinfo.id);
			std::string starts = "start: ";
			for (size_t i = 0; i < D; i++) {
				starts += " " + std::to_string(pinfo.starts[i]);
			}
			INFO(starts);
			std::string ns = "ns: ";
			for (size_t i = 0; i < D; i++) {
				ns += " " + std::to_string(pinfo.ns[i]);
			}
			INFO(ns);
			INFO("num_ghost_cells: " << pinfo.num_ghost_cells);
			for (int c = 0; c < vec->getNumComponents(); c++) {
				INFO("c: " << c);
				auto data = vec->getLocalData(c, pinfo.local_index);
				// check the ghost cells
				for (Corner<D> corner : Corner<D>::getValues()) {
					INFO("Corner: " << corner);

					if (pinfo.hasNbr(corner)) {
						switch (pinfo.getNbrType(corner)) {
							case NbrType::Normal: {
								INFO("NbrType: Normal");

								// value should be id of neighbor + index +c
								int  index   = 0;
								auto nbrinfo = pinfo.getNormalNbrInfo(corner);

								std::array<int, D> start = data.getGhostStart();
								std::array<int, D> end;
								end.fill(-1);
								nested_loop<D>(start, end, [&](std::array<int, D> &coord) {
									std::string coord_str = "coord: ";
									for (size_t i = 0; i < D; i++) {
										coord_str += " " + std::to_string(coord[i]);
									}
									INFO(coord_str);
									CHECK(data.getSliceOn(corner, coord)[{}] == nbrinfo.id + index + c);
									index++;
								});
							} break;
							case NbrType::Fine: {
								INFO("NbrType: Fine");

								int  ids     = 0;
								auto nbrinfo = pinfo.getFineNbrInfo(corner);

								// value should be id of neighbors + index
								int index = 0;

								std::array<int, D> start = data.getGhostStart();
								std::array<int, D> end;
								end.fill(-1);
								nested_loop<D>(start, end, [&](std::array<int, D> &coord) {
									std::string coord_str = "coord: ";
									for (size_t i = 0; i < D; i++) {
										coord_str += " " + std::to_string(coord[i]);
									}
									INFO(coord_str);
									CHECK(data.getSliceOn(corner, coord)[{}] == nbrinfo.ids[0] + index + c);
									index++;
								});
							} break;
							case NbrType::Coarse: {
								INFO("NbrType: Coarse");

								// value should be id of neighbor + index
								int  index   = 0;
								auto nbrinfo = pinfo.getCoarseNbrInfo(corner);

								std::array<int, D> start = data.getGhostStart();
								std::array<int, D> end;
								end.fill(-1);
								nested_loop<D>(start, end, [&](std::array<int, D> &coord) {
									std::string coord_str = "coord: ";
									for (size_t i = 0; i < D; i++) {
										coord_str += " " + std::to_string(coord[i]);
									}
									INFO(coord_str);
									CHECK(data.getSliceOn(corner, coord)[{}] == nbrinfo.id + index + c);
									index++;
								});

							} break;
						}
					} else {
						// values should be zero on ghost cells
						INFO("Physical Boundary");

						std::array<int, D> start = data.getGhostStart();
						std::array<int, D> end;
						end.fill(-1);
						nested_loop<D>(start, end, [&](std::array<int, D> &coord) {
							std::string coord_str = "coord: ";
							for (size_t i = 0; i < D; i++) {
								coord_str += " " + std::to_string(coord[i]);
							}
							INFO(coord_str);
							CHECK(data.getSliceOn(corner, coord)[{}] == 0);
						});
					}
				}
			}
		}
	}

	void
	checkEdges(std::shared_ptr<const Vector<D>> vec)
	{
		if constexpr (D == 3) {
			for (auto pinfo : this->domain->getPatchInfoVector()) {
				INFO("id: " << pinfo.id);
				std::string starts = "start: ";
				for (size_t i = 0; i < D; i++) {
					starts += " " + std::to_string(pinfo.starts[i]);
				}
				INFO(starts);
				std::string ns = "ns: ";
				for (size_t i = 0; i < D; i++) {
					ns += " " + std::to_string(pinfo.ns[i]);
				}
				INFO(ns);
				INFO("num_ghost_cells: " << pinfo.num_ghost_cells);
				for (int c = 0; c < vec->getNumComponents(); c++) {
					INFO("c: " << c);
					auto data = vec->getLocalData(c, pinfo.local_index);
					// check the ghost cells
					for (Edge e : Edge::getValues()) {
						INFO("Edge: " << e);

						if (pinfo.hasNbr(e)) {
							switch (pinfo.getNbrType(e)) {
								case NbrType::Normal: {
									INFO("NbrType: Normal");

									// value should be id of neighbor + index + c
									int  index   = 0;
									auto nbrinfo = pinfo.getNormalNbrInfo(e);
									for (int j = 0; j < pinfo.num_ghost_cells; j++) {
										for (int i = 0; i < pinfo.num_ghost_cells; i++) {
											auto slice = data.getSliceOn(e, {-i - 1, -j - 1});

											nested_loop<1>(
											slice.getStart(), slice.getEnd(),
											[&](std::array<int, 1> &coord) {
												std::string coord_str = "coord: ";
												for (size_t i = 0; i < 1; i++) {
													coord_str += " " + std::to_string(coord[i]);
												}
												INFO(coord_str);
												CHECK(slice[coord] == nbrinfo.id + index + c);
												index++;
											});
										}
									}

								} break;
								case NbrType::Fine: {
									INFO("NbrType: Fine");

									int  ids     = 0;
									auto nbrinfo = pinfo.getFineNbrInfo(e);
									for (size_t i = 0; i < nbrinfo.ids.size(); i++) {
										ids += nbrinfo.ids[i];
									}

									// value should be id of neighbors + index
									int index = 0;
									for (int j = 0; j < pinfo.num_ghost_cells; j++) {
										for (int i = 0; i < pinfo.num_ghost_cells; i++) {
											auto slice = data.getSliceOn(e, {-1 - i, -1 - j});
											nested_loop<1>(slice.getStart(), slice.getEnd(), [&](std::array<int, 1> &coord) {
												std::string coord_str = "coord: ";
												for (size_t i = 0; i < 1; i++) {
													coord_str += " " + std::to_string(coord[i]);
												}
												INFO(coord_str);
												CHECK(slice[coord] == ids + 2 * (index + c));
												index++;
											});
										}
									}
								} break;
								case NbrType::Coarse: {
									INFO("NbrType: Coarse");

									// value should be id of neighbor + index
									int  index   = 0;
									auto nbrinfo = pinfo.getCoarseNbrInfo(e);
									for (int j = 0; j < pinfo.num_ghost_cells; j++) {
										for (int i = 0; i < pinfo.num_ghost_cells; i++) {
											auto slice = data.getSliceOn(e, {-1 - i, -1 - j});
											nested_loop<1>(slice.getStart(), slice.getEnd(), [&](std::array<int, 1> &coord) {
												std::string coord_str = "coord: ";
												for (size_t i = 0; i < 1; i++) {
													coord_str += " " + std::to_string(coord[i]);
												}
												INFO(coord_str);
												CHECK(slice[coord] == nbrinfo.id + index + c);
												index++;
											});
										}
									}

								} break;
							}
						} else {
							// values should be zero on ghost cells
							INFO("Physical Boundary");

							for (int j = 0; j < pinfo.num_ghost_cells; j++) {
								for (int i = 0; i < pinfo.num_ghost_cells; i++) {
									auto slice = data.getSliceOn(e, {-1 - i, -1 - j});
									nested_loop<1>(slice.getStart(), slice.getEnd(), [&](std::array<int, 1> &coord) {
										INFO("coord: ");
										for (size_t i = 0; i < 1; i++) {
											INFO(coord[i]);
										}
										CHECK(slice[coord] == 0);
									});
								}
							}
						}
					}
				}
			}
		}
	}

	void checkFaces(std::shared_ptr<const Vector<D>> vec)
	{
		if constexpr (D == 3) {
			for (auto pinfo : this->domain->getPatchInfoVector()) {
				INFO("id: " << pinfo.id);
				std::string starts = "start: ";
				for (size_t i = 0; i < D; i++) {
					starts += " " + std::to_string(pinfo.starts[i]);
				}
				INFO(starts);
				std::string ns = "ns: ";
				for (size_t i = 0; i < D; i++) {
					ns += " " + std::to_string(pinfo.ns[i]);
				}
				INFO(ns);
				INFO("num_ghost_cells: " << pinfo.num_ghost_cells);
				for (int c = 0; c < vec->getNumComponents(); c++) {
					INFO("c: " << c);
					auto data = vec->getLocalData(c, pinfo.local_index);
					// check the ghost cells
					for (Side<D> s : Side<D>::getValues()) {
						INFO("Side: " << s);

						if (pinfo.hasNbr(s)) {
							switch (pinfo.getNbrType(s)) {
								case NbrType::Normal: {
									INFO("NbrType: Normal");

									// value should be id of neighbor + index +c
									int  index   = 0;
									auto nbrinfo = pinfo.getNormalNbrInfo(s);
									for (int i = 0; i < pinfo.num_ghost_cells; i++) {
										auto slice = data.getSliceOn(s, {-i - 1});

										nested_loop<D - 1>(
										slice.getStart(), slice.getEnd(),
										[&](std::array<int, D - 1> &coord) {
											std::string coord_str = "coord: ";
											for (size_t i = 0; i < D - 1; i++) {
												coord_str += " " + std::to_string(coord[i]);
											}
											INFO(coord_str);
											CHECK(slice[coord] == nbrinfo.id + index + c);
											index++;
										});
									}

								} break;
								case NbrType::Fine: {
									INFO("NbrType: Fine");

									int  ids     = 0;
									auto nbrinfo = pinfo.getFineNbrInfo(s);
									for (size_t i = 0; i < nbrinfo.ids.size(); i++) {
										ids += nbrinfo.ids[i];
									}

									// value should be id of neighbors + index
									int index = 0;
									for (int i = 0; i < pinfo.num_ghost_cells; i++) {
										auto slice = data.getSliceOn(s, {-i - 1});

										nested_loop<D - 1>(
										slice.getStart(), slice.getEnd(),
										[&](std::array<int, D - 1> &coord) {
											std::string coord_str = "coord: ";
											for (size_t i = 0; i < D - 1; i++) {
												coord_str += " " + std::to_string(coord[i]);
											}
											INFO(coord_str);
											CHECK(slice[coord] == ids + (1 << (D - 1)) * (index + c));
											index++;
										});
									}
								} break;
								case NbrType::Coarse: {
									INFO("NbrType: Coarse");

									// value should be id of neighbor + index
									int  index   = 0;
									auto nbrinfo = pinfo.getCoarseNbrInfo(s);
									for (int i = 0; i < pinfo.num_ghost_cells; i++) {
										auto slice = data.getSliceOn(s, {-i - 1});

										nested_loop<D - 1>(
										slice.getStart(), slice.getEnd(),
										[&](std::array<int, D - 1> &coord) {
											std::string coord_str = "coord: ";
											for (size_t i = 0; i < D - 1; i++) {
												coord_str += " " + std::to_string(coord[i]);
											}
											INFO(coord_str);
											CHECK(slice[coord] == nbrinfo.id + index + c);
											index++;
										});
									}

								} break;
							}
						} else {
							// values should be zero on ghost cells
							INFO("Physical Boundary");

							for (int i = 0; i < pinfo.num_ghost_cells; i++) {
								auto slice = data.getSliceOn(s, {-i - 1});

								nested_loop<D - 1>(slice.getStart(), slice.getEnd(),
								                   [&](std::array<int, D - 1> &coord) {
									                   INFO("coord: ");
									                   for (size_t i = 0; i < D - 1; i++) {
										                   INFO(coord[i]);
									                   }
									                   CHECK(slice[coord] == 0);
								                   });
							}
						}
					}
				}
			}
		}
	}

	void checkVector(std::shared_ptr<const Vector<D>> vec)
	{
		INFO("NumLocalPatches: " << this->domain->getNumLocalPatches());
		INFO("NumGlobalPatches: " << this->domain->getNumGlobalPatches());
		checkInterior(vec);
		switch (this->fill_type) {
			case GhostFillingType::Corners:
				checkCorners(vec);
			case GhostFillingType::Edges:
				checkEdges(vec);
			case GhostFillingType::Faces:
				checkFaces(vec);
		}
	}
};
} // namespace ThunderEgg
#endif