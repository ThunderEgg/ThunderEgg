#include "catch.hpp"
#include "utils/DomainReader.h"
#include <Thunderegg/DomainTools.h>
#include <Thunderegg/MPIGhostFiller.h>
#include <Thunderegg/ValVector.h>
#include <list>
using namespace std;
using namespace Thunderegg;

constexpr auto single_mesh_file  = "mesh_inputs/2d_uniform_2x2_mpi1.json";
constexpr auto refined_mesh_file = "mesh_inputs/2d_uniform_2x2_refined_nw_mpi1.json";
constexpr auto cross_mesh_file   = "mesh_inputs/2d_uniform_8x8_refined_cross_mpi1.json";

template <size_t D> class CallMockMPIGhostFiller : public MPIGhostFiller<D>
{
	private:
	mutable std::list<
	std::tuple<std::shared_ptr<const PatchInfo<D>>, const Side<D>, const NbrType, const Orthant<D>>>
	nbr_calls;

	void fillGhostCellsForNbrPatch(std::shared_ptr<const PatchInfo<D>> pinfo,
	                               const LocalData<D> local_data, const LocalData<D> nbr_data,
	                               const Side<D> side, const NbrType nbr_type,
	                               const Orthant<D> orthant) const override
	{
		called = true;
		nbr_calls.emplace_back(pinfo, side, nbr_type, orthant);
	}

	mutable std::list<std::shared_ptr<const PatchInfo<D>>> local_calls;

	void fillGhostCellsForLocalPatch(std::shared_ptr<const PatchInfo<D>> pinfo,
	                                 const LocalData<D>                  local_data) const override
	{
		called = true;
		local_calls.push_back(pinfo);
	}

	public:
	CallMockMPIGhostFiller(std::shared_ptr<const Domain<D>> domain_in, int side_cases_in)
	: MPIGhostFiller<D>(domain_in, side_cases_in)
	{
	}
	/**
	 * @brief was any of the functions called
	 */
	mutable bool called = false;

	void checkCalls()
	{
		INFO("NumLocalPatches: " << this->domain->getNumLocalPatches());
		INFO("NumGlobalPatches: " << this->domain->getNumGlobalPatches());

		// remove from this collection the calls
		auto remaining_nbr_calls   = nbr_calls;
		auto remaining_local_calls = local_calls;

		auto check_for_nbr_call
		= [&](const std::tuple<std::shared_ptr<const PatchInfo<D>>, const Side<D>, const NbrType,
		                       const Orthant<D>> &call) {
			  // check if call was made for neighbor
			  auto found_call
			  = std::find(remaining_nbr_calls.begin(), remaining_nbr_calls.end(), call);

			  CHECK(found_call != remaining_nbr_calls.end());

			  if (found_call != remaining_nbr_calls.end()) {
				  remaining_nbr_calls.erase(found_call);
			  }
		  };
		for (auto patch : this->domain->getPatchInfoVector()) {
			INFO("id: " << patch->id);
			std::string starts = "starts: ";
			for (size_t i = 0; i < D; i++) {
				starts = starts + " " + std::to_string(patch->starts[i]);
			}
			INFO(starts);

			// check for local call
			auto found_call
			= std::find(remaining_local_calls.begin(), remaining_local_calls.end(), patch);

			CHECK(found_call != remaining_local_calls.end());

			if (found_call != remaining_local_calls.end()) {
				remaining_local_calls.erase(found_call);
			}

			for (auto side : Side<D>::getValues()) {
				INFO("side: " << side);
				if (patch->hasNbr(side)) {
					switch (patch->getNbrType(side)) {
						case NbrType::Normal: {
							INFO("NbrType: Normal");

							auto call = make_tuple(patch, side, NbrType::Normal, Orthant<D>(0));
							check_for_nbr_call(call);

						} break;
						case NbrType::Fine: {
							INFO("NbrType: Fine");

							for (auto orthant : Orthant<D>::getValuesOnSide(side)) {
								INFO("Orthant: " << orthant.toInt());
								auto call = make_tuple(patch, side, NbrType::Fine, orthant);
								check_for_nbr_call(call);
							}
						} break;
						case NbrType::Coarse: {
							INFO("NbrType: Coarse");

							auto orthant = Orthant<D>::getValuesOnSide(
							side.opposite())[patch->getCoarseNbrInfo(side).orth_on_coarse.toInt()];

							INFO("Orthant: " << orthant.toInt());

							auto call = make_tuple(patch, side, NbrType::Coarse, orthant);
							check_for_nbr_call(call);

						} break;
						default:
							CHECK(false);
					}
				}
			}
		}
		CHECK(remaining_local_calls.empty());
		CHECK(remaining_nbr_calls.empty());
		INFO("REMAINING CALL")
		for (auto call : remaining_nbr_calls) {
			auto patch = std::get<0>(call);
			INFO("id: " << patch->id);
			std::string starts = "starts: ";
			for (size_t i = 0; i < D; i++) {
				starts = starts + " " + std::to_string(patch->starts[i]);
			}
			INFO(starts);
			INFO("side: " << std::get<1>(call));
			string nbrtype = "";
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
			INFO("orthant: " << std::get<3>(call).toInt());
			CHECK(false);
		}
	}
};
template <size_t D> class ExchangeMockMPIGhostFiller : public MPIGhostFiller<D>
{
	private:
	void fillGhostCellsForNbrPatch(std::shared_ptr<const PatchInfo<D>> pinfo,
	                               const LocalData<D> local_data, const LocalData<D> nbr_data,
	                               const Side<D> side, const NbrType nbr_type,
	                               const Orthant<D> orthant) const override
	{
		int index = 0;
		for (int i = 0; i < pinfo->num_ghost_cells; i++) {
			LocalData<D - 1> slice = nbr_data.getGhostSliceOnSide(side.opposite(), i + 1);
			nested_loop<D - 1>(slice.getStart(), slice.getEnd(),
			                   [&](const std::array<int, D - 1> &coord) {
				                   slice[coord] += pinfo->id + index;
				                   index++;
			                   });
		}
	}

	void fillGhostCellsForLocalPatch(std::shared_ptr<const PatchInfo<D>> pinfo,
	                                 LocalData<D>                        local_data) const override
	{
	}

	public:
	ExchangeMockMPIGhostFiller(std::shared_ptr<const Domain<D>> domain_in, int side_cases_in)
	: MPIGhostFiller<D>(domain_in, side_cases_in)
	{
	}

	void checkVector(std::shared_ptr<const Vector<D>> vec)
	{
		INFO("NumLocalPatches: " << this->domain->getNumLocalPatches());
		INFO("NumGlobalPatches: " << this->domain->getNumGlobalPatches());
		for (auto pinfo : this->domain->getPatchInfoVector()) {
			INFO("id: " << pinfo->id);
			INFO("start: ");
			for (size_t i = 0; i < D; i++) {
				INFO(pinfo->starts[i]);
			}
			INFO("ns: ");
			for (size_t i = 0; i < D; i++) {
				INFO(pinfo->ns[i]);
			}
			INFO("num_ghost_cells: " << pinfo->num_ghost_cells);

			auto data = vec->getLocalData(pinfo->local_index);
			// check that vector was not modified on the interior
			nested_loop<D>(data.getStart(), data.getEnd(), [&](const std::array<int, D> &coord) {
				std::string coord_str = "coord: ";
				for (size_t i = 0; i < D; i++) {
					coord_str += " " + std::to_string(coord[i]);
				}
				INFO(coord_str);
				CHECK(data[coord] == pinfo->id);
			});
			// check the ghost cells
			for (Side<D> s : Side<D>::getValues()) {
				INFO("Side: " << s);

				if (pinfo->hasNbr(s)) {
					switch (pinfo->getNbrType(s)) {
						case NbrType::Normal: {
							INFO("NbrType: Normal");

							// value should be id of neighbor + index
							int  index   = 0;
							auto nbrinfo = pinfo->getNormalNbrInfo(s);
							for (int i = 0; i < pinfo->num_ghost_cells; i++) {
								auto slice = data.getGhostSliceOnSide(s, i + 1);

								nested_loop<D - 1>(slice.getStart(), slice.getEnd(),
								                   [&](std::array<int, D - 1> &coord) {
									                   std::string coord_str = "coord: ";
									                   for (size_t i = 0; i < D - 1; i++) {
										                   coord_str
										                   += " " + std::to_string(coord[i]);
									                   }
									                   INFO(coord_str);
									                   CHECK(slice[coord] == nbrinfo.id + index);
									                   index++;
								                   });
							}

						} break;
						case NbrType::Fine: {
							INFO("NbrType: Fine");

							int  ids     = 0;
							auto nbrinfo = pinfo->getFineNbrInfo(s);
							for (int i = 0; i < nbrinfo.ids.size(); i++) {
								ids += nbrinfo.ids[i];
							}

							// value should be id of neighbors + index
							int index = 0;
							for (int i = 0; i < pinfo->num_ghost_cells; i++) {
								auto slice = data.getGhostSliceOnSide(s, i + 1);

								nested_loop<D - 1>(
								slice.getStart(), slice.getEnd(),
								[&](std::array<int, D - 1> &coord) {
									std::string coord_str = "coord: ";
									for (size_t i = 0; i < D - 1; i++) {
										coord_str += " " + std::to_string(coord[i]);
									}
									INFO(coord_str);
									CHECK(slice[coord] == ids + (1 << (D - 1)) * index);
									index++;
								});
							}
						} break;
						case NbrType::Coarse: {
							INFO("NbrType: Coarse");

							// value should be id of neighbor + index
							int  index   = 0;
							auto nbrinfo = pinfo->getCoarseNbrInfo(s);
							for (int i = 0; i < pinfo->num_ghost_cells; i++) {
								auto slice = data.getGhostSliceOnSide(s, i + 1);

								nested_loop<D - 1>(slice.getStart(), slice.getEnd(),
								                   [&](std::array<int, D - 1> &coord) {
									                   std::string coord_str = "coord: ";
									                   for (size_t i = 0; i < D - 1; i++) {
										                   coord_str
										                   += " " + std::to_string(coord[i]);
									                   }
									                   INFO(coord_str);
									                   CHECK(slice[coord] == nbrinfo.id + index);
									                   index++;
								                   });
							}

						} break;
					}
				} else {
					// values should be zero on ghost cells
					INFO("Physical Boundary")

					for (int i = 0; i < pinfo->num_ghost_cells; i++) {
						auto slice = data.getGhostSliceOnSide(s, i + 1);

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
};
TEST_CASE("No calls for 1 patch domain", "[MPIGhostFiller]")
{
	auto                  nx        = GENERATE(2);
	auto                  ny        = GENERATE(2);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(single_mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_coarse = domain_reader.getCoarserDomain();

	auto vec = ValVector<2>::GetNewVector(d_coarse);

	CallMockMPIGhostFiller<2> mgf(d_coarse, 1);

	mgf.fillGhost(vec);

	CHECK(mgf.called == true);
	mgf.checkCalls();
}
TEST_CASE("Calls for various domains 1-side cases", "[MPIGhostFiller]")
{
	auto mesh_file
	= GENERATE(as<std::string>{}, single_mesh_file, refined_mesh_file, cross_mesh_file);
	INFO("MESH: " << mesh_file);
	auto                  nx        = GENERATE(2);
	auto                  ny        = GENERATE(2);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine = domain_reader.getFinerDomain();

	auto vec = ValVector<2>::GetNewVector(d_fine);

	CallMockMPIGhostFiller<2> mgf(d_fine, 1);

	mgf.fillGhost(vec);

	CHECK(mgf.called == true);

	mgf.checkCalls();
}
TEST_CASE("Exchange for various domains 1-side cases", "[MPIGhostFiller]")
{
	auto mesh_file
	= GENERATE(as<std::string>{}, single_mesh_file, refined_mesh_file, cross_mesh_file);
	INFO("MESH: " << mesh_file);
	auto                  nx        = GENERATE(2);
	auto                  ny        = GENERATE(2);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine = domain_reader.getFinerDomain();

	auto vec = ValVector<2>::GetNewVector(d_fine);
	for (auto pinfo : d_fine->getPatchInfoVector()) {
		auto data = vec->getLocalData(pinfo->local_index);
		nested_loop<2>(data.getStart(), data.getEnd(),
		               [&](const std::array<int, 2> &coord) { data[coord] = pinfo->id; });
	}

	ExchangeMockMPIGhostFiller<2> mgf(d_fine, 1);

	mgf.fillGhost(vec);

	mgf.checkVector(vec);
}
TEST_CASE("Two Exchanges for various domains 1-side cases", "[MPIGhostFiller]")
{
	auto mesh_file
	= GENERATE(as<std::string>{}, single_mesh_file, refined_mesh_file, cross_mesh_file);
	INFO("MESH: " << mesh_file);
	auto                  nx        = GENERATE(2);
	auto                  ny        = GENERATE(2);
	int                   num_ghost = 1;
	DomainReader<2>       domain_reader(mesh_file, {nx, ny}, num_ghost);
	shared_ptr<Domain<2>> d_fine = domain_reader.getFinerDomain();

	auto vec = ValVector<2>::GetNewVector(d_fine);
	for (auto pinfo : d_fine->getPatchInfoVector()) {
		auto data = vec->getLocalData(pinfo->local_index);
		nested_loop<2>(data.getStart(), data.getEnd(),
		               [&](const std::array<int, 2> &coord) { data[coord] = pinfo->id; });
	}

	ExchangeMockMPIGhostFiller<2> mgf(d_fine, 1);

	mgf.fillGhost(vec);
	mgf.fillGhost(vec);

	mgf.checkVector(vec);
}