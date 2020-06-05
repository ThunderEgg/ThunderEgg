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
constexpr auto cross_mesh_file   = "mesh_inputs/2d_uniform_2x2_refined_cross_mpi1.json";

template <size_t D> class CallMockMPIGhostFiller : public MPIGhostFiller<D>
{
	private:
	mutable std::list<std::tuple<const std::vector<Side<D>>, const NbrType, const Orthant<D>>>
	nbr_calls;

	void fillGhostCellsForNbrPatch(const LocalData<D> local_data, LocalData<D> nbr_data,
	                               const std::vector<Side<D>> &side, const NbrType nbr_type,
	                               const Orthant<D> orthant) const override
	{
		called = true;
		nbr_calls.emplace_back(side, nbr_type, orthant);
	}

	mutable std::list<std::shared_ptr<const PatchInfo<D>>> local_calls;

	void fillGhostCellsForLocalPatch(std::shared_ptr<const PatchInfo<D>> pinfo,
	                                 LocalData<D>                        local_data) const override
	{
		called = true;
		local_calls.push_back(pinfo);
	}

	public:
	CallMockMPIGhostFiller(std::shared_ptr<Domain<D>> domain_in, int side_cases_in)
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
		= [&](const std::tuple<const std::vector<Side<D>>, const NbrType, const Orthant<D>> &call) {
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
			INFO("start: ");
			for (size_t i = 0; i < D; i++) {
				INFO(patch->starts[i]);
			}

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

							std::vector<Side<D>> side_vector;
							side_vector.push_back(side);

							auto call = make_tuple(side_vector, NbrType::Normal, Orthant<D>(0));
							check_for_nbr_call(call);

						} break;
						case NbrType::Fine: {
							INFO("NbrType: Fine");

							std::vector<Side<D>> side_vector;
							side_vector.push_back(side);

							for (auto orthant : Orthant<D>::getValuesOnSide(side)) {
								INFO("Orthant: " << orthant.toInt());
								auto call = make_tuple(side_vector, NbrType::Fine, orthant);
								check_for_nbr_call(call);
							}
						} break;
						case NbrType::Coarse: {
							INFO("NbrType: Coarse");

							std::vector<Side<D>> side_vector;
							side_vector.push_back(side);

							auto orthant = Orthant<D>::getValuesOnSide(
							side)[patch->getCoarseNbrInfo(side).orth_on_coarse.toInt()];

							INFO("Orthant: " << orthant.toInt());

							auto call = make_tuple(side_vector, NbrType::Fine, orthant);
							check_for_nbr_call(call);

						} break;
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

	CHECK(mgf.called == false);
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