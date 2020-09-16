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

#ifndef THUNDEREGG_SCHUR_PATCHIFACESCATTER_H
#define THUNDEREGG_SCHUR_PATCHIFACESCATTER_H
#include <ThunderEgg/RuntimeError.h>
#include <ThunderEgg/Schur/InterfaceDomain.h>
#include <ThunderEgg/ValVector.h>
namespace ThunderEgg
{
namespace Schur
{
/**
 * @brief Scatters between a global Schur compliment vector and a local patch iface vector
 *
 * The scatters functions are split with a Start and Finish, this allows for local computation to
 * occur while the communicating
 *
 * @tparam D the number of cartesian dimensions on a patch
 */
template <int D> class PatchIfaceScatter
{
	public:
	/**
	 * @brief Construct a new PatchIfaceScatter object
	 *
	 * @param iface_domain the InterfaceDomain
	 */
	PatchIfaceScatter(std::shared_ptr<InterfaceDomain<D>> iface_domain) {}
	/**
	 * @brief Get a nw local patch iface vector
	 *
	 * @return std::shared_ptr<Vector<D - 1>> the new vector
	 */
	std::shared_ptr<Vector<D - 1>> getNewLocalPatchIfaceVector()
	{
		return nullptr;
	}
	/**
	 * @brief Start the scatter from the global Schur compliment vector to the local patch iface
	 * vector
	 *
	 * Will throw an exception if any communcation is in progress
	 *
	 * @param global_vector the global Schur compliment vector
	 * @param local_patch_iface_vector the the local patch iface vector
	 */
	void scatterStart(std::shared_ptr<const Vector<D - 1>> global_vector,
	                  std::shared_ptr<Vector<D - 1>>       local_patch_iface_vector)
	{
	}
	/**
	 * @brief Finish the scatter from the global Schur compliment vector to the local patch iface
	 * vector
	 *
	 * Will throw an exception if different vectors are passed
	 * from when scatterStart was called
	 *
	 * @param global_vector the global Schur compliment vector
	 * @param local_patch_iface_vector the the local patch iface vector
	 */
	void scatterFinish(std::shared_ptr<const Vector<D - 1>> global_vector,
	                   std::shared_ptr<Vector<D - 1>>       local_patch_iface_vector)
	{
	}
};
extern template class PatchIfaceScatter<2>;
extern template class PatchIfaceScatter<3>;
} // namespace Schur
} // namespace ThunderEgg
#endif
