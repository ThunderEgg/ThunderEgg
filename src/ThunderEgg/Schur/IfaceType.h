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

#ifndef THUNDEREGG_SCHUR_IFACETYPE_H
#define THUNDEREGG_SCHUR_IFACETYPE_H
#include <ThunderEgg/Orthant.h>
#include <tuple>
namespace ThunderEgg
{
namespace Schur
{
/**
 * @brief An enum-style class that represents interface types
 *
 * @tparam D the number of cartesian dimensions on a patch
 */
template <size_t D> class IfaceType
{
	private:
	/**
	 * @brief the value of the enum
	 */
	unsigned char val = 10;
	/**
	 * @brief Orthant that the interface lies on
	 */
	Orthant<D - 1> orthant = Orthant<D - 1>::null();

	/**
	 * @brief Construct a new IfaceType object
	 *
	 * @param val the value
	 * @param orthant orthant on the coarser patch's side that the finer
	 * patch lies on.
	 */
	IfaceType(unsigned char val, Orthant<D - 1> orthant) : val(val), orthant(orthant) {}

	public:
	/**
	 * @brief An interface on a side of a patch with a neighbor at the same refinement level.
	 *
	 * @return IfaceType<D> the new IfaceType
	 */
	static IfaceType<D> Normal()
	{
		return IfaceType<D>(0, Orthant<D - 1>::null());
	}
	/**
	 * @brief Check if this type is Normal
	 *
	 * @return if it is Normal
	 */
	bool isNormal() const
	{
		return val == 0;
	}
	/**
	 * @brief An interface on a side of a patch with neighbors at a finer refinement level.
	 *
	 * This interface lines up with the cells on the coarser patch.
	 *
	 * @return IfaceType<D> the new IfaceType
	 */
	static IfaceType<D> CoarseToCoarse()
	{
		return IfaceType<D>(1, Orthant<D - 1>::null());
	}
	/**
	 * @brief Check if this type is Normal
	 *
	 * @return if it is Normal
	 */
	bool isCoarseToCoarse() const
	{
		return val == 1;
	}
	/**
	 * @brief An interface on a side of a patch with a neighbor at a coarser refinement level.
	 *
	 * This interface lines up with the cells on the coarser patch.
	 *
	 * @param orthant the orthant of the fine patch
	 * @return IfaceType<D> the new IfaceType
	 */
	static IfaceType<D> FineToCoarse(Orthant<D - 1> orth_on_coarse)
	{
		return IfaceType<D>(2, orth_on_coarse);
	}
	/**
	 * @brief Check if this type is Normal
	 *
	 * @return if it is Normal
	 */
	bool isFineToCoarse() const
	{
		return val == 2;
	}
	/**
	 * @brief An interface on a side of a patch with a neighbor at a coarser refinement level.
	 *
	 * This interface lines up with the cells on the finer patch.
	 * The orthant value should be set to the orthant on the coarser patch's side that the finer
	 * patch lies on.
	 *
	 * @param orthant the orthant of the fine patch
	 * @return IfaceType<D> the new IfaceType
	 */
	static IfaceType<D> FineToFine(Orthant<D - 1> orth_on_coarse)
	{
		return IfaceType<D>(3, orth_on_coarse);
	}
	/**
	 * @brief Check if this type is Normal
	 *
	 * @return if it is Normal
	 */
	bool isFineToFine() const
	{
		return val == 3;
	}
	/**
	 * @brief An interface on a side of a patch with neighbors at a finer refinement level.
	 *
	 * This interface lines up with the cells on the finer patch.
	 * The orthant value should be set to the orthant on the coarser patch's side that the finer
	 * patch lies on.
	 *
	 * @param orthant the orthant of the fine patch
	 * @return IfaceType<D> the new IfaceType
	 */
	static IfaceType<D> CoarseToFine(Orthant<D - 1> orthant)
	{
		return IfaceType<D>(4, orthant);
	}
	/**
	 * @brief Check if this type is Normal
	 *
	 * @return if it is Normal
	 */
	bool isCoarseToFine() const
	{
		return val == 4;
	}
	/**
	 * @brief Construct a new Iface Type object with the value set to 10 and the orthant set to null
	 */
	IfaceType() = default;
	/**
	 * @brief Get the Orthant that the finer patch lies on
	 */
	Orthant<D - 1> getOrthant() const
	{
		return orthant;
	}
	/**
	 * @brief Set the Orthant that the finer patch lies on
	 */
	void setOrthant(Orthant<D - 1> orthant)
	{
		this->orthant = orthant;
	}
	/**
	 * @brief Compare iface type values
	 */
	bool operator<(const IfaceType &b) const
	{
		return std::forward_as_tuple(val, orthant) < std::forward_as_tuple(b.val, b.orthant);
	}
};
} // namespace Schur
} // namespace ThunderEgg
#endif
