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
	char val = -1;
	/**
	 * @brief Orthant that the interface lies on
	 */
	Orthant<D - 1> orthant = Orthant<D - 1>::null();

	public:
	// enum definitions
	/**
	 * @brief An interface on a side of a patch with a neighbor at the same refinement level.
	 */
	static constexpr char normal = 0b000;
	/**
	 * @brief An interface on a side of a patch with neighbors at a finer refinement level.
	 *
	 * This interface lines up with the cells on the coarser patch.
	 */
	static constexpr char coarse_to_coarse = 0b001;
	/**
	 * @brief An interface on a side of a patch with a neighbor at a coarser refinement level.
	 *
	 * This interface lines up with the cells on the coarser patch.
	 */
	static constexpr char fine_to_coarse = 0b010;
	/**
	 * @brief An interface on a side of a patch with a neighbor at a coarser refinement level.
	 *
	 * This interface lines up with the cells on the finer patch.
	 * The orthant value should be set to the orthant on the coarser patch's side that the finer
	 * patch lies on.
	 */
	static constexpr char fine_to_fine = 0b011;
	/**
	 * @brief An interface on a side of a patch with neighbors at a finer refinement level.
	 *
	 * This interface lines up with the cells on the finer patch.
	 * The orthant value should be set to the orthant on the coarser patch's side that the finer
	 * patch lies on.
	 */
	static constexpr char coarse_to_fine = 0b100;
	/**
	 * @brief Construct a new Iface Type object with the value set to -1
	 */
	IfaceType() = default;
	/**
	 * @brief Construct a new Iface Type object with the specified value
	 *
	 * @param val the value
	 */
	IfaceType(const char val)
	{
		this->val = val;
	}
	/**
	 * @brief Construct a new Iface Type object with the specified value and orthant.
	 *
	 * @param val the value
	 * @param orthant orthant on the coarser patch's side that the finer
	 * patch lies on.
	 */
	IfaceType(const char val, const Orthant<D - 1> orthant)
	{
		this->val     = val;
		this->orthant = orthant;
	}
	/**
	 * @brief Return the integer value of this ifacetype
	 */
	char toInt()
	{
		return val;
	}
	/**
	 * @brief Get the Orthant that the finer patch lies on
	 */
	Orthant<D - 1> getOrthant()
	{
		return orthant;
	}
	/**
	 * @brief Set the Orthant that the finer patch lies on
	 */
	void setOrthant(const Orthant<D - 1> orthant)
	{
		this->orthant = orthant;
	}
	/**
	 * @brief Compare iface type values
	 */
	bool operator<(const IfaceType &b) const
	{
		int orth   = orthant.getIndex();
		int b_orth = b.orthant.getIndex();
		return std::tie(val, orth) < std::tie(b.val, b_orth);
	}
};
} // namespace Schur
} // namespace ThunderEgg
#endif
