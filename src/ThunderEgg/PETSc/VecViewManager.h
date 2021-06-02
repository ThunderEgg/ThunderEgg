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

#ifndef THUNDEREGG_PETSC_VECVIEWMANAGER_H
#define THUNDEREGG_PETSC_VECVIEWMANAGER_H
#include <ThunderEgg/RuntimeError.h>
#include <ThunderEgg/Vector.h>
#include <petscvec.h>

namespace ThunderEgg
{
namespace PETSc
{
/**
 * @brief Manages access to a PETSc Vec pointer to data
 *
 * This calls VecGetArray() and VecRestoreArray() for PETSc vectors
 */
class VecViewManager : public ViewManager
{
	private:
	/**
	 * @brief the underlying PETSc Vec
	 */
	Vec vec;
	/**
	 * @brief is the Vec read only
	 */
	bool is_read_only;
	/**
	 * @brief pointer to the vector view
	 */
	double *vec_view;

	public:
	/**
	 * @brief Construct a new VecViewManager object
	 *
	 * @param vec_in the PETSc Vec to manage
	 * @param is_read_only_in whether or not to treat this vector as read only
	 */
	VecViewManager(Vec vec_in, bool is_read_only_in) : vec(vec_in), is_read_only(is_read_only_in)
	{
		if (is_read_only) {
			VecGetArrayRead(vec, const_cast<const double **>(&vec_view));
		} else {
			int state;
			VecLockGet(vec, &state);
			if (state != 0) {
				throw RuntimeError("PETSc Vec is in state " + std::to_string(state));
			}
			VecGetArray(vec, &vec_view);
		}
	}
	VecViewManager(const VecViewManager &) = delete;
	VecViewManager &operator=(const VecViewManager &) = delete;
	VecViewManager(VecViewManager &&) noexcept        = delete;
	VecViewManager &operator=(VecViewManager &&) noexcept = delete;
	/**
	 * @brief Destroy the VecViewManager object
	 */
	~VecViewManager()
	{
		if (is_read_only) {
			VecRestoreArrayRead(vec, const_cast<const double **>(&vec_view));
		} else {
			VecRestoreArray(vec, &vec_view);
		}
	}
	/**
	 * @brief Get the pointer to the PETSc vector view
	 *
	 * @return double* the view
	 */
	double *getVecView() const
	{
		return vec_view;
	}
};
} // namespace PETSc
} // namespace ThunderEgg
#endif
