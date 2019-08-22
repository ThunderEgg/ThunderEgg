/***************************************************************************
 *  Thunderegg, a library for solving Poisson's equation on adaptively
 *  refined block-structured Cartesian grids
 *
 *  Copyright (C) 2019  Thunderegg Developers. See AUTHORS.md file at the
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

#ifndef THUNDEREGG_VARPOISSON_STARPATCHOPERATOR_H
#define THUNDEREGG_VARPOISSON_STARPATCHOPERATOR_H

#include <Thunderegg/DomainTools.h>
#include <Thunderegg/PatchOperator.h>
#include <Thunderegg/SchurHelper.h>

namespace Thunderegg
{
namespace VarPoisson
{
template <size_t D> class StarPatchOperator : public PatchOperator<D>
{
	private:
	std::shared_ptr<const Vector<D>>     coeffs;
	std::shared_ptr<const Vector<D - 1>> bc_coeffs;
	std::shared_ptr<const Vector<D - 1>> gamma_coeffs;

	public:
	StarPatchOperator(std::shared_ptr<const Vector<D>> h, std::shared_ptr<const Vector<D - 1>> h_bc,
	                  std::shared_ptr<SchurHelper<D>> sh)
	{
		coeffs    = h;
		bc_coeffs = h_bc;
	}
	void applyWithInterface(SchurInfo<D> &sinfo, const LocalData<D> u,
	                        std::shared_ptr<const Vector<D - 1>> gamma, LocalData<D> f) override
	{
		const LocalData<D>    c  = coeffs->getLocalData(sinfo.pinfo->local_index);
		std::array<double, D> h2 = sinfo.pinfo->spacings;
		for (size_t i = 0; i < D; i++) {
			h2[i] *= h2[i];
		}
		loop<0, D - 1>([&](int axis) {
			Side<D> lower_side = axis * 2;
			Side<D> upper_side = axis * 2 + 1;
			if (sinfo.pinfo->hasNbr(lower_side)) {
				LocalData<D - 1>       f_slice = f.getSliceOnSide(lower_side);
				const LocalData<D - 1> bnd
				= gamma->getLocalData(sinfo.getIfaceLocalIndex(lower_side));
				const LocalData<D - 1> mid   = u.getSliceOnSide(lower_side);
				const LocalData<D - 1> upper = u.getSliceOnSide(lower_side, 1);
				const LocalData<D - 1> c_bnd
				= gamma_coeffs->getLocalData(sinfo.getIfaceLocalIndex(lower_side));
				const LocalData<D - 1> c_mid   = c.getSliceOnSide(lower_side);
				const LocalData<D - 1> c_upper = c.getSliceOnSide(lower_side, 1);

				nested_loop<D - 1>(mid.getStart(), mid.getEnd(), [&](std::array<int, D - 1> coord) {
					f_slice[coord] += ((c_upper[coord] + c_mid[coord]) * (upper[coord] - mid[coord])
					                   - 4 * c_bnd[coord] * (mid[coord] - bnd[coord]))
					                  / (2 * h2[axis]);
				});
			} else if (sinfo.pinfo->isNeumann(lower_side)) {
				LocalData<D - 1>       f_slice = f.getSliceOnSide(lower_side);
				const LocalData<D - 1> mid     = u.getSliceOnSide(lower_side);
				const LocalData<D - 1> upper   = u.getSliceOnSide(lower_side, 1);

				nested_loop<D - 1>(mid.getStart(), mid.getEnd(), [&](std::array<int, D - 1> coord) {
					f_slice[coord] += (-mid[coord] + upper[coord]) / h2[axis];
				});
			} else {
				LocalData<D - 1>       f_slice = f.getSliceOnSide(lower_side);
				const LocalData<D - 1> mid     = u.getSliceOnSide(lower_side);
				const LocalData<D - 1> upper   = u.getSliceOnSide(lower_side, 1);
				const LocalData<D - 1> c_bnd
				= bc_coeffs->getLocalData(sinfo.pinfo->getBCLocalIndex(lower_side));
				const LocalData<D - 1> c_mid   = c.getSliceOnSide(lower_side);
				const LocalData<D - 1> c_upper = c.getSliceOnSide(lower_side, 1);
				nested_loop<D - 1>(mid.getStart(), mid.getEnd(), [&](std::array<int, D - 1> coord) {
					f_slice[coord]
					+= ((c_upper[coord] + c_mid[coord]) * upper[coord]
					    - (4 * c_bnd[coord] + c_upper[coord] + c_mid[coord]) * mid[coord])
					   / (2 * h2[axis]);
				});
			}
			// middle
			{
				std::array<int, D> start = f.getStart();
				std::array<int, D> end   = f.getEnd();
				start[axis] += 1;
				end[axis] -= 1;
				int stride   = u.getStrides()[axis];
				int c_stride = c.getStrides()[axis];
				nested_loop<D>(start, end, [&](std::array<int, D> coord) {
					const double *ptr     = u.getPtr(coord);
					const double *c_ptr   = c.getPtr(coord);
					double        lower   = *(ptr - stride);
					double        mid     = *ptr;
					double        upper   = *(ptr + stride);
					double        c_lower = *(c_ptr - c_stride);
					double        c_mid   = *c_ptr;
					double        c_upper = *(c_ptr + c_stride);
					f[coord]
					+= ((c_upper + c_mid) * (upper - mid) - (c_lower + c_mid) * (mid - lower))
					   / (2 * h2[axis]);
				});
			}
			// east
			if (sinfo.pinfo->hasNbr(upper_side)) {
				LocalData<D - 1>       f_slice = f.getSliceOnSide(upper_side);
				const LocalData<D - 1> lower   = u.getSliceOnSide(upper_side, 1);
				const LocalData<D - 1> mid     = u.getSliceOnSide(upper_side);
				const LocalData<D - 1> bnd
				= gamma->getLocalData(sinfo.getIfaceLocalIndex(upper_side));

				nested_loop<D - 1>(mid.getStart(), mid.getEnd(), [&](std::array<int, D - 1> coord) {
					f_slice[coord] += (lower[coord] - 3 * mid[coord] + 2 * bnd[coord]) / h2[axis];
				});
			} else if (sinfo.pinfo->isNeumann(upper_side)) {
				LocalData<D - 1>       f_slice = f.getSliceOnSide(upper_side);
				const LocalData<D - 1> lower   = u.getSliceOnSide(upper_side, 1);
				const LocalData<D - 1> mid     = u.getSliceOnSide(upper_side);

				nested_loop<D - 1>(mid.getStart(), mid.getEnd(), [&](std::array<int, D - 1> coord) {
					f_slice[coord] += (lower[coord] - mid[coord]) / h2[axis];
				});
			} else {
				LocalData<D - 1>       f_slice = f.getSliceOnSide(upper_side);
				const LocalData<D - 1> lower   = u.getSliceOnSide(upper_side, 1);
				const LocalData<D - 1> mid     = u.getSliceOnSide(upper_side);
				const LocalData<D - 1> c_lower = c.getSliceOnSide(upper_side, 1);
				const LocalData<D - 1> c_mid   = c.getSliceOnSide(upper_side);
				const LocalData<D - 1> c_bnd
				= bc_coeffs->getLocalData(sinfo.pinfo->getBCLocalIndex(upper_side));

				nested_loop<D - 1>(mid.getStart(), mid.getEnd(), [&](std::array<int, D - 1> coord) {
					f_slice[coord]
					+= ((c_lower[coord] + c_mid[coord]) * lower[coord]
					    - (4 * c_bnd[coord] + c_lower[coord] + c_mid[coord]) * mid[coord])
					   / (2 * h2[axis]);
				});
			}
		});
	}
	void addInterfaceToRHS(SchurInfo<D> &sinfo, std::shared_ptr<const Vector<D - 1>> gamma,
	                       LocalData<D> f) override
	{
		throw 3;
		for (Side<D> s : Side<D>::getValues()) {
			if (sinfo.pinfo->hasNbr(s)) {
				const LocalData<D - 1> gamma_view
				= gamma->getLocalData(sinfo.getIfaceLocalIndex(s));

				LocalData<D - 1> slice = f.getSliceOnSide(s);

				double h2 = pow(sinfo.pinfo->spacings[s.axis()], 2);

				nested_loop<D - 1>(gamma_view.getStart(), gamma_view.getEnd(),
				                   [&](std::array<int, D - 1> coord) {
					                   slice[coord] -= 2.0 / h2 * gamma_view[coord];
				                   });
			}
		}
	}
	void apply(const SchurInfo<D> &sinfo, const LocalData<D> u, LocalData<D> f) override
	{
		throw 3;
		std::array<double, D> h2 = sinfo.pinfo->spacings;
		for (size_t i = 0; i < D; i++) {
			h2[i] *= h2[i];
		}
		{
			constexpr int axis       = 0;
			Side<D>       lower_side = axis * 2;
			Side<D>       upper_side = axis * 2 + 1;
			if (sinfo.pinfo->isNeumann(lower_side)) {
				LocalData<D - 1>       f_slice = f.getSliceOnSide(lower_side);
				const LocalData<D - 1> mid     = u.getSliceOnSide(lower_side);
				const LocalData<D - 1> upper   = u.getSliceOnSide(lower_side, 1);

				nested_loop<D - 1>(mid.getStart(), mid.getEnd(), [&](std::array<int, D - 1> coord) {
					f_slice[coord] = (-mid[coord] + upper[coord]) / h2[axis];
				});
			} else {
				LocalData<D - 1>       f_slice = f.getSliceOnSide(lower_side);
				const LocalData<D - 1> mid     = u.getSliceOnSide(lower_side);
				const LocalData<D - 1> upper   = u.getSliceOnSide(lower_side, 1);

				nested_loop<D - 1>(mid.getStart(), mid.getEnd(), [&](std::array<int, D - 1> coord) {
					f_slice[coord] = (-3 * mid[coord] + upper[coord]) / h2[axis];
				});
			}
			// middle
			{
				std::array<int, D> start = f.getStart();
				std::array<int, D> end   = f.getEnd();
				start[axis] += 1;
				end[axis] -= 1;
				int stride = u.getStrides()[axis];
				nested_loop<D>(start, end, [&](std::array<int, D> coord) {
					const double *ptr   = u.getPtr(coord);
					double        lower = *(ptr - stride);
					double        mid   = *ptr;
					double        upper = *(ptr + stride);
					f[coord]            = (lower - 2 * mid + upper) / h2[axis];
				});
			}
			// east
			if (sinfo.pinfo->isNeumann(upper_side)) {
				LocalData<D - 1>       f_slice = f.getSliceOnSide(upper_side);
				const LocalData<D - 1> lower   = u.getSliceOnSide(upper_side, 1);
				const LocalData<D - 1> mid     = u.getSliceOnSide(upper_side);

				nested_loop<D - 1>(mid.getStart(), mid.getEnd(), [&](std::array<int, D - 1> coord) {
					f_slice[coord] = (lower[coord] - mid[coord]) / h2[axis];
				});
			} else {
				LocalData<D - 1>       f_slice = f.getSliceOnSide(upper_side);
				const LocalData<D - 1> lower   = u.getSliceOnSide(upper_side, 1);
				const LocalData<D - 1> mid     = u.getSliceOnSide(upper_side);

				nested_loop<D - 1>(mid.getStart(), mid.getEnd(), [&](std::array<int, D - 1> coord) {
					f_slice[coord] = (lower[coord] - 3 * mid[coord]) / h2[axis];
				});
			}
		}
		for (size_t axis = 1; axis < D; axis++) {
			Side<D> lower_side = axis * 2;
			Side<D> upper_side = axis * 2 + 1;
			if (sinfo.pinfo->isNeumann(lower_side)) {
				LocalData<D - 1>       f_slice = f.getSliceOnSide(lower_side);
				const LocalData<D - 1> mid     = u.getSliceOnSide(lower_side);
				const LocalData<D - 1> upper   = u.getSliceOnSide(lower_side, 1);

				nested_loop<D - 1>(mid.getStart(), mid.getEnd(), [&](std::array<int, D - 1> coord) {
					f_slice[coord] += (-mid[coord] + upper[coord]) / h2[axis];
				});
			} else {
				LocalData<D - 1>       f_slice = f.getSliceOnSide(lower_side);
				const LocalData<D - 1> mid     = u.getSliceOnSide(lower_side);
				const LocalData<D - 1> upper   = u.getSliceOnSide(lower_side, 1);

				nested_loop<D - 1>(mid.getStart(), mid.getEnd(), [&](std::array<int, D - 1> coord) {
					f_slice[coord] += (-3 * mid[coord] + upper[coord]) / h2[axis];
				});
			}
			// middle
			{
				std::array<int, D> start = f.getStart();
				std::array<int, D> end   = f.getEnd();
				start[axis] += 1;
				end[axis] -= 1;
				int stride = u.getStrides()[axis];
				nested_loop<D>(start, end, [&](std::array<int, D> coord) {
					const double *ptr   = u.getPtr(coord);
					double        lower = *(ptr - stride);
					double        mid   = *ptr;
					double        upper = *(ptr + stride);
					f[coord] += (lower - 2 * mid + upper) / h2[axis];
				});
			}
			// east
			if (sinfo.pinfo->isNeumann(upper_side)) {
				LocalData<D - 1>       f_slice = f.getSliceOnSide(upper_side);
				const LocalData<D - 1> lower   = u.getSliceOnSide(upper_side, 1);
				const LocalData<D - 1> mid     = u.getSliceOnSide(upper_side);

				nested_loop<D - 1>(mid.getStart(), mid.getEnd(), [&](std::array<int, D - 1> coord) {
					f_slice[coord] += (lower[coord] - mid[coord]) / h2[axis];
				});
			} else {
				LocalData<D - 1>       f_slice = f.getSliceOnSide(upper_side);
				const LocalData<D - 1> lower   = u.getSliceOnSide(upper_side, 1);
				const LocalData<D - 1> mid     = u.getSliceOnSide(upper_side);

				nested_loop<D - 1>(mid.getStart(), mid.getEnd(), [&](std::array<int, D - 1> coord) {
					f_slice[coord] += (lower[coord] - 3 * mid[coord]) / h2[axis];
				});
			}
		}
	}
	static void addDrichletBCToRHS(std::shared_ptr<Domain<D>> domain, std::shared_ptr<Vector<D>> f,
	                               std::shared_ptr<Vector<D>>                           h,
	                               std::function<double(const std::array<double, D> &)> gfunc,
	                               std::function<double(const std::array<double, D> &)> hfunc)
	{
		/*
		for (int i = 0; i < vec->getNumLocalPatches(); i++) {
		    LocalData<D> f_ld  = f->getLocalData(i);
		    LocalData<D> h_ld  = h->getLocalData(i);
		    auto         pinfo = domain->getPatchInfoVector()[i];
		    for (int axis = 0; axis < D; axis++) {
		        Side<D> lower = 2 * axis;
		        Side<D> upper = 2 * axis + 1;
		        // west
		        if (!pinfo->hasNbr(lower)) {
		            LocalData<D> f_slice = f->getLocalData(i);
		            LocalData<D> h_slice = h_ld for (int zi = 0; zi < pinfo->ns[2]; zi++)
		            nested_loop<D>(f_slice.getStart(), f_slice.getEnd(),
		                           [&](const std::array<int, D> &coord) {
		                               f_slice[0] getRealCoord(pinfo, coord, real_coord);
		                               ld[coord] = func(real_coord);
		                           });
		        }
		        // east
		        if (!pinfo->hasNbr(upper)) {
		            for (int zi = 0; zi < pinfo->ns[2]; zi++) {
		                for (int yi = 0; yi < pinfo->ns[1]; yi++) {
		                    double x, y, z;
		                    getXYZ(pinfo, pinfo->ns[0], yi, zi, x, y, z);
		                    f_vals[index(pinfo, pinfo->ns[0] - 1, yi, zi)] -= nfunx(x, y, z) / h_x;
		                }
		            }
		        }
		    }
		}
		*/
	}
};
extern template class StarPatchOperator<2>;
extern template class StarPatchOperator<3>;
} // namespace VarPoisson
} // namespace Thunderegg
#endif