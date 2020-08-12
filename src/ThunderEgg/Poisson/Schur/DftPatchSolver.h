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

#ifndef THUNDEREGG_POISSON_SCHUR_DFTPATCHSOLVER_H
#define THUNDEREGG_POISSON_SCHUR_DFTPATCHSOLVER_H
#include <ThunderEgg/Schur/PatchSolver.h>
#include <ThunderEgg/Schur/SchurHelper.h>
#include <ThunderEgg/ValVector.h>
#include <bitset>
#include <map>
#include <valarray>

extern "C" void dgemv_(char &, int &, int &, double &, double *, int &, double *, int &, double &,
                       double *, int &);

namespace ThunderEgg
{
namespace Poisson
{
namespace Schur
{
enum class DftType { DCT_II, DCT_III, DCT_IV, DST_II, DST_III, DST_IV };
#ifndef DOMAINK
#define DOMAINK
template <size_t D> struct DomainK {
	unsigned long neumann = 0;
	double        h_x     = 0;

	DomainK() {}
	DomainK(const ThunderEgg::Schur::SchurInfo<D> &sinfo)
	{
		this->neumann = sinfo.pinfo->neumann.to_ulong();
		this->h_x     = sinfo.pinfo->spacings[0];
	}
	friend bool operator<(const DomainK &l, const DomainK &r)
	{
		return std::tie(l.neumann, l.h_x) < std::tie(r.neumann, r.h_x);
	}
};
#endif

template <size_t D> class DftPatchSolver : public ThunderEgg::Schur::PatchSolver<D>
{
	private:
	std::shared_ptr<ThunderEgg::Schur::SchurHelper<D>>                          sh;
	std::array<int, D>                                                          ns;
	int                                                                         n;
	int                                                                         patch_stride;
	bool                                                                        initialized = false;
	static bool                                                                 compareDomains();
	double                                                                      lambda;
	std::map<DomainK<D>, std::array<std::shared_ptr<std::valarray<double>>, D>> plan1;
	std::map<DomainK<D>, std::array<std::shared_ptr<std::valarray<double>>, D>> plan2;
	std::shared_ptr<ValVector<D>>                                               f_copy;
	std::shared_ptr<ValVector<D>>                                               tmp;
	std::shared_ptr<ValVector<D>>                                               local_tmp;
	std::map<DomainK<D>, std::valarray<double>>                                 eigen_vals;
	std::array<std::shared_ptr<std::valarray<double>>, 6>                       transforms
	= {{nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}};

	std::array<std::shared_ptr<std::valarray<double>>, D> plan(DftType *types);

	std::shared_ptr<std::valarray<double>> getTransformArray(DftType type);
	void execute_plan(std::array<std::shared_ptr<std::valarray<double>>, D> plan, LocalData<D> in,
	                  LocalData<D> out, const bool inverse);

	public:
	DftPatchSolver(std::shared_ptr<ThunderEgg::Schur::SchurHelper<D>> sh, double lambda = 0);
	void solve(ThunderEgg::Schur::SchurInfo<D> &sinfo, std::shared_ptr<const Vector<D>> f,
	           std::shared_ptr<Vector<D>> u, std::shared_ptr<const Vector<D - 1>> gamma);
	void solve(std::shared_ptr<const Vector<D>> f, std::shared_ptr<Vector<D>> u,
	           std::shared_ptr<const Vector<D - 1>> gamma) override
	{
		for (auto &sinfo : sh->getSchurInfoVector()) {
			solve(*sinfo, f, u, gamma);
		}
	}
	void addPatch(ThunderEgg::Schur::SchurInfo<D> &sinfo) override;
};

template <size_t D>
inline DftPatchSolver<D>::DftPatchSolver(std::shared_ptr<ThunderEgg::Schur::SchurHelper<D>> sh,
                                         double                                             lambda)
{
	this->sh     = sh;
	ns           = sh->getDomain()->getNs();
	n            = ns[0];
	patch_stride = sh->getDomain()->getNumCellsInPatch();
	this->lambda = lambda;
}
template <size_t D> inline void DftPatchSolver<D>::addPatch(ThunderEgg::Schur::SchurInfo<D> &sinfo)
{
	using namespace std;
	if (!initialized) {
		initialized = true;
		std::array<int, D> lengths;
		lengths.fill(n);
		f_copy = std::make_shared<ValVector<D>>(MPI_COMM_SELF, lengths, 0, 1);
		tmp    = std::make_shared<ValVector<D>>(MPI_COMM_SELF, lengths, 0, 1);
		if (!(D % 2)) {
			local_tmp = std::make_shared<ValVector<D>>(MPI_COMM_SELF, lengths, 0, 1);
		}
	}

	DftType transforms[D];
	DftType transforms_inv[D];
	if (!plan1.count(sinfo)) {
		for (size_t i = 0; i < D; i++) {
			// x direction
			if (sinfo.pinfo->isNeumann(Side<D>(2 * i))
			    && sinfo.pinfo->isNeumann(Side<D>(2 * i + 1))) {
				transforms[i]     = DftType::DCT_II;
				transforms_inv[i] = DftType::DCT_III;
			} else if (sinfo.pinfo->isNeumann(Side<D>(2 * i))) {
				transforms[i]     = DftType::DCT_IV;
				transforms_inv[i] = DftType::DCT_IV;
			} else if (sinfo.pinfo->isNeumann(Side<D>(2 * i + 1))) {
				transforms[i]     = DftType::DST_IV;
				transforms_inv[i] = DftType::DST_IV;
			} else {
				transforms[i]     = DftType::DST_II;
				transforms_inv[i] = DftType::DST_III;
			}
		}

		plan1[sinfo] = plan(transforms);
		plan2[sinfo] = plan(transforms_inv);
	}

	if (!eigen_vals.count(sinfo)) {
		valarray<double> &denom = eigen_vals[sinfo];
		denom.resize(patch_stride);

		valarray<double> ones(pow(n, D - 1));
		ones = 1;

		for (size_t i = 0; i < D; i++) {
			valarray<size_t> sizes(D - 1);
			sizes = n;
			valarray<size_t> strides(D - 1);
			for (size_t axis = 1; axis < D; axis++) {
				strides[axis - 1] = pow(n, (i + axis) % D);
			}
			double h = sinfo.pinfo->spacings[0];

			if (sinfo.pinfo->isNeumann(Side<D>(i * 2))
			    && sinfo.pinfo->isNeumann(Side<D>(i * 2 + 1))) {
				for (int xi = 0; xi < n; xi++) {
					denom[gslice(xi * pow(n, i), sizes, strides)]
					-= 4 / (h * h) * pow(sin(xi * M_PI / (2 * n)), 2) * ones;
				}
			} else if (sinfo.pinfo->isNeumann(Side<D>(i * 2))
			           || sinfo.pinfo->isNeumann(Side<D>(i * 2 + 1))) {
				for (int xi = 0; xi < n; xi++) {
					denom[gslice(xi * pow(n, i), sizes, strides)]
					-= 4 / (h * h) * pow(sin((xi + 0.5) * M_PI / (2 * n)), 2) * ones;
				}
			} else {
				for (int xi = 0; xi < n; xi++) {
					denom[gslice(xi * pow(n, i), sizes, strides)]
					-= 4 / (h * h) * pow(sin((xi + 1) * M_PI / (2 * n)), 2) * ones;
				}
			}
		}

		denom += lambda;
	}
}

template <size_t D>
inline void
DftPatchSolver<D>::solve(ThunderEgg::Schur::SchurInfo<D> &sinfo, std::shared_ptr<const Vector<D>> f,
                         std::shared_ptr<Vector<D>> u, std::shared_ptr<const Vector<D - 1>> gamma)
{
	using namespace std;
	const LocalData<D> f_view      = f->getLocalData(sinfo.pinfo->local_index);
	LocalData<D>       f_copy_view = f_copy->getLocalData(0);
	LocalData<D>       tmp_view    = tmp->getLocalData(0);

	std::array<int, D> start, end;
	start.fill(0);
	end.fill(n - 1);

	nested_loop<D>(start, end,
	               [&](std::array<int, D> coord) { f_copy_view[coord] = f_view[coord]; });

	for (Side<D> s : Side<D>::getValues()) {
		std::array<int, D - 1> start, end;
		start.fill(0);
		end.fill(n - 1);

		if (sinfo.pinfo->hasNbr(s)) {
			const LocalData<D - 1> gamma_view = gamma->getLocalData(sinfo.getIfaceLocalIndex(s));
			LocalData<D - 1>       slice      = f_copy->getLocalData(0).getSliceOnSide(s);
			double                 h2         = pow(sinfo.pinfo->spacings[s.getAxisIndex()], 2);
			nested_loop<D - 1>(start, end, [&](std::array<int, D - 1> coord) {
				slice[coord] -= 2.0 / h2 * gamma_view[coord];
			});
		}
	}

	execute_plan(plan1[sinfo], f_copy_view, tmp_view, false);

	tmp->vec /= eigen_vals[sinfo];

	if (sinfo.pinfo->neumann.all()) {
		tmp->vec[0] = 0;
	}

	LocalData<D> u_view = u->getLocalData(sinfo.pinfo->local_index);

	execute_plan(plan2[sinfo], tmp_view, u_view, false);

	double scale = pow(2.0 / n, D);
	nested_loop<D>(start, end, [&](std::array<int, D> coord) { u_view[coord] *= scale; });
}
template <size_t D>
inline std::array<std::shared_ptr<std::valarray<double>>, D> DftPatchSolver<D>::plan(DftType *types)
{
	std::array<std::shared_ptr<std::valarray<double>>, D> retval;
	for (size_t i = 0; i < D; i++) {
		retval[i] = getTransformArray(types[i]);
	}
	return retval;
}
template <size_t D>
inline std::shared_ptr<std::valarray<double>> DftPatchSolver<D>::getTransformArray(DftType type)
{
	using namespace std;
	int idx = static_cast<int>(type);

	shared_ptr<valarray<double>> matrix_ptr;

	if (transforms[idx] == nullptr) {
		matrix_ptr.reset(new valarray<double>(n * n));
		valarray<double> &matrix = *matrix_ptr;
		switch (type) {
			case DftType::DCT_II:
				for (int j = 0; j < n; j++) {
					for (int i = 0; i < n; i++) {
						matrix[i * n + j] = cos(M_PI / n * (i * (j + 0.5)));
					}
				}
				break;
			case DftType::DCT_III:
				for (int i = 0; i < n; i++) {
					matrix[i * n] = 0.5;
				}
				for (int j = 1; j < n; j++) {
					for (int i = 0; i < n; i++) {
						matrix[i * n + j] = cos(M_PI / n * ((i + 0.5) * j));
					}
				}
				break;
			case DftType::DCT_IV:
				for (int j = 0; j < n; j++) {
					for (int i = 0; i < n; i++) {
						matrix[i * n + j] = cos(M_PI / n * ((i + 0.5) * (j + 0.5)));
					}
				}
				break;
			case DftType::DST_II:
				for (int i = 0; i < n; i++) {
					for (int j = 0; j < n; j++) {
						matrix[i * n + j] = sin(M_PI / n * ((i + 1) * (j + 0.5)));
					}
				}
				break;
			case DftType::DST_III:
				for (int i = 0; i < n; i += 2) {
					matrix[i * n + n - 1] = 0.5;
				}
				for (int i = 1; i < n; i += 2) {
					matrix[i * n + n - 1] = -0.5;
				}
				for (int i = 0; i < n; i++) {
					for (int j = 0; j < n - 1; j++) {
						matrix[i * n + j] = sin(M_PI / n * ((i + 0.5) * (j + 1)));
					}
				}
				break;
			case DftType::DST_IV:
				for (int j = 0; j < n; j++) {
					for (int i = 0; i < n; i++) {
						matrix[i * n + j] = sin(M_PI / n * ((i + 0.5) * (j + 0.5)));
					}
				}
				break;
		}
	} else {
		matrix_ptr = transforms[idx];
	}
	return matrix_ptr;
}
template <size_t D>
inline void
DftPatchSolver<D>::execute_plan(std::array<std::shared_ptr<std::valarray<double>>, D> plan,
                                LocalData<D> in, LocalData<D> out, const bool inverse)
{
	LocalData<D> prev_result = in;

	std::array<int, D> start;
	start.fill(0);
	std::array<int, D> end = in.getLengths();
	for (size_t i = 0; i < D; i++) {
		end[i]--;
	}

	for (size_t dim = 0; dim < D; dim++) {
		int old_end                   = end[dim];
		end[dim]                      = 0;
		std::valarray<double> &matrix = *plan[dim];

		LocalData<D> new_result;
		if (D % 2) {
			if (dim % 2) {
				new_result = in;
			} else {
				new_result = out;
			}
		} else {
			if (dim == D - 1) {
				new_result = out;
			} else if (dim == D - 2) {
				new_result = local_tmp->getLocalData(0);
			} else if (dim % 2) {
				new_result = in;
			} else {
				new_result = out;
			}
		}

		int pstride = prev_result.getStrides()[dim];
		int nstride = new_result.getStrides()[dim];

		char   T    = 'T';
		double one  = 1;
		double zero = 0;
		nested_loop<D>(start, end, [&](std::array<int, D> coord) {
			dgemv_(T, n, n, one, &matrix[0], n, &prev_result[coord], pstride, zero,
			       &new_result[coord], nstride);
		});

		prev_result = new_result;
		end[dim]    = old_end;
	}
}
extern template class DftPatchSolver<2>;
extern template class DftPatchSolver<3>;
} // namespace Schur
} // namespace Poisson
} // namespace ThunderEgg
#endif
