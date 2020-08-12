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

#ifndef THUNDEREGG_POISSON_SCHUR_FFTWPATCHSOLVER_H
#define THUNDEREGG_POISSON_SCHUR_FFTWPATCHSOLVER_H
#include <ThunderEgg/Poisson/Schur/StarPatchOperator.h>
#include <ThunderEgg/Schur/PatchSolver.h>
#include <ThunderEgg/Schur/SchurHelper.h>
#include <ThunderEgg/ValVector.h>
#include <bitset>
#include <fftw3.h>
#include <map>
namespace ThunderEgg
{
namespace Poisson
{
namespace Schur
{
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

template <size_t D> class FftwPatchSolver : public ThunderEgg::Schur::PatchSolver<D>
{
	private:
	int                                                n;
	bool                                               initialized = false;
	static bool                                        compareDomains();
	double                                             lambda;
	std::map<DomainK<D>, fftw_plan>                    plan1;
	std::map<DomainK<D>, fftw_plan>                    plan2;
	std::shared_ptr<ValVector<D>>                      f_copy;
	std::shared_ptr<ValVector<D>>                      tmp;
	std::shared_ptr<ValVector<D>>                      sol;
	std::array<int, D + 1>                             npow;
	std::map<DomainK<D>, std::valarray<double>>        denoms;
	std::shared_ptr<ThunderEgg::Schur::SchurHelper<D>> sh;

	public:
	FftwPatchSolver(std::shared_ptr<ThunderEgg::Schur::SchurHelper<D>> sh, double lambda = 0);
	~FftwPatchSolver();
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
FftwPatchSolver<D>::FftwPatchSolver(std::shared_ptr<ThunderEgg::Schur::SchurHelper<D>> sh,
                                    double                                             lambda)
{
	this->sh     = sh;
	n            = sh->getLengths()[0];
	this->lambda = lambda;
	for (auto &sinfo : sh->getSchurInfoVector()) {
		addPatch(*sinfo);
	}
}
template <size_t D> FftwPatchSolver<D>::~FftwPatchSolver()
{
	for (auto p : plan1) {
		fftw_destroy_plan(p.second);
	}
	for (auto p : plan2) {
		fftw_destroy_plan(p.second);
	}
}
template <size_t D> void FftwPatchSolver<D>::addPatch(ThunderEgg::Schur::SchurInfo<D> &sinfo)
{
	using namespace std;
	if (!initialized) {
		initialized = true;
		for (size_t i = 0; i <= D; i++) {
			npow[i] = (int) std::pow(n, i);
		}
		std::array<int, D> lengths;
		lengths.fill(n);
		f_copy = std::make_shared<ValVector<D>>(MPI_COMM_SELF, lengths, 0, 1);
		tmp    = std::make_shared<ValVector<D>>(MPI_COMM_SELF, lengths, 0, 1);
		sol    = std::make_shared<ValVector<D>>(MPI_COMM_SELF, lengths, 0, 1);
	}

	int           ns[D];
	fftw_r2r_kind transforms[D];
	fftw_r2r_kind transforms_inv[D];
	if (!plan1.count(sinfo)) {
		for (size_t i = 0; i < D; i++) {
			ns[D - 1 - i] = n;
			// x direction
			if (sinfo.pinfo->isNeumann(Side<D>(2 * i))
			    && sinfo.pinfo->isNeumann(Side<D>(2 * i + 1))) {
				transforms[D - 1 - i]     = FFTW_REDFT10;
				transforms_inv[D - 1 - i] = FFTW_REDFT01;
			} else if (sinfo.pinfo->isNeumann(Side<D>(2 * i))) {
				transforms[D - 1 - i]     = FFTW_REDFT11;
				transforms_inv[D - 1 - i] = FFTW_REDFT11;
			} else if (sinfo.pinfo->isNeumann(Side<D>(2 * i + 1))) {
				transforms[D - 1 - i]     = FFTW_RODFT11;
				transforms_inv[D - 1 - i] = FFTW_RODFT11;
			} else {
				transforms[D - 1 - i]     = FFTW_RODFT10;
				transforms_inv[D - 1 - i] = FFTW_RODFT01;
			}
		}

		plan1[sinfo] = fftw_plan_r2r(D, ns, &f_copy->vec[0], &tmp->vec[0], transforms,
		                             FFTW_MEASURE | FFTW_DESTROY_INPUT);
		plan2[sinfo] = fftw_plan_r2r(D, ns, &tmp->vec[0], &sol->vec[0], transforms_inv,
		                             FFTW_MEASURE | FFTW_DESTROY_INPUT);
	}

	if (!denoms.count(sinfo)) {
		valarray<double> &denom = denoms[sinfo];
		denom.resize(pow(n, D));

		valarray<double> ones(pow(n, D - 1));
		ones = 1;

		for (size_t i = 0; i < D; i++) {
			valarray<size_t> sizes(D - 1);
			sizes = n;
			valarray<size_t> strides(D - 1);
			for (size_t sinfo = 1; sinfo < D; sinfo++) {
				strides[sinfo - 1] = pow(n, (i + sinfo) % D);
			}
			double h = sinfo.pinfo->spacings[i];

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
void FftwPatchSolver<D>::solve(ThunderEgg::Schur::SchurInfo<D> &sinfo,
                               std::shared_ptr<const Vector<D>> f, std::shared_ptr<Vector<D>> u,
                               std::shared_ptr<const Vector<D - 1>> gamma)
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

	StarPatchOperator<D> op(sh);
	op.addInterfaceToRHS(sinfo, gamma, f_copy->getLocalData(0));

	fftw_execute(plan1[sinfo]);

	tmp->vec /= denoms[sinfo];

	if (sinfo.pinfo->neumann.all()) {
		tmp->vec[0] = 0;
	}

	fftw_execute(plan2[sinfo]);

	sol->vec /= pow(2.0 * n, D);

	LocalData<D> u_view   = u->getLocalData(sinfo.pinfo->local_index);
	LocalData<D> sol_view = sol->getLocalData(0);
	nested_loop<D>(start, end, [&](std::array<int, D> coord) { u_view[coord] = sol_view[coord]; });
}
extern template class FftwPatchSolver<2>;
extern template class FftwPatchSolver<3>;
} // namespace Schur
} // namespace Poisson
} // namespace ThunderEgg
#endif