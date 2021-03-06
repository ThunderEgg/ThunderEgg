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

#include "PolyChebPrec.h"
#include <ThunderEgg/ValVector.h>
#include <iostream>
using namespace std;
using namespace ThunderEgg;
using namespace ThunderEgg::Experimental;
PolyChebPrec::PolyChebPrec(std::shared_ptr<Domain<3>>          domain,
                           std::shared_ptr<InterfaceDomain<3>> sh)
{
	this->sh     = sh;
	this->domain = domain;
}
void PolyChebPrec::apply(std::shared_ptr<const Vector<2>> x, std::shared_ptr<Vector<2>> b) const
{
	std::shared_ptr<Vector<3>> f = ValVector<3>::GetNewVector(domain);
	std::shared_ptr<Vector<3>> u = ValVector<3>::GetNewVector(domain);

	// std::shared_ptr<Vector<2>> bk  = sh->getNewGlobalInterfaceVector();
	// std::shared_ptr<Vector<2>> bk1 = sh->getNewGlobalInterfaceVector();
	// std::shared_ptr<Vector<2>> bk2 = sh->getNewGlobalInterfaceVector();

	for (int i = coeffs.size() - 1; i > 0; i--) {
		// solver->solve(f, u, bk1);
		// interp->interpolateToInterface(u, bk);
		// bk->scaleThenAddScaled(4 / interval, -2, bk1);
		// bk->addScaled(coeffs[i], x, -1, bk2);
		// auto tmp = bk2;
		// bk2      = bk1;
		// bk1      = bk;
		// bk       = tmp;
	}
	// solver->solve(f, u, bk1);
	// interp->interpolateToInterface(u, b);
	// b->scaleThenAddScaled(2 / interval, -1, bk1);
	// b->addScaled(coeffs[0], x, -1, bk2);
}
