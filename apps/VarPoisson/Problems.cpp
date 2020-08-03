/***************************************************************************
 *  ThunderEgg, a library for solving Poisson's equation on adaptively
 *  refined block-structured Cartesian grids
 *
 *  Copyright (C) 2019-2020  ThunderEgg Developers. See AUTHORS.md file at
 *  the top-level directory.
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
#include "Problems.h"
#include <cmath>
#include <locale>
void GetProblem(std::string problem, function2d &u, function2d &f, function2d &beta)
{
	std::locale loc;
	problem = std::tolower(problem, loc);
	if (problem == "exp") {
		u = [](const std::array<double, 2> &coord) {
			double x = coord[0];
			double y = coord[1];
			return (1 - x) * x * (1 - y) * y * exp(x * y);
		};
		f = [](const std::array<double, 2> &coord) {
			double x = coord[0];
			double y = coord[1];
			return exp(x * y)
			       * (pow(x, 5) * (-1 + y) * pow(y, 2) + y * (-2 + (5 - 3 * y) * y)
			          - pow(x, 4) * y * (4 + (-7 + y) * y)
			          - 2 * x * (1 + 2 * (-2 + y) * (-1 + y) * pow(y, 2))
			          - pow(x, 2) * (-5 + 8 * y + (-6 + y) * (-1 + y) * pow(y, 3))
			          + pow(x, 3) * (-3 + y * (12 - 6 * y - pow(y, 3) + pow(y, 4))));
		};
		beta = [](const std::array<double, 2> &coord) {
			double x = coord[0];
			double y = coord[1];
			return 1;
		};
	} else if (problem == "trig") {
		u = [](const std::array<double, 2> &coord) {
			double x = coord[0];
			double y = coord[1];
			return sinl(M_PI * y) * cosl(2 * M_PI * x);
		};
		f = [](const std::array<double, 2> &coord) {
			double x = coord[0];
			double y = coord[1];
			return -5 * M_PI * M_PI * sinl(M_PI * y) * cosl(2 * M_PI * x);
		};
		beta = [](const std::array<double, 2> &coord) { return 1; };
	}
}