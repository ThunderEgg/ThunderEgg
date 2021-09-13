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

#ifndef CLAWWRITER_H
#define CLAWWRITER_H
#include <ThunderEgg/Domain.h>
#include <ThunderEgg/Vector.h>
#include <list>
class ClawWriter
{
	private:
	ThunderEgg::Domain<2>            domain;
	std::list<ThunderEgg::Vector<2>> vectors;
	void                             writePatch(const ThunderEgg::PatchInfo<2> &d, std::ostream &os);

	public:
	ClawWriter(const ThunderEgg::Domain<2> &domain);
	void addVector(const ThunderEgg::Vector<2> &vec);
	void write();
};
#endif
