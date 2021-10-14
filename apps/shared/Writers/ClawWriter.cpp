/***************************************************************************
 *  ThunderEgg, a library for solvers on adaptively refined block-structured 
 *  Cartesian grids.
 *
 *  Copyright (c) 2017-2021 Scott Aiton
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

#include "ClawWriter.h"
#include <fstream>
using namespace std;
using namespace ThunderEgg;
ClawWriter::ClawWriter(const Domain<2> &domain)
: domain(domain)
{
}
void ClawWriter::addVector(const Vector<2> &vec)
{
	vectors.push_back(vec);
}
void ClawWriter::write()
{
	ofstream     t_file("fort.t0000");
	const string tab = "\t";
	t_file << 0.0 << tab << "time" << endl;
	t_file << vectors.size() << tab << "meqn" << endl;
	t_file << domain.getNumLocalPatches() << tab << "ngrids" << endl;
	t_file << 2 << tab << "num_aux" << endl;
	t_file << 2 << tab << "num_dim" << endl;
	t_file.close();
	ofstream q_file("fort.q0000");

	q_file.precision(10);
	q_file << scientific;
	for (auto &pinfo : domain.getPatchInfoVector()) {
		writePatch(pinfo, q_file);
	}
	q_file.close();
}
void ClawWriter::writePatch(const PatchInfo<2> &pinfo, std::ostream &os)
{
	const string tab = "\t";
	os << pinfo.id << tab << "grid_number" << endl;
	os << pinfo.refine_level << tab << "AMR_level" << endl;
	os << 0 << tab << "block_number" << endl;
	os << 0 << tab << "mpi_rank" << endl;
	os << pinfo.ns[0] << tab << "mx" << endl;
	os << pinfo.ns[1] << tab << "my" << endl;
	os << pinfo.starts[0] << tab << "xlow" << endl;
	os << pinfo.starts[1] << tab << "ylow" << endl;
	os << pinfo.spacings[0] << tab << "dx" << endl;
	os << pinfo.spacings[1] << tab << "dy" << endl;
	os << endl;
	list<ComponentView<const double, 2>> lds;
	for (const auto &vec : vectors) {
		lds.push_back(vec.getComponentView(0, pinfo.local_index));
	}
	for (int y = 0; y < pinfo.ns[1]; y++) {
		for (int x = 0; x < pinfo.ns[0]; x++) {
			auto lds_iter       = lds.begin();
			auto one_before_end = --lds.end();
			while (lds_iter != one_before_end) {
				os << (*lds_iter)[{x, y}] << tab;
				lds_iter++;
			}
			os << (*lds_iter)[{x, y}] << endl;
		}
		os << endl;
	}
}
