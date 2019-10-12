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

#ifndef VTKWRITER2D_H
#define VTKWRITER2D_H
#include <Thunderegg/Domain.h>
#include <Thunderegg/Vector.h>
#include <map>
#include <petscvec.h>
#include <set>
#include <string>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkImageData.h>
#include <vtkMPIController.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkMultiPieceDataSet.h>
#include <vtkSmartPointer.h>
#include <vtkXMLMultiBlockDataWriter.h>
#include <vtkXMLPMultiBlockDataWriter.h>
class VtkWriter2d
{
	private:
	std::shared_ptr<Thunderegg::Domain<2>>            dc;
	std::string                                       file_name;
	static vtkSmartPointer<vtkMultiProcessController> controller;
	std::map<int, vtkSmartPointer<vtkImageData>>      images;
	std::set<vtkSmartPointer<vtkDoubleArray>>         arrays;
	vtkSmartPointer<vtkXMLPMultiBlockDataWriter>      writer;
	vtkSmartPointer<vtkMultiBlockDataSet>             block;
	vtkSmartPointer<vtkMultiPieceDataSet>             data;

	public:
	VtkWriter2d(std::shared_ptr<Thunderegg::Domain<2>> dc, std::string file_name);
	void add(std::shared_ptr<const Thunderegg::Vector<2>> u, std::string name);
	void write();
};
#endif
