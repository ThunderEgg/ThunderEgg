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

#include "VtkWriter.h"
using namespace std;
using namespace ThunderEgg;
vtkSmartPointer<vtkMultiProcessController> VtkWriter::controller;
VtkWriter::VtkWriter(shared_ptr<Domain<3>> dc, string file_name)
{
	if (controller == nullptr) {
		controller = vtkSmartPointer<vtkMPIController>::New();
		controller->Initialize(0, 0, 1);
	}
	this->dc        = dc;
	this->file_name = file_name;
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	writer = vtkSmartPointer<vtkXMLPMultiBlockDataWriter>::New();
	writer->SetController(controller);
	block = vtkSmartPointer<vtkMultiBlockDataSet>::New();
	data  = vtkSmartPointer<vtkMultiPieceDataSet>::New();
	data->SetNumberOfPieces(dc->getNumGlobalCells());
	block->SetNumberOfBlocks(1);
	block->SetBlock(0, data);

	string out = file_name + ".vtmb";
	writer->SetFileName(out.c_str());
	writer->SetNumberOfPieces(1);
	writer->SetStartPiece(0);
	writer->SetInputData(block);
	if (rank == 0) { writer->SetWriteMetaFile(1); }
}
void VtkWriter::add(Vec u, string name)
{
	// create MultiPieceDataSet and fill with patch information

	double *u_view;
	VecGetArray(u, &u_view);

	for (auto &pinfo : dc->getPatchInfoVector()) {
		PatchInfo<3> &d     = *pinfo;
		double        h_x   = d.spacings[0];
		double        h_y   = d.spacings[1];
		double        h_z   = d.spacings[2];
		int           start = d.local_index * d.ns[0] * d.ns[1] * d.ns[2];

		// create image object
		vtkSmartPointer<vtkImageData> image = images[d.id];
		if (image == nullptr) {
			image        = vtkSmartPointer<vtkImageData>::New();
			images[d.id] = image;
			image->SetOrigin(d.starts[0], d.starts[1], d.starts[2]);
			image->SetSpacing(h_x, h_y, h_z);
			image->SetExtent(d.starts[0], d.starts[0] + d.spacings[0] * d.ns[0], d.starts[1],
			                 d.starts[1] + d.ns[1] * d.spacings[1], d.starts[2],
			                 d.starts[2] + d.spacings[2] * d.ns[2]);
			image->PrepareForNewData();
			image->SetDimensions(d.ns[0] + 1, d.ns[1] + 1, d.ns[2] + 1);
			// add image to dataset
			data->SetPiece(d.global_index, image);
		}

		// create solution vector
		vtkSmartPointer<vtkDoubleArray> solution = vtkSmartPointer<vtkDoubleArray>::New();
		solution->SetNumberOfComponents(1);
		solution->SetNumberOfValues(dc->getNumCellsInPatch());
		solution->SetName(name.c_str());
		for (int i = 0; i < dc->getNumCellsInPatch(); i++) {
			solution->SetValue(i, u_view[start + i]);
		}
		image->GetCellData()->AddArray(solution);
	}

	VecRestoreArray(u, &u_view);
}
void VtkWriter::write()
{
	// write data
	writer->Update();
	writer->Write();
}
