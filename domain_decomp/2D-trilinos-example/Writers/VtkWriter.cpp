#include "VtkWriter.h"
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkImageData.h>
#include <vtkMPIController.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkMultiPieceDataSet.h>
#include <vtkSmartPointer.h>
#include <vtkXMLMultiBlockDataWriter.h>
#include <vtkXMLPMultiBlockDataWriter.h>

using namespace std;
VtkWriter::VtkWriter(DomainCollection &dc) { this->dc = dc; }
void VtkWriter::write(string file_name, vector_type &u, vector_type &error, vector_type &resid)
{
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	vtkMultiProcessController *Controller;
	Controller = vtkMPIController::New();
	Controller->Initialize(0, 0, 1);
	set<vtkSmartPointer<vtkImageData>>   images;
	set<vtkSmartPointer<vtkDoubleArray>> arrays;
	// create MultiPieceDataSet and fill with patch information
	vtkSmartPointer<vtkXMLPMultiBlockDataWriter> writer
	= vtkSmartPointer<vtkXMLPMultiBlockDataWriter>::New();
	writer->SetController(Controller);
	vtkSmartPointer<vtkMultiBlockDataSet> block = vtkSmartPointer<vtkMultiBlockDataSet>::New();
	vtkSmartPointer<vtkMultiPieceDataSet> data  = vtkSmartPointer<vtkMultiPieceDataSet>::New();

	data->SetNumberOfPieces(dc.num_global_domains);

	auto u_view     = u.get1dView();
	auto error_view = error.get1dView();
	auto resid_view = resid.get1dView();

	for (auto &p : dc.domains) {
		Domain &d     = p.second;
		int     n     = d.n;
		double  h_x   = d.x_length / n;
		double  h_y   = d.y_length / n;
		int     start = d.id_local * n * n;

		// create image object
		vtkSmartPointer<vtkImageData> image = vtkSmartPointer<vtkImageData>::New();
		images.insert(image);
		image->SetOrigin(d.x_start, d.y_start, 0);
		image->SetSpacing(h_x, h_y, 0);
		image->SetExtent(d.x_start, d.x_start + d.x_length, d.y_start, d.y_start + d.y_length, 0,
		                 0);
		image->PrepareForNewData();
		image->SetDimensions(n + 1, n + 1, 1);

		// create solution vector
		vtkSmartPointer<vtkDoubleArray> solution = vtkSmartPointer<vtkDoubleArray>::New();
		arrays.insert(solution);
		solution->SetNumberOfComponents(1);
		solution->SetNumberOfValues(n * n);
		solution->SetName("Solution");
		for (int i = 0; i < n * n; i++) {
			solution->SetValue(i, u_view[start + i]);
		}
		image->GetCellData()->AddArray(solution);

		// create error vector
		vtkSmartPointer<vtkDoubleArray> error = vtkSmartPointer<vtkDoubleArray>::New();
		arrays.insert(error);
		error->SetNumberOfComponents(1);
		error->SetNumberOfValues(n * n);
		error->SetName("Error");
		for (int i = 0; i < n * n; i++) {
			error->SetValue(i, error_view[start + i]);
		}
		image->GetCellData()->AddArray(error);

		// create residual vector
		vtkSmartPointer<vtkDoubleArray> resid = vtkSmartPointer<vtkDoubleArray>::New();
		arrays.insert(resid);
		resid->SetNumberOfComponents(1);
		resid->SetNumberOfValues(n * n);
		resid->SetName("Residual");
		for (int i = 0; i < n * n; i++) {
			resid->SetValue(i, resid_view[start + i] * h_x * h_y);
		}
		image->GetCellData()->AddArray(resid);

		// add image to dataset
		data->SetPiece(d.id_global, image);
	}

	block->SetNumberOfBlocks(1);
	block->SetBlock(0, data);

	string out = file_name + ".vtmb";
	writer->SetFileName(out.c_str());
	writer->SetNumberOfPieces(1);
	writer->SetStartPiece(0);
	writer->SetInputData(block);
	if (rank == 0) {
		writer->SetWriteMetaFile(1);
	}
	writer->Update();

	// write data
	writer->Write();
	Controller->Delete();
}
