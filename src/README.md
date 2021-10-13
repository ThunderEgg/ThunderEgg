A Quick Overview                         {#mainpage}
================

Documentation that will appear on the main page

Mesh Metadata
=============

The metadata for a tree based mesh is stored in a \ref ThunderEgg::Domain "Domain" object. 
The Domain object stores a collection of \ref ThunderEgg::PatchInfo "PatchInfo" objects.
The PatchInfo contains the metadata associated with a specific patch, 
such as information about neighboring patches, the number of cells in a patch, processor-local index, etc.

Vectors
=======
%ThunderEgg expects data to be laid out in a series of logically Cartesian patches (i.e. patch-level data is in a strided layout).
A \ref ThunderEgg::Vector "Vector" object represents the data for an entire Domain. 
A Vector contains data for each cell of the patch, and each cell can have multiple components. 
Users can retrieve a \ref ThunderEgg::View "PatchView" for each patch.

An example of looping over a 2D Vector's data:

	Domain<2> domain;
	Vector<2> x;
	for(int i = 0; i < domain.getNumLocalPatches(); i++){
		PatchInfo<2> patch_info = domain.getPatchInfoVector()[i];
		int nx = patch_info.ns[0];
		int ny = patch_info.ns[1];

		PatchView<double, 2> view = x.getPatchView(i);
		for(int c = 0; c < x.getNumComponents(); c++){
			for(int y = 0; c < ny; y++){
				for(int x = 0; x < nx; x++){
					view(x,y,c) = 1.23;
				}
			}
		}
	}

Patch Solver
============
patch operator something
patch solver something

Geometric Multigrid
===================
if you have a patch operator and a patch solver you can create a Geometric Multigrid preconditioner.
