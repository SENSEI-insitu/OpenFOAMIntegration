/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "senseiFunctionObject.H"
#include "addToRunTimeSelectionTable.H"

#include <svtkDoubleArray.h>
//#include <svtkUnstructuredGridWriter.h>
#include <svtkMultiBlockDataSet.h>
//#include <svtkXMLMultiBlockDataWriter.h>
#include <svtkCellType.h>
#include <svtkCellData.h>
#include <svtkPointData.h>


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(senseiFunctionObject, 0);
    addToRunTimeSelectionTable(functionObject, senseiFunctionObject, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::senseiFunctionObject::initialize()
{
    initialised_ = true;
}

void Foam::functionObjects::senseiFunctionObject::suspendFPE()
{
    // Get the current FPE environment and hold FP exceptions.
    feholdexcept(&fp_env);
}
 
void Foam::functionObjects::senseiFunctionObject::restoreFPE()
{
    // Restore the last FPE environment.
    fesetenv(&fp_env);
}

void Foam::functionObjects::senseiFunctionObject::addFields(
	svtkUnstructuredGrid *ug)
{
    Log << "addFields(): begin" << endl;

    for ( const word& fieldName : fields_ )
    {
        if ( mesh_.foundObject< volVectorField >( fieldName ) )
        {
            volVectorField& vvecField =
                mesh_.lookupObjectRef< volVectorField >( fieldName );

	    Log << "    adding vector field: " << fieldName
		<< ", " << vvecField.size() << endl;

	    double *vvdata = const_cast<double*>(&vvecField[0][0]);

	    svtkSmartPointer<svtkDoubleArray> vv =
		svtkSmartPointer<svtkDoubleArray>::New();
	    vv->SetName(fieldName.c_str());
	    vv->SetNumberOfComponents(3);
	    vv->SetNumberOfTuples(vvecField.size());
	    vv->SetArray(vvdata, 3 * vvecField.size(), 1 /*save*/);
	    ug->GetCellData()->AddArray(vv);
        }

        if ( mesh_.foundObject< volScalarField >( fieldName ) )
        {
            volScalarField& vscaField =
                mesh_.lookupObjectRef< volScalarField >( fieldName );

	    Log << "    adding scalar field: " << fieldName
		<< ", " << vscaField.size() << endl;

	    double *vsdata = const_cast<double*>(&vscaField[0]);

	    svtkSmartPointer<svtkDoubleArray> vs =
		svtkSmartPointer<svtkDoubleArray>::New();
	    vs->SetName(fieldName.c_str());
	    vs->SetNumberOfComponents(1);
	    vs->SetNumberOfTuples(vscaField.size());
	    vs->SetArray(vsdata, vscaField.size(), 1 /*save*/);
	    ug->GetCellData()->AddArray(vs);
        }
    }

    Log << "addFields(): end" << endl;
}

void Foam::functionObjects::senseiFunctionObject::addMesh(
	svtkUnstructuredGrid *ug)
{
    // From Foam::vtkPVFoam::volumeVTKMesh
    Log << "addMesh(): begin" << endl;

    const cellModel& tet = *(cellModeller::lookup("tet"));
    const cellModel& pyr = *(cellModeller::lookup("pyr"));
    const cellModel& prism = *(cellModeller::lookup("prism"));
    const cellModel& wedge = *(cellModeller::lookup("wedge"));
    const cellModel& hex = *(cellModeller::lookup("hex"));

    const cellShapeList& cellShapes = mesh_.cellShapes();

    // face owner is needed to determine the face orientation
    const labelList& owner = mesh_.faceOwner();

    Log << "    mesh nCells  = " << mesh_.nCells() << endl
	<< "         nPoints = " << mesh_.nPoints() << endl;

    // Add OpenFOAM mesh vertices to VTK
    double *ptdata = const_cast<double*>(&mesh_.points()[0][0]);

    svtkSmartPointer<svtkDoubleArray> xyz =
	svtkSmartPointer<svtkDoubleArray>::New();
    xyz->SetNumberOfComponents(3);
    xyz->SetArray(ptdata, 3 * mesh_.nPoints(), 1 /* save */);

    svtkSmartPointer<svtkPoints> vtkpoints = svtkSmartPointer<svtkPoints>::New();
    vtkpoints->SetData(xyz);

    Log << "    converting cells" << endl;

    ug->Allocate(mesh_.nCells());

    // Set counters for additional points and additional cells
    label addPointi = 0, addCelli = 0;

    // Create storage for points - needed for mapping from OpenFOAM to VTK
    // data types - max 'order' = hex = 8 points
    svtkIdType nodeIds[8];

    // face-stream for a polyhedral cell
    // [numFace0Pts, id1, id2, id3, numFace1Pts, id1, id2, id3, ...]
    DynamicList<svtkIdType> faceStream(256);

#undef GAA_DEBUG
#if GAA_DEBUG
    svtkIdType np = mesh_.nPoints(); //gaatmp
#endif
    forAll(cellShapes, celli)
    {
        const cellShape& cellShape = cellShapes[celli];
        const cellModel& cellModel = cellShape.model();

        if (cellModel == tet)
        {
            for (int j = 0; j < 4; j++)
                nodeIds[j] = cellShape[j];
#if GAA_DEBUG
	    cout << "tet: " <<
		nodeIds[0] << ", " <<
		nodeIds[1] << ", " <<
		nodeIds[2] << ", " <<
		nodeIds[3] << "\n";
	    if (nodeIds[0] >= np || nodeIds[1] >= np || nodeIds[2] >= np
		    || nodeIds[3] >= np)
		cout << "ERROR: point out of range!\n";
#endif
            ug->InsertNextCell(SVTK_TETRA, 4, nodeIds);
        }
        else if (cellModel == pyr)
        {
            for (int j = 0; j < 5; j++)
                nodeIds[j] = cellShape[j];
#if GAA_DEBUG
	    cout << "pyr: " <<
		nodeIds[0] << ", " <<
		nodeIds[1] << ", " <<
		nodeIds[2] << ", " <<
		nodeIds[3] << ", " <<
		nodeIds[4] << "\n";
	    if (nodeIds[0] >= np || nodeIds[1] >= np || nodeIds[2] >= np
		    || nodeIds[3] >= np || nodeIds[4] >= np)
		cout << "ERROR: point out of range!\n";
#endif
            ug->InsertNextCell(SVTK_PYRAMID, 5, nodeIds);
        }
        else if (cellModel == prism)
        {
            // VTK has a different node order for SVTK_WEDGE
            // their triangles point outwards!
            nodeIds[0] = cellShape[0];
            nodeIds[1] = cellShape[2];
            nodeIds[2] = cellShape[1];
            nodeIds[3] = cellShape[3];
            nodeIds[4] = cellShape[5];
            nodeIds[5] = cellShape[4];

#if GAA_DEBUG
	    cout << "prism: " <<
		nodeIds[0] << ", " <<
		nodeIds[1] << ", " <<
		nodeIds[2] << ", " <<
		nodeIds[3] << ", " <<
		nodeIds[4] << ", " <<
		nodeIds[5] << "\n";
	    if (nodeIds[0] >= np || nodeIds[1] >= np || nodeIds[2] >= np
		    || nodeIds[3] >= np || nodeIds[4] >= np || nodeIds[5] >= np)
		cout << "ERROR: point out of range!\n";
#endif
            ug->InsertNextCell(SVTK_WEDGE, 6, nodeIds);
        }
        else if (cellModel == wedge)
        {
            // Treat as squeezed hex
            nodeIds[0] = cellShape[0];
            nodeIds[1] = cellShape[1];
            nodeIds[2] = cellShape[2];
            nodeIds[3] = cellShape[2];
            nodeIds[4] = cellShape[3];
            nodeIds[5] = cellShape[4];
            nodeIds[6] = cellShape[5];
            nodeIds[7] = cellShape[6];

#if GAA_DEBUG
	    cout << "wedge (degen hex): " <<
		nodeIds[0] << ", " <<
		nodeIds[1] << ", " <<
		nodeIds[2] << ", " <<
		nodeIds[3] << ", " <<
		nodeIds[4] << ", " <<
		nodeIds[5] << ", " <<
		nodeIds[6] << ", " <<
		nodeIds[7] << "\n";
	    if (nodeIds[0] >= np || nodeIds[1] >= np || nodeIds[2] >= np
		    || nodeIds[3] >= np || nodeIds[4] >= np || nodeIds[5] >= np
		    || nodeIds[6] >= np || nodeIds[7] >= np)
		cout << "ERROR: point out of range!\n";
#endif
            ug->InsertNextCell(SVTK_HEXAHEDRON, 8, nodeIds);
        }
        else if (cellModel == hex)
        {
            for (int j = 0; j < 8; j++)
                nodeIds[j] = cellShape[j];
#if GAA_DEBUG
	    cout << "hex: " <<
		nodeIds[0] << ", " <<
		nodeIds[1] << ", " <<
		nodeIds[2] << ", " <<
		nodeIds[3] << ", " <<
		nodeIds[4] << ", " <<
		nodeIds[5] << ", " <<
		nodeIds[6] << ", " <<
		nodeIds[7] << "\n";
	    if (nodeIds[0] >= np || nodeIds[1] >= np || nodeIds[2] >= np
		    || nodeIds[3] >= np || nodeIds[4] >= np || nodeIds[5] >= np
		    || nodeIds[6] >= np || nodeIds[7] >= np)
		cout << "ERROR: point out of range!\n";
#endif
            ug->InsertNextCell(SVTK_HEXAHEDRON, 8, nodeIds);
        }
        else
        {
            // Polyhedral cell - use SVTK_POLYHEDRON
            const labelList& cFaces = mesh_.cells()[celli];

            svtkIdType nFaces = cFaces.size();
            svtkIdType nLabels = nFaces;

            // count size for face stream
            forAll(cFaces, cFacei)
            {
                const face& f = mesh_.faces()[cFaces[cFacei]];
                nLabels += f.size();
            }

            // build face-stream
            // [numFace0Pts, id1, id2, id3, numFace1Pts, id1, id2, id3, ...]
            // point Ids are global
            faceStream.clear();
            faceStream.reserve(nLabels + nFaces);

            forAll(cFaces, cFacei)
            {
                const face& f = mesh_.faces()[cFaces[cFacei]];
                const bool isOwner = (owner[cFaces[cFacei]] == celli);
                const label nFacePoints = f.size();

                // number of labels for this face
                faceStream.append(nFacePoints);

                if (isOwner)
                {
                    forAll(f, fp)
                    {
                        faceStream.append(f[fp]);
                    }
                }
                else
                {
                    // fairly immaterial if we reverse the list
                    // or use face::reverseFace()
                    forAllReverse(f, fp)
                    {
                        faceStream.append(f[fp]);
                    }
                }
            }

#if GAA_DEBUG
	    {
		cout << "polyhedron: (" << nFaces << " faces) ";
		for (int i = 0; i < faceStream.size() - 1; i++)
		    cout << faceStream[i] << ", ";
		cout << faceStream[faceStream.size() - 1] << "\n";
		for (int i = 0; i < faceStream.size() - 1; i++)
		{
		    if (faceStream[i] >= np)
		    {
			cout << "ERROR: point out of range!\n";
			break;
		    }
		}
	    }
#endif
            ug->InsertNextCell(SVTK_POLYHEDRON, nFaces, faceStream.data());
        }
    }

    ug->SetPoints(vtkpoints);

    Log << "addMesh(): end" << endl;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::senseiFunctionObject::senseiFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    initialised_(false),
    fields_(),
    config_file_(dict.lookupOrDefault<string>("config_file",
		"constant/sensei_config.xml"))
{
    Log << "senseiFunctionObject(): begin" << endl;

    if (!Pstream::parRun())
    {
        Log << "    this is not a parallel run" << endl;

        int argc = 0; char ** argv= new char*[0];

        if (MPI_Init(&argc, &argv) != MPI_SUCCESS)
        {
            FatalError
                << "could not initialize mpi\n"
                << exit(FatalError);
        }

        delete[] argv;
    }

    if (!read(dict))
    {
        FatalError
            << type() << " with errors"
            << exit(FatalError);
    }

    // Create the data adaptor, send it the communicator, initialize with
    // any starting data (grid, etc.)
    data_adaptor_ = svtkSmartPointer<sensei::SVTKDataAdaptor>::New();
    //data_adaptor_->SetCommunicator(PstreamGlobals::MPI_COMM_FOAM);

    // Create the analysis adaptor, send it the communicator, initialize
    // with the configuration file.
    analysis_adaptor_ = svtkSmartPointer<sensei::ConfigurableAnalysis>::New();
    //analysis_adaptor_->SetCommunicator(PstreamGlobals::MPI_COMM_FOAM);
    suspendFPE();
    analysis_adaptor_->Initialize(config_file_.c_str());
    restoreFPE();

    Log << "senseiFunctionObject(): end" << endl;
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::senseiFunctionObject::~senseiFunctionObject()
{
    Log << "~senseiFunctionObject(): begin" << endl;
    // Clean up.
    suspendFPE();
    analysis_adaptor_->Finalize();
    restoreFPE();
    Log << "~senseiFunctionObject(): end" << endl;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::senseiFunctionObject::read(const dictionary& dict)
{
    Log << "read(): begin" << endl;

    fvMeshFunctionObject::read(dict);

    //dict.readEntry( "fields", fields_ );
    fields_ = List<word>(dict.lookup("fields"));

    if (fields_.empty())
    {
	Log << "    no fields specified" << endl;
        return false;
    }
    else
    {
	Log << "    fields:" << endl;
        for( const word& i : fields_ )
	    Log << "        " << i << endl;
    }

    initialised_ = false;

    Log << "read(): end" << endl;
    return true;
}

bool Foam::functionObjects::senseiFunctionObject::execute()
{
    Log << "execute(): begin" << endl;

    if (! initialised_)
	initialize();

    // create the cells from the face list
    svtkSmartPointer<svtkUnstructuredGrid> ug =
	svtkSmartPointer<svtkUnstructuredGrid>::New();
    addMesh(ug);
    addFields(ug);

#if 0
    {
	int par_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &par_rank);

	std::string srank = std::to_string(par_rank);
	std::string sstep = std::to_string(time_.timeIndex());

	std::string fname = "gaa_" + sstep + "_" + srank + ".vtk";

	// write file (for debugging)
	svtkSmartPointer<svtkUnstructuredGridWriter> writer =
	    svtkSmartPointer<svtkUnstructuredGridWriter>::New();
	writer->SetFileName(fname.c_str());
	writer->SetInputData(ug);
	writer->Write();
    }
#endif

    int par_rank, par_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &par_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &par_size);

    svtkSmartPointer<svtkMultiBlockDataSet> mb =
	svtkSmartPointer<svtkMultiBlockDataSet>::New();
    mb->SetNumberOfBlocks(par_size);
    mb->SetBlock(par_rank, ug);

#if 0
    {
	int par_rank, par_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &par_rank);

	std::string srank = std::to_string(par_rank);
	std::string sstep = std::to_string(time_.timeIndex());

	std::string fname = "gaa_" + sstep + "_" + srank + ".vtm";

	// write file (for debugging)
	svtkSmartPointer<svtkXMLMultiBlockDataWriter> writer =
	    svtkSmartPointer<svtkXMLMultiBlockDataWriter>::New();
	writer->SetFileName(fname.c_str());
	writer->SetInputData(mb);
	writer->Write();
    }
#endif

    int tstep = time_.timeIndex();
    double tval = time_.timeOutputValue();

    // update the data adaptor with current info
    data_adaptor_->SetDataTime(tval);
    data_adaptor_->SetDataTimeStep(tstep);
    data_adaptor_->SetDataObject(mesh_.name(), mb);
    //data_adaptor_->SetDataObject("mesh", mb);
    Log << "    adding mesh: " << mesh_.name() << endl;

    // execute the in situ analysis
    suspendFPE();
    analysis_adaptor_->Execute(data_adaptor_, nullptr);
    restoreFPE();

    // free up memory
    data_adaptor_->SetDataObject(mesh_.name(), nullptr);
    //data_adaptor_->SetDataObject("mesh", nullptr);

    Log << "execute(): end" << endl;
    return true;
}

bool Foam::functionObjects::senseiFunctionObject::write()
{
    Log << "write(): begin"<< endl;
    Log << "write(): end"<< endl;
    return true;
}


// ************************************************************************* //
