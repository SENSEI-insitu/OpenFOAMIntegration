/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

Class
    Foam::functionObjects::senseiFunctionObject

Description
    TODO - update

Usage
    TODO - update
    \table
        Property          | Description               | Required | Default
        type              | type name: senseiFunctionObject      | yes |
        restartOnRestart  | Restart the averaging on restart     | no  | no
        restartOnOutput   | Restart the averaging on output      | no  | no
        periodicRestart   | Periodically restart the averaging   | no  | no
        restartPeriod     | Periodic restart period              | conditional |
        fields            | list of fields and averaging options | yes |
    \endtable


Note
    TODO - update

See also
    TODO - update
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObject

SourceFiles
    senseiFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef SENSEI_FUNCTION_OBJECT_H
#define SENSEI_FUNCTION_OBJECT_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "wallPolyPatch.H"
#include "volFieldsFwd.H"
#include "HashSet.H"

#include <ConfigurableAnalysis.h>
#include <SVTKDataAdaptor.h>

#include <svtkSmartPointer.h>
#include <svtkUnstructuredGrid.h>
#include <svtkMultiBlockDataSet.h> //gaatmp

#include <fenv.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

// Forward declaration of classes
//class fieldAverageItem;

/*---------------------------------------------------------------------------*\
                        Class senseiFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class senseiFunctionObject
:
    public fvMeshFunctionObject
{
protected:

    // Protected data

        //- Initialised flag
        bool initialised_;

	word config_file_;

	svtkSmartPointer<sensei::SVTKDataAdaptor> data_adaptor_;
	svtkSmartPointer<sensei::ConfigurableAnalysis> analysis_adaptor_;

    // OpenFOAM specifics

        //- names the volume fields that should be processed
        wordList fields_;

        //- connectivity for volume
        labelList volConn_;

        //- neighbours
        labelList neigh_;

        //- boundaries enabled
        bool boundaries_;

        //- face type list for boundaries
        Map< labelList > bndTypes_;
        //- face size (length) for boundaries
        Map< labelList > bndSizes_;
        //- connectivity for boundaries
        Map< labelList > bndConns_;

	//- floating point environment
	fenv_t fp_env;

    // Protected Member Functions

	//Foam::label calcVolConnectivitySize();
	//Foam::label calcBndConnectivitySize(const fvPatch& p);

	//void buildVolConnectivity();
	//void buildBndConnectivity();
	//void buildConnectivity();

	//void addCoordinates(svtkUnstructuredGrid *ug);
	//void addConnectivity(svtkUnstructuredGrid *ug);

	void addMesh(svtkUnstructuredGrid *ug);
	void addFields(svtkUnstructuredGrid *ug);

	void suspendFPE();
	void restoreFPE();

    // Initialisation routines

	//- Reset lists (clear existing values) and initialize averaging.
	//  Check requested field averages are valid, populate field lists
	void initialize();

public:

    //- Runtime type information
    TypeName("senseiFunctionObject");


    // Constructors

        //- Construct from Time and dictionary
        senseiFunctionObject
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        senseiFunctionObject(const senseiFunctionObject&) = delete;


    //- Destructor
    virtual ~senseiFunctionObject();


    // Member Functions

        //- Read the field average data
        virtual bool read(const dictionary&);

        //- Calculate the field averages
        virtual bool execute();

        //- Write the field averages
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const senseiFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

#endif
