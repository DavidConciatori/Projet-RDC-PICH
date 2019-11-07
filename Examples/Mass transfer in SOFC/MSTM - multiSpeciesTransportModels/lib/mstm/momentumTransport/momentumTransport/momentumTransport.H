/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::momentumTransport

Description
    An abstract base class for momentum (and mass) transports equations
    based on a thermo package.

SourceFiles
    momentumTransportI.H
    momentumTransport.C
    newDiffusivityModel.C

\*---------------------------------------------------------------------------*/

#ifndef momentumTransport_H
#define momentumTransport_H

#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "fvMatricesFwd.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

#include "volFields.H"
#include "surfaceFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "IOdictionary.H"

#include "hsCombustionThermoSOFC.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class momentumTransport Declaration
\*---------------------------------------------------------------------------*/

class momentumTransport
{

protected:

    // Protected data

        //- Thermo package
        hsCombustionThermoSOFC& thermo_;
          
        //- Velocity field [m/s]
        volVectorField U_;

        //- Flux based on mass average velocity [kg/s]
        surfaceScalarField phi_;

        //- Density of the fluid [kg/m^3]
        volScalarField rho_;
        
        //- Surface mass source (i.e. mass source through a porous boundary
        //  where mass enters without changing the momentum) [kg/s]
        surfaceScalarField surfMassSource_;     

        //- Volumetric momentum source (i.e. momentum source
	//  due to the reactions in a multispecies mixture) [???]
//        volScalarField momentumSource_;

        
        scalar convergenceCriterion_;
        
        //- Number of non orthogonal mesh corrections
        label nNonOrthCorr_;


    // Private Member Functions

        //- Disallow copy construct
        momentumTransport(const momentumTransport&);

        //- Disallow default bitwise assignment
        void operator=(const momentumTransport&);


public:

    //- Runtime type information
    TypeName("momentumTransport");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            momentumTransport,
            dictionary,
            (
                hsCombustionThermoSOFC& thermo,
                const fvMesh& mesh
            ),
            (thermo, mesh)
        );


    // Selectors

        //- Return a reference to the selected momentum transport model
        //  selecting it in the standard transportProperties dictionary 
        static autoPtr<momentumTransport> New
        (
            hsCombustionThermoSOFC& thermo,
            const fvMesh& mesh
        );

        
        //- Return a reference to the selected momentum transport model
        //  selecting it by name 
        static autoPtr<momentumTransport> New
        (
            hsCombustionThermoSOFC& thermo,
            const fvMesh& mesh,
            const word momentumTransportName
        );


    // Constructors

        //- Construct from components
        momentumTransport
        (
            hsCombustionThermoSOFC&,
            const fvMesh& mesh
        );


    // Destructor

        virtual ~momentumTransport()
        {}


    // Member Functions
    
        //- Return const-access to the velocity
        inline const volVectorField& U() const;

        //- Return the velocity
        inline volVectorField& U();
    
        //- Return const-access to the flux based on mass average velocity
        inline const surfaceScalarField& phi() const;

        //- Return the flux based on mass average velocity
        inline surfaceScalarField& phi();

        //- Return const-access to the density
        inline const volScalarField& rho() const;

        //- Return the density
        inline volScalarField& rho();

        //- Return const-access to the surface mass source
        inline const surfaceScalarField& surfMassSource() const;
    
        //- Return the surface mass source
        inline surfaceScalarField& surfMassSource();
        
        //- Print on screen the mass balance through the desired boundaries
        void printMassFlux(wordList boundaryName);

        void momentumSolve(label nitMax);
        
        //- Solve a single step of the momentum equation
        virtual scalar momentumLoop() = 0;
        
        //- Write on file the fluid variables
        void write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "momentumTransportI.H"

#endif

// ************************************************************************* //
