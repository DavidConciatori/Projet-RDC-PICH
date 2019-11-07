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
    Foam::SIMPLE

Description
    SIMPLE algorithm for coupling velocity and pressure to solve
    the momentum equation

SourceFiles
    SIMPLE.C

\*---------------------------------------------------------------------------*/

#ifndef SIMPLE_H
#define SIMPLE_H

#include "momentumTransport.H"
#include "RASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class SIMPLE Declaration
\*---------------------------------------------------------------------------*/

class SIMPLE
:
    public momentumTransport
{

private:
  
    // Private date
    
        //- Ras turbulence model
        autoPtr<compressible::RASModel> turbulence_;
  
public:

    //- Runtime type information
    TypeName("SIMPLE");


    // Constructors

        //- construct from components
        SIMPLE
        (
            hsCombustionThermoSOFC& thermo,
            const fvMesh& mesh
        );


    // Destructor

        ~SIMPLE()
        {}


    // Member Functions

        //- Solve a single step of the momentum equation
        scalar momentumLoop();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
