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
    Foam::Darcy

Description
    Simplified form of the momentum equation in porous media that only
    considers the Darcy drag term.
    
SourceFiles
    Darcy.C

\*---------------------------------------------------------------------------*/

#ifndef Darcy_H
#define Darcy_H

#include "momentumTransport.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Darcy Declaration
\*---------------------------------------------------------------------------*/

class Darcy
:
    public momentumTransport
{

protected:  
  
    // Protected data
  
        //- Permeability [m^2]
        const dimensionedScalar B0_;  
  
        //- Switch on/off the Darcy term in the species equations
        bool DarcySwitch_;  
  

public:

    //- Runtime type information
    TypeName("Darcy");


    // Constructors

        //- construct from components
        Darcy
        (
            hsCombustionThermoSOFC& thermo,
            const fvMesh& mesh
        );


    // Destructor

        ~Darcy()
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
