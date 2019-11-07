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
    Foam::dustyGas

Description
    Dusty gas multispecies mass transport model.

SourceFiles
    dustyGas.C

\*---------------------------------------------------------------------------*/

#ifndef dustyGas_H
#define dustyGas_H

#include "multiSpeciesTransportModel.H"
#include "multiComponentMixtureMolarBases.H"
#include "Knudsen.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class MaxwellStefan Declaration
\*---------------------------------------------------------------------------*/

class dustyGas
:
    public multiSpeciesTransportModel
{
  
protected:
  
    // Protected data
  
        //- Permeability [m^2]
        const dimensionedScalar B0_;  
  
        //- Switch on/off the Darcy term in the species equations
        bool DarcySwitch_;
        
        //- Flux given by the Darcy term
        surfaceScalarField phiDarcy_;

      
    // Protected member functions	
    
        //- Update the coefficients of the species equations
        void updateCoefficients();       


public:     
    
    // Runtime type information
    TypeName("dustyGas");

    
    // Constructors

        //- Construct from components
        dustyGas
        (
            const fvMesh& mesh
        );

    //- Destructor
    ~dustyGas()
    {}


    // Member Functions
    
        //- Solve a single step of the species equations
        scalar massTransportLoop();

        //- Solve the mass transport
        void correct();
    
        //- Return the thermodynamic data of the species
        const PtrList<gasThermoPhysicsSOFC>& speciesThermo() const
        {
            return reinterpret_cast
                <const multiComponentMixtureMolarBases<gasThermoPhysicsSOFC>&>
                (mix()).speciesData();
        }        
    };

    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}	// End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
