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
    Foam::Fick

Description
    Fick multispecies mass transport model.
  
SourceFiles
    Fick.C

\*---------------------------------------------------------------------------*/

#ifndef Fick_H
#define Fick_H

#include "multiSpeciesTransportModel.H"
#include "multiComponentMixtureMassBases.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  
/*---------------------------------------------------------------------------*\
                             Class Fick Declaration
\*---------------------------------------------------------------------------*/

class Fick
:
    public multiSpeciesTransportModel
{  

protected:
  
    // Protected data
      
        //- Name of the species not directly calculated
        word inertSpecie_;

    // Protected member functions
    
        //- Update the coefficients of the species equations
        virtual void updateCoefficients();
    

public:

    //- Runtime type information
    TypeName("Fick");

    
    // Constructors

        //- Construct from components
        Fick
        (
            const fvMesh& mesh
        );


      //- Destructor
      ~Fick()
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
                <const multiComponentMixtureMassBases<gasThermoPhysicsSOFC>&>
                (mix()).speciesData();
        }
        
        //- Return const-access to the inert specie
        virtual const word inertSpecie() const
        {
            return inertSpecie_;
        }
  };

  
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}	// End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

