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
    Foam::SOFC::concentrationOverpotentials

Description
    Concentration overpotential for a group of species

SourceFiles
    concentrationOverpotentials.C

\*---------------------------------------------------------------------------*/

#ifndef concentrationOverpotentials_H
#define concentrationOverpotentials_H

#include "channel.H"
#include "concentrationOverpotential.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  
namespace SOFC
{ 

/*---------------------------------------------------------------------------*\
                           Class concentrationOverpotentials Declaration
\*---------------------------------------------------------------------------*/

class concentrationOverpotentials

{

public:

    // Public data
  
        //- List of reactions
        wordList reactions_;
  
        //- List of concentration overpotential
        PtrList<concentrationOverpotential> etaConc_;
  

      //- Runtime type information
      TypeName("concentrationOverpotentials");


    // Constructors

        //- Construct from components
        concentrationOverpotentials(wordList reactions, const dictionary& dic);
        

        // Destructor
        virtual ~concentrationOverpotentials()
        {}


    // Member Functions
    
        //- Return the Nernst voltage
        scalar NernstVoltage(scalar T);
    
        //- Return scalarField of moral fraction by name
        scalarField& findSpecies(word name, PtrList<volScalarField>& x);
    
        //- Initializate the molar fractions array
        void initX(PtrList<volScalarField>& x, channel& ch);

        //- Evaluate the inverse of the function and return the result (return eta)
        scalarField eta(scalarField& T);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace SOFC

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
