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
    Foam::SOFC::activationOverpotential

Description
    Activation overpotential for single electrochemical semi-reaction
    
SourceFiles
    activationOverpotential.C

\*---------------------------------------------------------------------------*/

#ifndef activationOverpotential_H
#define activationOverpotential_H

#include "ButtlerVolmer.H"
#include "exchangeCurrent.H"
#include "dictionary.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  
namespace SOFC
{ 

/*---------------------------------------------------------------------------*\
                           Class activationOverpotential Declaration
\*---------------------------------------------------------------------------*/

class activationOverpotential

{
  
    // Private data
    
        //- Species names
        wordList species_;
        
        //- Molar fractions
        PtrList<scalarField> x_;
        
        //- Exchange current density model
        exchangeCurrent i0model;

        
public:
    
    // Public data
    
        //- ButtlerVolmer functions
        ButtlerVolmer BV_;
    
        //- Exchange current density
        scalarField i0_;    
    
        //- Current density
        scalarField i_;
  

      //- Runtime type information
      TypeName("activationOverpotential");


    // Constructors

        //- Construct from components
        activationOverpotential(const dictionary& dic);
       

        // Destructor
        virtual ~activationOverpotential()
        {}


    // Member Functions
    
        //- Return the species of the reaction
        const wordList species() const
        { 
            return species_;  
        }
        
        //- Update the overpotential
        void update(scalarField& T);
        
        //- Return the current density
        scalarField i(const scalarField& T, const scalarField& eta) const;

        //- Read the molar Fractions
        void readX(int index, scalarField& x);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace SOFC

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
