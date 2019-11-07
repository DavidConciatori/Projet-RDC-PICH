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
    Foam::SOFC::concentrationOverpotential

Description
    Concentration overpotential for a signle species

SourceFiles
    concentrationOverpotential.C

\*---------------------------------------------------------------------------*/

#ifndef concentrationOverpotential_H
#define concentrationOverpotential_H

#include "dictionary.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  
namespace SOFC
{ 

/*---------------------------------------------------------------------------*\
                           Class concentrationOverpotential Declaration
\*---------------------------------------------------------------------------*/

class concentrationOverpotential

{
    // Private data
    
        //- Perfect gas constant
        static const scalar R;  

        //- Faraday constant
        static const scalar F;    
    
        //- Species names
        wordList species_;
        
        //- Reference molar fractions
        scalarList xRef_;
        
        //- Concentration overpotential coefficients
        scalar n_;
        scalarList stechCoeff_;
        
        //- Molar fractions
        PtrList<scalarField> x_;

        
public:
    
      //- Runtime type information
      TypeName("concentrationOverpotential");


    // Constructors

        //- Construct from components
        concentrationOverpotential(const dictionary& dic);
        

        // Destructor
        virtual ~concentrationOverpotential()
        {}


    // Member Functions
    
        //- Return the species of the reaction
        const wordList& species() const
        { 
            return species_;  
        }
        
        //- Return the number of electrons
        const scalar& n() const
        { 
            return n_;  
        }
        
        //- Return the stechiometric coefficients
        const scalarList& stechCoeff() const
        { 
            return stechCoeff_;  
        }
        
        //- Return the Nernst voltage
        scalar NernstVoltage(scalar T);
        
        //- Return the overpotential
        scalarField eta(const scalarField& T) const;

        //- Read the molar fractions
        void readX(int index, scalarField& x, scalar xRef);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace SOFC

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
