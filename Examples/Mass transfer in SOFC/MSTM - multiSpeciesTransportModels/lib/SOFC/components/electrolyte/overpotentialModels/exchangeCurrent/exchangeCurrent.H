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
    Foam::SOFC::exchangeCurrent

Description
    Exchange current density for an electrochemical reaction
    
SourceFiles
    exchangeCurrent.C

\*---------------------------------------------------------------------------*/

#ifndef exchangeCurrent_H
#define exchangeCurrent_H

#include "scalar.H"
#include "scalarList.H"
#include "scalarField.H"
#include "IOstreams.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  
namespace SOFC
{ 

/*---------------------------------------------------------------------------*\
                           Class exchangeCurrent Declaration
\*---------------------------------------------------------------------------*/

class exchangeCurrent

{
    // Private data

        //- Gas constant
        static const scalar R;

        //- Species
        wordList species_;
        
        //- Exchange current density coefficients
        scalar gamma_, Eact_;
        scalarList exponent_;

        //- Molar fractions
        PtrList<scalarField> x_;


public:

      //- Runtime type information
      TypeName("exchangeCurrent");


    // Constructors

        //- Construct from components
        exchangeCurrent(wordList species, scalar gamma, scalar Eact, scalarList exponent)
        :
            species_(species),
            gamma_(gamma),
            Eact_(Eact),
            exponent_(exponent)
        {
            if (species_.size() > 2)
            {
                FatalErrorIn
                (
                    "exchangeCurrent(scalar gamma, scalar Eact, scalarList exponent)"
                )   << "Species overflow (reaction with more than 2 species)"
                    << abort(FatalError);
            }
            if (species.size() == 1)
            {
              exponent_.setSize(2);
              exponent_[1] = 1;
            }
        }

        //- Construct from dictionary
        exchangeCurrent(const dictionary& dic)
        :
            species_(dic.lookup("species")),
            gamma_(readScalar(dic.lookup("gamma"))),
            Eact_(readScalar(dic.lookup("Eact"))),
            exponent_(dic.lookup("speciesExponent"))
        {
            if (species_.size() > 2)
            {
                FatalErrorIn
                (
                    "exchangeCurrent(scalar gamma, scalar Eact, scalarList exponent)"
                )   << "Species overflow (reaction with more than 2 species)"
                    << abort(FatalError);
            }  
        }

        
        // Destructor
        virtual ~exchangeCurrent()
        {}


    // Member Functions

        //- Initializate the molar fractions array
        void initX(PtrList<scalarField>& x)
        {   
            forAll(x, i)
            {
              x_.set(i, x[i]);
            }
        }

        //- Return the species of the reaction
        const wordList species() const
        { 
            return species_;  
        }

        //- Evaluate the function and return the result for one species (return i0)
        inline scalar calc(scalar T, scalar x0) const
        { 
            return gamma_*exp(-1000*Eact_/(R*T)) * pow( x0, exponent_[0]);  
        }   
        
        //- Evaluate the function and return the result for two species (return i0)
        inline scalar calc(scalar T, scalar x0, scalar x1) const
        { 
            return gamma_*exp(-1000*Eact_/(R*T)) * pow( x0, exponent_[0]) * pow( x1, exponent_[1]);  
        }   
                
        //- Evaluate the function and return the result for one species (return i0 like a scalarField)
        scalarField calc(scalarField T, scalarField x0) const
        {
            return gamma_*exp(-1000*Eact_/(R*T)) * pow( x0, exponent_[0]);  
        }
        
        //- Evaluate the function and return the result for two species (return i0 like a scalarField)
        scalarField calc(scalarField T, scalarField x0, scalarField x1) const
        {
            return gamma_*exp(-1000*Eact_/(R*T)) * pow( x0, exponent_[0]) * pow( x1, exponent_[1]);  
        }
        
        //- Evaluate the function and return the result for one species (return i0 like a scalarField)
        scalarField calc(scalarField T, PtrList<scalarField> x) const
        {
            scalarField i0 = gamma_*exp(-1000*Eact_/(R*T));
            forAll(x, speciesI)
            {
                i0 *= pow( x[speciesI], exponent_[speciesI]);
            }
            return i0;  
        }
        
        
    // I-O

        //- Write the function coefficients
        void writeData(Ostream& os) const
        {
            os  << gamma_ << token::SPACE
                << Eact_ << token::SPACE
                << exponent_;
        }


    // Ostream Operator

        friend Ostream& operator<<(Ostream& os, const exchangeCurrent& f)
        {
            f.writeData(os);
            return os;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace SOFC

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
