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
    Foam::SOFC::activationOverpotentials

Description
    Activation overpotential for a set of electrochemical reactions
    
SourceFiles
    activationOverpotentials.C

\*---------------------------------------------------------------------------*/

#ifndef activationOverpotentials_H
#define activationOverpotentials_H

#include "activationOverpotential.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  
namespace SOFC
{ 

/*---------------------------------------------------------------------------*\
                           Class activationOverpotentials Declaration
\*---------------------------------------------------------------------------*/

class activationOverpotentials

{
    // Private data

        //- Numerical tolerance
        static const scalar tol_;

        //- Maximum number of iterations
        static const scalar maxIter_;        


public:

    // Public data
    
        //- List of reactions
        wordList reactions_;
  
        //- List of ButtlerVolmer functions
        PtrList<activationOverpotential> etaAct_;
  

      //- Runtime type information
      TypeName("activationOverpotentials");


    // Constructors

        //- Construct from components
        activationOverpotentials(wordList reactions, const dictionary& dic);
        

        // Destructor
        virtual ~activationOverpotentials()
        {}


    // Member Functions
    
        //- Return scalarField of moral fraction by name
        scalarField& findSpecies(word name, PtrList<volScalarField>& x);
        
        //- Update the overpotential
        void update(scalarField& T);
    
        //- Initializate the molar fractions
        void initX(PtrList<volScalarField>& x);
          
        //- Evaluate the inverse of the function and return the result (return eta)
        inline scalar eta(scalar T, scalar I, int index, scalar etaStart = 0) const
        {
            scalar etaOld = etaStart;
            scalar etaNew = etaStart;
            scalar etaTol = tol_;
            int    iter = 0;          
          
            do
            {
                etaOld = etaNew;
                scalar alpha = 0;
                scalar beta = 0;
                
                forAll(etaAct_, j)
                {
                    alpha += etaAct_[j].BV_.i(etaAct_[j].i0_[index], T, etaOld);
                    beta += etaAct_[j].BV_.derivate(etaAct_[j].i0_[index], T, etaOld);
                }
                                
                etaNew = etaOld - (alpha - I)/(beta+Foam::SMALL);

                if (iter++ > maxIter_)
                {
                    FatalErrorIn
                    (
                        "activationOverpotentials::eta(scalar i0, scalar T, "
                        " scalar I, scalar etaStart = 0) const"
                    )   << "Maximum number of iterations exceeded"
                        << abort(FatalError);
                }
            } while (mag(etaNew - etaOld) > etaTol);  
          
            return etaNew;
        }   
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace SOFC

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
