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
    Foam::SOFC::ButtlerVolmer

Description
    Buttler-Volmer function 

SourceFiles
    ButtlerVolmer.C

\*---------------------------------------------------------------------------*/

#ifndef ButtlerVolmer_H
#define ButtlerVolmer_H

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
                           Class ButtlerVolmer Declaration
\*---------------------------------------------------------------------------*/

class ButtlerVolmer

{
    // Private data

        //- Faraday constant
        static const scalar F;
        
        //- Perfect gas constant
        static const scalar R;
        
        //- Numerical tolerance
        static const scalar tol_;
        
        //- Maximum number of iterations
        static const scalar maxIter_;        

        // Buttler-Volmer coefficients
        scalar theta1_, theta2_;


public:

      //- Runtime type information
      TypeName("ButtlerVolmer");


    // Constructors

        //- Construct from components
        ButtlerVolmer(scalar theta1, scalar theta2)
        :
            theta1_(theta1),
            theta2_(theta2)
        {}
        
        //- Construct from scalarList
        ButtlerVolmer(scalarList sl)
        :
            theta1_(sl[0]),
            theta2_(sl[1])
        {}

        //- Construct from Istream
        ButtlerVolmer(Istream& is)
        :
            theta1_(readScalar(is)),
            theta2_(readScalar(is))
        {}
   
   
        // Destructor
        virtual ~ButtlerVolmer()
        {}


    // Member Functions

        //- Evaluate the function and return the result (return i)
        inline scalar i(scalar i0, scalar T, scalar eta) const
        {
            return i0*(Foam::exp((theta1_*F)/(R*T)*eta)-Foam::exp((-theta2_*F)/(R*T)*eta));  
        }
        
        //- Evaluate the function and return the result (return i like a scalarField)
        scalarField i(const scalarField& i0, const scalarField& T, const scalarField& eta) const
        {
            return i0*(Foam::exp((theta1_*F)/(R*T)*eta)-Foam::exp((-theta2_*F)/(R*T)*eta));  
        }  
     
        //- Evaluate the derivate of the function and return the result 
        inline scalar derivate(scalar i0, scalar T, scalar eta) const
        {
            return i0*((theta1_*F)/(R*T)*Foam::exp((theta1_*F)/(R*T)*eta)+(theta2_*F)/(R*T)*Foam::exp((-theta2_*F)/(R*T)*eta));  
        }

        //- Evaluate the inverse of the function and return the result (return eta)
        inline scalar eta(scalar i0, scalar T, scalar I, scalar etaStart = 0) const
        {
            scalar etaOld = etaStart;
            scalar etaNew = etaStart;
            scalar etaTol = tol_;
            int    iter = 0;          
          
            do
            {
                etaOld = etaNew;                
                etaNew = etaOld - (i(i0, T, etaOld) - I)/(derivate(i0, T, etaOld)+Foam::SMALL);

                if (iter++ > maxIter_)
                {
                    FatalErrorIn
                    (
                        "ButtlerVolmer::eta(scalar i0, scalar T, "
                        " scalar I, scalar etaStart = 0) const"
                    )   << "Maximum number of iterations exceeded"
                        << abort(FatalError);
                }

            } while (mag(etaNew - etaOld) > etaTol);  
          
            return etaNew;
        }


    // I-O

        //- Write the function coefficients
        void writeData(Ostream& os) const
        {
            os  << theta1_ << token::SPACE
                << theta2_;
        }


    // Ostream Operator

        friend Ostream& operator<<(Ostream& os, const ButtlerVolmer& f)
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
