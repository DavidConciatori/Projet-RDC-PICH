/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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
    Foam::SOFCteamFunc3

Description
    Interpolation function given by    
    @f[
        a \ctod exp(-b / T)
    @f]
    
SourceFiles
    SOFCteamFunc3.C    

\*---------------------------------------------------------------------------*/

#ifndef SOFCteamFunc3_H
#define SOFCteamFunc3_H

#include "SOFCteamFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class SOFCteamFunc3 Declaration
\*---------------------------------------------------------------------------*/

class SOFCteamFunc3
:
    public SOFCteamFunction
{

public:

    //- SOFC team function 3 coefficients
    scalar a_, b_;
    
    //- Runtime type information
    TypeName("SOFCteamFunc3");


    // Constructors

        //- Construct from components
        SOFCteamFunc3(scalar a, scalar b)
        :
            a_(a),
            b_(b)
        {}

        //- Construct from Istream
        SOFCteamFunc3(Istream& is)
        :
            a_(readScalar(is)),
            b_(readScalar(is))
        {}


    // Member Functions

        //- Evaluate the function and return the result
        inline scalar f(scalar T) const
        {
            return a_ * exp(-b_/T);
        }

        //- Return the function coefficients
        inline List<scalar> coeff()
        {
          List<scalar> cf(2);
          cf[0] = a_;	 cf[1] = b_;
          return cf;
        }

        //- Write the function coefficients
        void writeData(Ostream& os) const
        {
            os  << a_ << token::SPACE
                << b_ ;
        }


    // Ostream Operator

        friend Ostream& operator<<(Ostream& os, const SOFCteamFunc3& f)
        {
            f.writeData(os);
            return os;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
