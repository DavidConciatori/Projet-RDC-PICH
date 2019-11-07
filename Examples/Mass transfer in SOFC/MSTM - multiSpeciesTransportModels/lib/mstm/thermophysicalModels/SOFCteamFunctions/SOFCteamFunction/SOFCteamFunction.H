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
    Foam::SOFCteamFunction

Description
    Abstract base class for thermo-physical functions
    
SourceFiles
    SOFCteamFunction.C

\*---------------------------------------------------------------------------*/

#ifndef SOFCteamFunction_H
#define SOFCteamFunction_H

#include "scalar.H"
#include "IOstreams.H"
#include "typeInfo.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class SOFCteamFunction Declaration
\*---------------------------------------------------------------------------*/

class SOFCteamFunction
{

public:

    //- Runtime type information
    TypeName("SOFCteamFunction");
    
    bool assigned_;


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            SOFCteamFunction,
            Istream,
            (Istream& is),
            (is)
        );


    // Constructors

        //- Construct null
        SOFCteamFunction()
        {
            assigned_ = 1;
        }

        //- Return a pointer to a new SOFCteamFunction created from input
        static autoPtr<SOFCteamFunction> New(Istream& is);


    // Destructor

        virtual ~SOFCteamFunction()
        {}


    // Member Functions

        //- Evaluate the function and return the result
        virtual scalar f(scalar T) const = 0;

        //- Write the function coefficients
        virtual void writeData(Ostream& os) const = 0;


    // Ostream Operator

        friend Ostream& operator<<(Ostream& os, const SOFCteamFunction& f)
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
