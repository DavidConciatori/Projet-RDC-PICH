/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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
    Foam::Crofer22APU

Description
    Crofer22APU class

SourceFiles
    Crofer22APU.C

\*---------------------------------------------------------------------------*/

#ifndef Crofer22APU_H
#define Crofer22APU_H

#include "material.H"
#include "SOFCteamFunc4.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Crofer22APU Declaration
\*---------------------------------------------------------------------------*/

class Crofer22APU
:
    public material
{
    // Private data

	//- Electrical conductivity [S/m]
	SOFCteamFunc4 sigma_;

	
public:

    // Runtime type information
    TypeName("Crofer22APU");
   
    
    // Constructors
    
	//- Construct null
	Crofer22APU();
	
	//- Construct from components
	Crofer22APU
	(
	    const material& m,
	    const SOFCteamFunc4& sigma
	);
	
	//- Construct from Istream
	Crofer22APU(Istream& is);


    // Member Functions


        //- Electrical conductivity [S/m]
        inline scalar sigma(scalar T) const
	{
	  return sigma_.f(T);
	}
	
	
    // I-O

        //- Write the function coefficients
        void writeData(Ostream& os) const
        {
            material::writeData(os); os << nl;
            sigma_.writeData(os); os << nl;
        }	
	
        //- Ostream Operator
        friend Ostream& operator<<(Ostream& os, const Crofer22APU& m)
        {
            m.writeData(os);
            return os;
        }	
};
    
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}	// End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
