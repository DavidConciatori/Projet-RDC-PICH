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
    Foam::material

Description
    The thermophysical properties of a material

SourceFiles
    material.C

\*---------------------------------------------------------------------------*/

#ifndef material_H
#define material_H

#include "scalar.H"
#include "IOstreams.H"
#include "typeInfo.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class material Declaration
\*---------------------------------------------------------------------------*/

class material

{
    // Private data
	
        // Density [kg/m^3]
        scalar rho_;

        //- Specific heat capacity [J/(kg.K)]
        scalar cp_;

        //- Thermal conductivity [W/(m.K)]
        scalar K_;

	
public:

    //- Runtime type information
    TypeName("material");

    // Declare run-time constructor selection tables
    
	declareRunTimeSelectionTable
	(
	    autoPtr,
	    material,
	    ,
	    (),
	    ()
	);

	declareRunTimeSelectionTable
	(
	    autoPtr,
	    material,
	    Istream,
	    (Istream& is),
	    (is)
	);


    // Constructors

        //- Construct from components 
        material
        (
            scalar rho,
	    scalar cp,
            scalar K
        )
        :
            rho_(rho),
	    cp_(cp),
	    K_(K)
        {}
	
        //- Construct from Istream
        material(Istream& is)
        :
            rho_(readScalar(is)),
            cp_(readScalar(is)),
            K_(readScalar(is))
        {}

        //- Return a pointer to a new solid created from input
        static autoPtr<material> New(word name);

    // Destructor

        virtual ~material()
        {}


    // Member Functions

        // Physical constants which define the material

            //- Density [kg/m^3]
            inline scalar rho() const
	    {
		return rho_;
	    }

	    //- Specific heat capacity [J/(kg.K)]
	    inline scalar cp() const
	    {
		return cp_;
	    }
	    
	    //- Thermal conductivity [W/(m.K)]
	    inline scalar K() const
	    {
		return K_;
	    }
	       
	    
	// Physical constants function of temperature

	    //- Electrical conductivity [S/m]
	    virtual scalar sigma(scalar T) const = 0;  
	    
	    tmp<volScalarField> sigma(volScalarField& T);
	
	    
    // I-O

        //- Write the material data
        virtual void writeData(Ostream& os) const
        {
            os  << rho_ << token::SPACE
                << cp_ << token::SPACE
		<< K_ << token::SPACE;
        }

        //- Ostream Operator
        friend Ostream& operator<<(Ostream& os, const material& s)
        {
            s.writeData(os);
            return os;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
