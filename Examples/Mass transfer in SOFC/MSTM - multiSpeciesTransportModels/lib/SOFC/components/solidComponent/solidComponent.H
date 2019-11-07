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
    Foam::SOFC::CFDsolidComponent

Description
    A base class for solid component (electrode, electrolite) of
    solid oxide fuel cell
    
SourceFiles
    solidComponent.C

\*---------------------------------------------------------------------------*/

#ifndef solidComponent_H
#define solidComponent_H

#include "component.H"
#include "material.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace SOFC
{
  
/*---------------------------------------------------------------------------*\
                           Class component Declaration
\*---------------------------------------------------------------------------*/  

class solidComponent
:

    virtual public component
{

public:

    // Public data  

        //- Solid
        autoPtr<material> solid_;
        
        //- Current density
        volScalarField i_;

        
    //- Constructor
    
        //- Construct from components
        solidComponent
        (
            const word elementName,
           regionPropertiesSOFC& regions
        );

        
    //- Destructor
        
        virtual ~solidComponent() {}


    // Member Functions

        //- Write on file
        void write();

        //- Return the solid
        autoPtr<material>& solid()
        {
            return(solid_);
        }

        //- Return the current density
        volScalarField& i()
        {
            return(i_);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace SOFC

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
