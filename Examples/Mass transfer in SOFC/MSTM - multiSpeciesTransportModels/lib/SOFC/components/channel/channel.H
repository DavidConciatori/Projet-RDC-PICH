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
    Foam::SOFC::channel

Description
    channel (anode side and cathode side) class

SourceFiles
    channel.C

\*---------------------------------------------------------------------------*/

#ifndef channel_H
#define channel_H

#include "fluidComponent.H"
#include "solidComponent.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace SOFC
{
  
/*---------------------------------------------------------------------------*\
                           Class component Declaration
\*---------------------------------------------------------------------------*/  

class channel
:
    virtual public fluidComponent
{

public:
          
    //- Constructor
    
        //- Construct from components
        channel
        (
          const word name,
          regionPropertiesSOFC& regions
        );

        
    //- Destructor
        
        ~channel() {}


    // Member Functions
    
        //- Set the boundary condition of mass transport variables
        //- from neighbours components
        void setMassTransportBC(fluidComponent& electrode);

        //- Solve the mass transport problem
        void solveMassTransport(fluidComponent& electrode);
        
        
    // I-O        
        
        void write()
        {
            if (activated_)
            {
              fluidComponent::write();
            }
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace SOFC

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
