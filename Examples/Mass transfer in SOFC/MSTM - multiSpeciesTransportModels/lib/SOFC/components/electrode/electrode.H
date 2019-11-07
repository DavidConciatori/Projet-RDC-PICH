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
    Foam::SOFC::electrode

Description
    Electrode (anodie side and cathode side) class

SourceFiles
    electrode.C

\*---------------------------------------------------------------------------*/

#ifndef electrode_H
#define electrode_H

#include "fluidComponent.H"
#include "solidComponent.H"
#include "electrolyte.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace SOFC
{
  
/*---------------------------------------------------------------------------*\
                           Class component Declaration
\*---------------------------------------------------------------------------*/  

class electrode
:
    public solidComponent,
    
    public fluidComponent
{

public:
              
    //- Constructor
    
        //- Construct from components
        electrode
        (
            const word name,
            regionPropertiesSOFC& regions,
            bool  reactive=0
        );

        
    //- Destructor
        
        ~electrode() {}


    // Member Functions

        //- Set the molar flux boundary values reading the current density field
        //- from the electrolyte
        void setMolarFluxBC
        (
            electrolyte& electrolyteComponent,
            word fromPatchName,
            word toPatchName
        );

        //- Solve the mass transport problem
        void solveMassTransport
        (
            fluidComponent& channel,
            electrolyte& electrolyte,
            word fromPatchName
        );

        //- Set the boundary condition of mass transport variables
        //- from neighbours components
        void setMassTransportBC
        (
            fluidComponent& channel,
            electrolyte& electrolyte,
            word fromPatchName
        );

        
    // I-O
        
        void write()
        {
            if (activated_)
            {
              fluidComponent::write();
              solidComponent::write();
            }
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace SOFC

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
