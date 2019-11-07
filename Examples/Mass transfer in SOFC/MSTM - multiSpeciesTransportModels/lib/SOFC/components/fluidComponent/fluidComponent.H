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
    Foam::SOFC::fluidComponent

Description
    A base class for fluid component (channel) of solid oxide fuel cell
    
SourceFiles
    fluidComponent.C

\*---------------------------------------------------------------------------*/

#ifndef fluidComponent_H
#define fluidComponent_H

#include "component.H"
#include "multiSpeciesTransportModel.H"
#include "basicMultiComponentMixtureSOFC.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace SOFC
{
  
/*---------------------------------------------------------------------------*\
                           Class component Declaration
\*---------------------------------------------------------------------------*/  

class fluidComponent
:

    virtual public component
{

protected:

    // Protected data
    
        //- Multi species transport model
        autoPtr<multiSpeciesTransportModel> mstm_;
        

public:

    //- Constructor
    
        //- Construct from components
        fluidComponent
        (
            const word elementName,
            regionPropertiesSOFC& regions
        );
        
        
    //- Destructor
        
        virtual ~fluidComponent() {}


    // Member Functions

        //- Return the species names
        const wordList speciesName()
        {
            return (mix().species());
        }

        //- Return a pointer to the mixture
        basicMultiComponentMixtureSOFC& mix()
        {
            return (mstm().mix());
        }

        //- Return a pointer to the pressure field
        volScalarField& p()
        {
            return (mstm().thermo().p());
        }

        //- Return a pointer to the multi species transport model
        multiSpeciesTransportModel& mstm()
        {
            return mstm_();
        }
        
        //- Update the boundary condition for mixture mass flux
        void updadeBoundaryFlux(label speciesI, word patchName);

        //- Set the boundary conditions for the global flow
        //- depending on the multispecies transport
        virtual void updateBoundaryConditions(word diffusiveBoundary);
        

    // I-O
        
        void write();
        
        //- Print on screen the mass balance of the mixture
        //  through a given boundaries list
        void printMassFlux(wordList boundaryName = wordList());
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace SOFC

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
