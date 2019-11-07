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

\*---------------------------------------------------------------------------*/

#include "fluidComponent.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SOFC::fluidComponent::fluidComponent
(
    const Foam::word elementName,
    Foam::regionPropertiesSOFC& regions
)
:
    component(elementName, regions)
{      
    mstm_ = multiSpeciesTransportModel::New(regions.mesh(elementName));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::SOFC::fluidComponent::updadeBoundaryFlux
(
    Foam::label speciesI,
    Foam::word patchName
)
{
    label index = mesh().boundaryMesh().findPatchID(patchName);
    scalarField& boundaryN = mix().N(speciesI).boundaryField()[index];
    scalarField& boundaryn = mix().n(speciesI).boundaryField()[index];
    scalar Wi(mstm().speciesThermo()[speciesI].W()); 
    
    forAll(boundaryn, i)
    {
      boundaryn[i] = boundaryN[i] * Wi;
    } 
}


void Foam::SOFC::fluidComponent::updateBoundaryConditions
(
    Foam::word diffusiveBoundary
)
{ 
    label index = mesh().boundaryMesh().findPatchID(diffusiveBoundary);
    if(index != -1)
    {
        mstm().mt().surfMassSource().boundaryField()[index] =
            -mstm().mt().phi().boundaryField()[index] ;
        forAll(mstm().mix().species(), speciesI)
        {
            mstm().mt().surfMassSource().boundaryField()[index] +=
                mstm().mix().n(speciesI).boundaryField()[index];
        }
    }
    else
    {
        Info << "\n\nBoundary " << diffusiveBoundary << "is not present\n\n";
        Foam::FatalError.exit();
    }
}


void Foam::SOFC::fluidComponent::printMassFlux(Foam::wordList boundaryName)
{
    if (activated_)
    {
        Info << "\n\nMass flux for the " << name() << "\n\n";
        
        if(boundaryName.size() == 0)
        {
            boundaryName.setSize(mesh().boundaryMesh().size());
            forAll(boundaryName, i)
            {
                boundaryName[i] = mstm().mt().phi().mesh().boundaryMesh()[i].name();
            }
        }
      
        forAll(mstm().mix().species(), speciesI)
        {
            scalar sumPhi = 0;
            Info << "\n\n";
            forAll(boundaryName, boundaryI)
            {
                label index = mesh().boundaryMesh().findPatchID(boundaryName[boundaryI]);
                scalar massFlow = gSum(mstm().mix().n(speciesI).boundaryField()[index]);

                Info << mstm().mix().species()[speciesI]
                    << " mass flow at " 
                    << boundaryName[boundaryI]
                    << ": " << massFlow
                    << "\n";

                sumPhi += massFlow;
            }
            Info << mix().species()[speciesI] << " mass flow balance: " << sumPhi;
        }

        mstm().mt().printMassFlux(boundaryName);
        
        Info << "\n\n";
    }
}


void Foam::SOFC::fluidComponent::write()
{
    mstm().write();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
