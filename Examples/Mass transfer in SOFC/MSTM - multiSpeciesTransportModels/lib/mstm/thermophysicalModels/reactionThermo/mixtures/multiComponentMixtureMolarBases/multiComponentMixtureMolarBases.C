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

#include "multiComponentMixtureMolarBases.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ThermoType>
const ThermoType& Foam::multiComponentMixtureMolarBases<ThermoType>::constructSpeciesData
(
    const dictionary& thermoDict
)
{
    forAll(species_, i)
    {
        speciesData_.set
        (
            i,
            new ThermoType(thermoDict.lookup(species_[i]))
        );
    }

    return speciesData_[0];
}


template<class ThermoType>
void Foam::multiComponentMixtureMolarBases<ThermoType>::constructMixtureComposition
(
    const fvMesh& mesh
)
{
    forAll(species_, i)
    {        
        IOobject header
        (
            "x_" + species_[i],
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ
        );

        // check if field exists and can be read
        if (header.headerOk())
        {
            x_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "x_" + species_[i],
                        mesh.time().timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh
                )
            );
        }
        else
        {
            volScalarField xDefault
            (
                IOobject
                (
                    "xDefault",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );

            x_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "x_" + species_[i],
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                xDefault
                )
            );
        }

        y_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "y_" + species_[i],
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("y", dimless, 0.0)
            )
        );

        N_.set
        (
            i,
            new surfaceScalarField
            (
                IOobject
                (
                    "N_" + species_[i],
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("N", dimMoles/dimTime, 0.0)
            )
        );
        
        n_.set
        (
            i,
            new surfaceScalarField
            (
                IOobject
                (
                    "n_" + species_[i],
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("n", dimMass/dimTime, 0.0)
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::multiComponentMixtureMolarBases<ThermoType>::multiComponentMixtureMolarBases
(
    const dictionary& thermoDict,
    const wordList& specieNames,
    const HashPtrTable<ThermoType>& specieThermoData,
    const fvMesh& mesh
)
:
    basicMultiComponentMixtureSOFC(thermoDict, specieNames, mesh),
    speciesData_(species_.size()),
    mixture_("mixture", *specieThermoData[specieNames[0]])
{
    forAll(species_, i)
    {
        speciesData_.set
        (
            i,
            new ThermoType(*specieThermoData[species_[i]])
        );
    }

    constructMixtureComposition(mesh);
    correct();
}


template<class ThermoType>
Foam::multiComponentMixtureMolarBases<ThermoType>::multiComponentMixtureMolarBases
(
    const dictionary& thermoDict,
    const fvMesh& mesh
)
:
    basicMultiComponentMixtureSOFC(thermoDict, thermoDict.lookup("species"), mesh),
    speciesData_(species_.size()),
    mixture_("mixture", constructSpeciesData(thermoDict))
{
    constructMixtureComposition(mesh);
    correct();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::multiComponentMixtureMolarBases<ThermoType>::correct()
{  
    // Update the molar weight of the mixture
    scalarField& WCells = W_.internalField();
    forAll(WCells, celli)
    {
        WCells[celli] = x_[0].internalField()[celli] * speciesData_[0].W();
        for(label i = 1; i < x_.size(); i++)
        {
            WCells[celli] +=
                x_[i].internalField()[celli] * speciesData_[i].W();
        }
    }
    forAll(W_.boundaryField(), patchi)
    {
        fvPatchScalarField& pW = W_.boundaryField()[patchi];
        forAll(pW, facei)
        {
            pW[facei] =
                x_[0].boundaryField()[patchi][facei] * speciesData_[0].W();
            for(label i = 1; i < x_.size(); i++)
            {
                pW[facei] +=
                    x_[i].boundaryField()[patchi][facei] * speciesData_[i].W();
            }
        }
    }
    
    // Update the mass fractions
    forAll(y_, speciesI)
    {        
        const scalarField& xCells = x_[speciesI].internalField();
        const scalarField& WCells = W_.internalField();            
        scalarField& yCells = y_[speciesI].internalField();
        forAll(yCells, celli)
        {
            yCells[celli] =
                xCells[celli] / WCells[celli] * speciesData_[speciesI].W();
        }
        forAll(y_[speciesI].boundaryField(), patchi)
        {
            const fvPatchScalarField& px =
                x_[speciesI].boundaryField()[patchi];
            const fvPatchScalarField& pW = W_.boundaryField()[patchi];            
            fvPatchScalarField& py = y_[speciesI].boundaryField()[patchi];
            forAll(px, facei)
            {
                py[facei] =
                    px[facei] / pW[facei] * speciesData_[speciesI].W();
            }
        }
    }
}


template<class ThermoType>
void Foam::multiComponentMixtureMolarBases<ThermoType>::correctMassFluxes()
{ 
    forAll(n_, i)
    {        
        forAll(n_[i].internalField(), cellI)
        {
            n_[i].internalField()[cellI] =
            N_[i].internalField()[cellI] * speciesData_[i].W();
        }
        forAll(n_[i].boundaryField(), boundaryI)
        {
            forAll(n_[i].boundaryField()[boundaryI], faceI)
            {
                n_[i].boundaryField()[boundaryI][faceI] =
                    N_[i].boundaryField()[boundaryI][faceI]
                    * speciesData_[i].W();
            }
        }
    }
}


template<class ThermoType>
const ThermoType& Foam::multiComponentMixtureMolarBases<ThermoType>::cellMixture
(
    const label celli
) const
{
    mixture_ = y_[0][celli]/speciesData_[0].W()*speciesData_[0];

    for (label n=1; n<y_.size(); n++)
    {
        mixture_ += y_[n][celli]/speciesData_[n].W()*speciesData_[n];
    }

    return mixture_;
}


template<class ThermoType>
const ThermoType& Foam::multiComponentMixtureMolarBases<ThermoType>::patchFaceMixture
(
    const label patchi,
    const label facei
) const
{
    mixture_ =
        y_[0].boundaryField()[patchi][facei]
       /speciesData_[0].W()*speciesData_[0];

    for (label n=1; n<y_.size(); n++)
    {
        mixture_ +=
            y_[n].boundaryField()[patchi][facei]
           /speciesData_[n].W()*speciesData_[n];
    }

    return mixture_;
}


template<class ThermoType>
void Foam::multiComponentMixtureMolarBases<ThermoType>::read
(
    const dictionary& thermoDict
)
{
    forAll(species_, i)
    {
        speciesData_[i] = ThermoType(thermoDict.lookup(species_[i]));
    }
}


// ************************************************************************* //
