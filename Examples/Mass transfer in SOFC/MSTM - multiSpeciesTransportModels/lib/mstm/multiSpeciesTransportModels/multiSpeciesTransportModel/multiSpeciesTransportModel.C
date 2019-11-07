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

#include "multiSpeciesTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //  
  
    defineTypeNameAndDebug(multiSpeciesTransportModel, 0);
    defineRunTimeSelectionTable(multiSpeciesTransportModel, dictionary);
    
    const dimensionedScalar multiSpeciesTransportModel::RR
        ("R", dimensionSet(1, 2, -2, -1, -1), 8314.51); 
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::multiSpeciesTransportModel::constructP(const fvMesh& mesh)
{
    p_.setSize(ns_); 
    forAll(p_, i)
    {
        const volScalarField::GeometricBoundaryField& xbf =
            mix().x(i).boundaryField();

        p_.set
        (
            i, new volScalarField
            (
                IOobject
                (
                    "p_" + mix().species()[i],
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimPressure,
                xbf.types()
            )
        );
  
        forAll(p_[i].internalField(), cellI)
        {
            p_[i].internalField()[cellI] = 
                mix().x(i).internalField()[cellI] *
                thermo().p().internalField()[cellI];
        }
        forAll(p_[i].boundaryField(), patchI)
        {
            forAll(p_[i].boundaryField()[patchI], faceI)
            {
                p_[i].boundaryField()[patchI][faceI] = 
                    mix().x(i).boundaryField()[patchI][faceI] *
                    thermo().p().boundaryField()[patchI][faceI];
            }
        }
    } 
}  
  
  
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiSpeciesTransportModel::multiSpeciesTransportModel
(
    const fvMesh& mesh
)
:
    IOdictionary
    (
        IOobject
        (
            "transportProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
  
    convergenceCriterion_
    (
        mesh.solutionDict().subDict("SIMPLE")
            .lookupOrDefault<scalar>("convergence", 0)
    ),
    
    nitMax_
    (
        mesh.solutionDict().subDict("SIMPLE")
            .lookupOrDefault<scalar>("nitMax", 1000)
    )
{   
    // Pre-read the species
    wordList speciesName_
    (
        IOdictionary
        (
            IOobject
            (
                "thermophysicalProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        ).lookup("species")
    );

    ns_ = speciesName_.size();

    // Resize the pointer lists
    Gamma_.setSize(ns_);
    phi_N_.setSize(ns_);
    phi_P_.setSize(ns_);
    
    // Create the source fields
    Sy_.setSize(ns_);
    Sx_.setSize(ns_);
    forAll(Sy_, i)
    {
        Sy_.set
        (
            i, new volScalarField
            (
                IOobject
                (
                    "Sy_" + speciesName_[i],
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("Sy", dimMass/dimTime/dimVol, 0)
            )
        );

        Sx_.set
        (
            i, new volScalarField
            (
                IOobject
                (
                    "Sx_" + speciesName_[i],
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("Sx", dimMoles/dimTime/dimVol, 0)
            )
        );
    }
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::multiSpeciesTransportModel>
Foam::multiSpeciesTransportModel::New
(
    const fvMesh& mesh
)
{
    word modelName;

    // Enclose the creation of the dictionary to ensure it is deleted
    // before the multiSpeciesTransportModel is created otherwise the
    // dictionary is entered in the database twice
    {
        IOdictionary dict
        (
            IOobject
            (
                "transportProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        dict.lookup("multispeciesTransportModel") >> modelName;
    }
  
    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "DiffusivityModel::New(const volVectorField&, "
            "const surfaceScalarField&)"
        )   << "Unknown diffusivityModel type "
            << modelName << endl << endl
            << "Valid  diffusivityModels are : " << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
  }

  return autoPtr<multiSpeciesTransportModel>
      (cstrIter()(mesh));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::word Foam::multiSpeciesTransportModel::inertSpecie() const
{
    FatalErrorIn("Foam::multiSpeciesTransportModel::inertSpecie() const")
      << "attempt to access null inertSpecie, probably using wrong boundary condition"
      << abort(FatalError);

    return word(NULL);
}  


void Foam::multiSpeciesTransportModel::write()
{
    mix().write();
    thermo().p().write();
    mt().write();
    forAll(p_, i)
    {
        p_[i].write();
    }
}  


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
