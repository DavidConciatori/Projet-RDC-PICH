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

#include "hsCombustionThermoSOFC.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::hsCombustionThermoSOFC> Foam::hsCombustionThermoSOFC::New
(
    const fvMesh& mesh
)
{
    word hsCombustionThermoSOFCTypeName;

    // Enclose the creation of the thermophysicalProperties to ensure it is
    // deleted before the turbulenceModel is created otherwise the dictionary
    // is entered in the database twice
    {
        IOdictionary thermoDict
        (
            IOobject
            (
                "thermophysicalProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        thermoDict.lookup("thermoType") >> hsCombustionThermoSOFCTypeName;
    }

    Info<< "Selecting thermodynamics package " << hsCombustionThermoSOFCTypeName
        << endl;

    fvMeshConstructorTable::iterator cstrIter =
        fvMeshConstructorTablePtr_->find(hsCombustionThermoSOFCTypeName);

    if (cstrIter == fvMeshConstructorTablePtr_->end())
    {
        FatalErrorIn("hsCombustionThermoSOFC::New(const fvMesh&)")
            << "Unknown hsCombustionThermoSOFC type "
            << hsCombustionThermoSOFCTypeName << nl << nl
            << "Valid hsCombustionThermoSOFC types are:" << nl
            << fvMeshConstructorTablePtr_->toc() << nl
            << exit(FatalError);
    }

    return autoPtr<hsCombustionThermoSOFC>(cstrIter()(mesh));
}


Foam::autoPtr<Foam::hsCombustionThermoSOFC> Foam::hsCombustionThermoSOFC::New
(
    const fvMesh& mesh,
    const word& thermoTypeName
)
{
    word hsCombustionThermoSOFCTypeName(thermoTypeName);

    Info<< "Selecting thermodynamics package " << hsCombustionThermoSOFCTypeName
        << endl;

    fvMeshConstructorTable::iterator cstrIter =
        fvMeshConstructorTablePtr_->find(hsCombustionThermoSOFCTypeName);

    if (cstrIter == fvMeshConstructorTablePtr_->end())
    {
        FatalErrorIn("hsCombustionThermoSOFC::New(const fvMesh&)")
            << "Unknown hsCombustionThermoSOFC type "
            << hsCombustionThermoSOFCTypeName << nl << nl
            << "Valid hsCombustionThermoSOFC types are:" << nl
            << fvMeshConstructorTablePtr_->toc() << nl
            << exit(FatalError);
    }

    return autoPtr<hsCombustionThermoSOFC>(cstrIter()(mesh));
}


Foam::autoPtr<Foam::hsCombustionThermoSOFC>
Foam::hsCombustionThermoSOFC::NewType
(
    const fvMesh& mesh,
    const word& thermoType
)
{
    word hsCombustionThermoSOFCTypeName;

    // Enclose the creation of the thermophysicalProperties to ensure it is
    // deleted before the turbulenceModel is created otherwise the dictionary
    // is entered in the database twice
    {
        IOdictionary thermoDict
        (
            IOobject
            (
                "thermophysicalProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        thermoDict.lookup("thermoType") >> hsCombustionThermoSOFCTypeName;

        if (hsCombustionThermoSOFCTypeName.find(thermoType) == string::npos)
        {
            wordList allModels = fvMeshConstructorTablePtr_->toc();
            DynamicList<word> validModels;
            forAll(allModels, i)
            {
                if (allModels[i].find(thermoType) != string::npos)
                {
                    validModels.append(allModels[i]);
                }
            }

            FatalErrorIn
            (
                "autoPtr<hsCombustionThermoSOFC> hsCombustionThermoSOFC::NewType"
                "("
                    "const fvMesh&, "
                    "const word&"
                ")"
            )   << "Inconsistent thermo package selected:" << nl << nl
                << hsCombustionThermoSOFCTypeName << nl<< nl << "Please select a "
                << "thermo package based on " << thermoType
                << ". Valid options include:" << nl << validModels << nl
                << exit(FatalError);
        }
    }

    Info<< "Selecting thermodynamics package " 
        << hsCombustionThermoSOFCTypeName << endl;

    fvMeshConstructorTable::iterator cstrIter =
        fvMeshConstructorTablePtr_->find(hsCombustionThermoSOFCTypeName);

    if (cstrIter == fvMeshConstructorTablePtr_->end())
    {
        FatalErrorIn("hsCombustionThermoSOFC::New(const fvMesh&)")
            << "Unknown hsCombustionThermoSOFC type "
            << hsCombustionThermoSOFCTypeName << nl << nl
            << "Valid hsCombustionThermoSOFC types are:" << nl
            << fvMeshConstructorTablePtr_->toc() << nl
            << exit(FatalError);
    }

    return autoPtr<hsCombustionThermoSOFC>(cstrIter()(mesh));
}


// ************************************************************************* //
