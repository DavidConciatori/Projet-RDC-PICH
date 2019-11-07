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

#include "component.H"
#include "FaradayFvPatchField.H"
#include "electrolyteThermalFluxFvPatchField.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::SOFC::component::check
(
    label& fromIndex,
    label& toIntex,
    word& fromPatchName,
    word& toPatchName
)
{
    labelList index(2);  index[0] = fromIndex;      index[1] = toIntex;
    wordList patch(2);   patch[0] = fromPatchName;  patch[1] = toPatchName;

    forAll(index, I)
    {
        if(index[I] < 0)
        {
            Info << "No boundary called " << patch[I] << "\n";
            Foam::FatalError.exit();
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SOFC::component::component
(
    const Foam::word name,
    Foam::regionPropertiesSOFC& regions
)
:
    IOdictionary
    (
        IOobject
        (
            "elementProperties",
            regions.mesh(name).time().constant(),
            regions.mesh(name),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(regions.mesh(name)),
    name_(name),
    activated_(1)
{
    Info << "\n\nCreatin the " << name << "\n\n";
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
