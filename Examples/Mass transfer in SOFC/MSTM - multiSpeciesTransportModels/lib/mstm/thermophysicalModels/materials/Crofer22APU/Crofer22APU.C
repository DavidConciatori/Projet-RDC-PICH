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

\*---------------------------------------------------------------------------*/

#include "Crofer22APU.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Crofer22APU, 0);
    addToRunTimeSelectionTable(material, Crofer22APU,);
    addToRunTimeSelectionTable(material, Crofer22APU, Istream);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Crofer22APU::Crofer22APU()
:
    material
    (
        7700,	// Density [kg/m^3]
        660,	// Specific heat capacity [J/(kg.K)]
        24	// Thermal conductivity [W/(m.K)]
    ),
    sigma_
    (
        4.475e5,
	5.632e-1
    )
{}


Foam::Crofer22APU::Crofer22APU
(
    const material& m,
    const SOFCteamFunc4& sigma
)
:
    material(m),
    sigma_(sigma)
{}


Foam::Crofer22APU::Crofer22APU(Istream& is)
:
    material(is),
    sigma_(is)
{}


// ************************************************************************* //
