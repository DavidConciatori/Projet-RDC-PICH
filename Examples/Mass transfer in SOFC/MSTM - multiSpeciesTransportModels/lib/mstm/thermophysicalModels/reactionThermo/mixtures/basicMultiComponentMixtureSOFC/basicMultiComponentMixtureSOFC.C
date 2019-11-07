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

#include "basicMultiComponentMixtureSOFC.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicMultiComponentMixtureSOFC::basicMultiComponentMixtureSOFC
(
    const dictionary& thermoDict,
    const wordList& specieNames,
    const fvMesh& mesh
)
:
    species_(specieNames),
    W_
    (
        IOobject
        (
            "W",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("W", dimMass/dimMoles, 0.0)
    ),    
    
    y_(species_.size()),
    
    x_(species_.size()),
    
    n_(species_.size()),
    
    N_(species_.size())
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::basicMultiComponentMixtureSOFC::write()
{
    forAll(y_, i)
    {
        y_[i].write();
        x_[i].write();
        n_[i].write();
        N_[i].write();
    }
} 


// ************************************************************************* //
