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

#include "activationOverpotential.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  
namespace SOFC
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

    defineTypeNameAndDebug(activationOverpotential, 0);    
}

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SOFC::activationOverpotential::activationOverpotential
(
    const dictionary& dic
)
:
    species_(dic.lookup("species")),
    
    i0model(dic),
    
    BV_(scalarList(dic.lookup("theta")))
{
    x_.setSize(species_.size());
}
  
  
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::SOFC::activationOverpotential::update(Foam::scalarField& T)
{
    if(species_.size() == 1)
    {
        i0_ = i0model.calc(T, x_[0]);
    }
    if(species_.size() == 2)
    {
        i0_ = i0model.calc(T, x_[0], x_[1]);
    }
}
  
  
Foam::scalarField Foam::SOFC::activationOverpotential::i
(
    const Foam::scalarField& T,
    const Foam::scalarField& eta
) const
{  
     return BV_.i(i0_, T, eta);       
}


void Foam::SOFC::activationOverpotential::readX
(
    int index,
    Foam::scalarField& x
)
{
    x_.set(index, x);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
