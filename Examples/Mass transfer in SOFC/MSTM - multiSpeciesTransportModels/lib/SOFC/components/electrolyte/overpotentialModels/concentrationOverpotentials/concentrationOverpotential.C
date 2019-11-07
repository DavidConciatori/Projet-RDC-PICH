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

#include "concentrationOverpotential.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  
namespace SOFC
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

    defineTypeNameAndDebug(concentrationOverpotential, 0);

    const  scalar concentrationOverpotential::R = 8314.51;
    const  scalar concentrationOverpotential::F = 96485.3415 * 1000;
}

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SOFC::concentrationOverpotential::concentrationOverpotential
(
    const dictionary& dic
)
:
    species_(dic.lookup("species")),
    
    n_(readScalar(dic.lookup("n"))),
    
    stechCoeff_(dic.lookup("stechCoeff"))
{
    x_.setSize(species_.size());
    xRef_.setSize(species_.size());
}
  
  
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
  
Foam::scalar Foam::SOFC::concentrationOverpotential::NernstVoltage
(
    Foam::scalar Tref
)
{
    // Computational Fluid-Dynamics Modeling of H2-fed Solid Oxide Fuel Cells 
    // María García-Camprubí, Hrvoje Jasak and Norberto Fueyo
    scalar E0 = (1.271 - 2.731e-4*Tref) / 2;
    scalar Nernst = E0 ;
    
    forAll(species_, speciesI)
    {
        Nernst -= R*Tref/(n_*F)*log(pow(xRef_[speciesI],stechCoeff_[speciesI]));
    }
    
    return Nernst;
}


Foam::scalarField Foam::SOFC::concentrationOverpotential::eta
(
    const Foam::scalarField& T
) const
{  
    scalarField etaConc = 0*T;        
    forAll(species(), speciesI)
    {
        etaConc -= R*T/(n_*F)*log(pow(xRef_[speciesI]/x_[speciesI],stechCoeff_[speciesI]));
    }
    
    return etaConc;      
}


void Foam::SOFC::concentrationOverpotential::readX
(
    int index,
    Foam::scalarField& x,
    Foam::scalar xRef
)
{
    x_.set(index, x);
    xRef_[index] = xRef;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
