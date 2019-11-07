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

#include "ohmicOverpotential.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  
namespace SOFC
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

    defineTypeNameAndDebug(ohmicOverpotential, 0);
}

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SOFC::ohmicOverpotential::ohmicOverpotential
(
    Foam::SOFC::solidComponent& elec,
    Foam::SOFC::solidComponent& anode,
    Foam::SOFC::solidComponent& cathode
)
:
    anode_(anode),
  
    electrolyte_(elec),
  
    cathode_(cathode)
{
    IOdictionary electrolyteDictionary
    (
        IOobject
        (
            "electrolyte",
            elec.mesh().time().constant(),
            elec.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
      
    scalar dAThickness = readScalar(electrolyteDictionary.lookup("dAThickness"));
    
    thickness_.setSize(3);
    thickness_[0] = thickness(anode);
    thickness_[1] = thickness(electrolyte_) + dAThickness;
    thickness_[2] = thickness(cathode);
}
  

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::SOFC::ohmicOverpotential::update(Foam::scalarField& T)
{
    totalSigma_.setSize(T.size());
    forAll(totalSigma_, i)
    {
      totalSigma_[i] = thickness_[0]/anode_.solid()->sigma(T[i]);
      totalSigma_[i] += thickness_[1]/electrolyte_.solid()->sigma(T[i]);
      totalSigma_[i] += thickness_[2]/cathode_.solid()->sigma(T[i]);  
    }
}


Foam::scalar Foam::SOFC::ohmicOverpotential::thickness
(
    Foam::SOFC::solidComponent& component
)
{
    const fvMesh& mesh = component.mesh();
    
    scalarList deltaThick(3);
    
    deltaThick[0] = max(mesh.Cf().component(vector::X)).value()-min(mesh.Cf().component(vector::X)).value();
    deltaThick[1] = max(mesh.Cf().component(vector::Y)).value()-min(mesh.Cf().component(vector::Y)).value();
    deltaThick[2] = max(mesh.Cf().component(vector::Z)).value()-min(mesh.Cf().component(vector::Z)).value();
     
    scalar thikness = 100;
    
    forAll(deltaThick, i)
    {
      if (deltaThick[i] > 1e-10)
      {
        thikness = min(thikness,deltaThick[i]);
      }
    }
    
    return thikness;    
  }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
