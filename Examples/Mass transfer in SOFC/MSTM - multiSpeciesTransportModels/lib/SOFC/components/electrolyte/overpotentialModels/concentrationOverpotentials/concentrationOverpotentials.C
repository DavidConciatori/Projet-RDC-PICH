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

#include "concentrationOverpotentials.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  
namespace SOFC
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

    defineTypeNameAndDebug(concentrationOverpotentials, 0);
}

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SOFC::concentrationOverpotentials::concentrationOverpotentials
(
    Foam::wordList reactions,
    const Foam::dictionary& dic
)
:
    reactions_(reactions)
{
    etaConc_.setSize(reactions.size());
    forAll(etaConc_, reactionI)
    {
        etaConc_.set
        (
            reactionI,
            new concentrationOverpotential(dic.subDict(reactions[reactionI]))
        );
    }
}
  

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
  
Foam::scalar Foam::SOFC::concentrationOverpotentials::NernstVoltage
(
    Foam::scalar Tref
)
{
    scalar Nernst = 0;
    forAll(etaConc_, reactionI)
    {
      Nernst += etaConc_[reactionI].NernstVoltage(Tref);
    }
    
    return Nernst;
}


Foam::scalarField Foam::SOFC::concentrationOverpotentials::eta
(
    Foam::scalarField& T
)
{
    scalarField etaConc = 0*T;
    forAll(etaConc_, reactionI)
    {
        etaConc += etaConc_[reactionI].eta(T);
    }
    
    return etaConc;
}


void Foam::SOFC::concentrationOverpotentials::initX
(
    Foam::PtrList<Foam::volScalarField>& x,
    Foam::SOFC::channel& ch
)
{   
    label inletIndex = ch.mesh().boundaryMesh().findPatchID("inlet");
    scalarField A = ch.mesh().magSf().boundaryField()[inletIndex];
    
    forAll(etaConc_, reactionI)
    {    
        forAll(etaConc_[reactionI].species(), speciesI)
        {
            scalar inletValue = gSum(A*ch.mstm().mix().x(etaConc_[reactionI]
                .species()[speciesI]).boundaryField()[inletIndex])/gSum(A);

            etaConc_[reactionI].readX(speciesI, findSpecies(etaConc_[reactionI]
                .species()[speciesI], x), inletValue);
        }
    }
}


Foam::scalarField& Foam::SOFC::concentrationOverpotentials::findSpecies
(
    Foam::word name,
    Foam::PtrList<Foam::volScalarField>& x
)
{ 
    forAll (x, speciesI)
    {
        if (x[speciesI].name() == ("x_"+name))
        {
            return (x[speciesI].internalField());
        }
    }
    
    Info << "\nFatal Error!!!\n" <<
    "The species " << name << " is not present\n" <<"(electrolyte)\n";
    Foam::FatalError.exit();
    
    return (x[0].internalField());
  }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
