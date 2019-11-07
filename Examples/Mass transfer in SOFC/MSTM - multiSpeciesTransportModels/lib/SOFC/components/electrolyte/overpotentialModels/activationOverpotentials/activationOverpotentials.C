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

#include "activationOverpotentials.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  
namespace SOFC
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

    defineTypeNameAndDebug(activationOverpotentials, 0);

    const  scalar activationOverpotentials::tol_(1.0e-8);
    const  scalar activationOverpotentials::maxIter_(1000);
}

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SOFC::activationOverpotentials::activationOverpotentials
(
    Foam::wordList reactions,
    const Foam::dictionary& dic
)
:
    reactions_(reactions)
{
    etaAct_.setSize(reactions.size());
    forAll(etaAct_, reactionI)
    {
        etaAct_.set
        (
            reactionI,
            new activationOverpotential(dic.subDict(reactions[reactionI]))
        );
    }
}
  
  
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::SOFC::activationOverpotentials::update(Foam::scalarField& T)
{
    forAll(etaAct_, reactionI)
    {
      etaAct_[reactionI].update(T);
    }
}
  
  
void Foam::SOFC::activationOverpotentials::initX
(
    Foam::PtrList<Foam::volScalarField>& x
)
{   
    forAll(etaAct_, reactionI)
    {    
        forAll(etaAct_[reactionI].species(), speciesI)
        {
            etaAct_[reactionI].readX(speciesI,
                findSpecies(etaAct_[reactionI].species()[speciesI], x));
        }
    }
}


Foam::scalarField& Foam::SOFC::activationOverpotentials::findSpecies
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
