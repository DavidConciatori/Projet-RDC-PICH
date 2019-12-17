/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation 
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "krBazant.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diffusionCoefficientModels
{
defineTypeNameAndDebug(krBazant, 0);

addToRunTimeSelectionTable
(
    diffusionCoefficientModel,
    krBazant,
    dictionary
);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diffusionCoefficientModels::krBazant::krBazant
(
    const word& name,
    const dictionary& transportProperties,
    const volScalarField& S
)
    :
    diffusionCoefficientModel(name, transportProperties,S),

    krBazantCoeffs_(transportProperties.subDict(typeName + "Coeffs")),
    n_
    (
        IOobject
        (
            "n",
            S_.time().timeName(),
            S_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        S.mesh(),
        dimensionedScalar("n",dimless,krBazantCoeffs_.lookupOrDefault<scalar>("n",0))
    ),

    Sc_
    (
        IOobject
        (
            "Sc",
            S_.time().timeName(),
            S_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        S.mesh(),
        dimensionedScalar("Sc",dimless,krBazantCoeffs_.lookupOrDefault<scalar>("Sc",0))
    ),

    alpha0_
    (
        IOobject
        (
            "alpha0",
            S_.time().timeName(),
            S_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        S.mesh(),
        dimensionedScalar("alpha0",dimless,krBazantCoeffs_.lookupOrDefault<scalar>("alpha0",0))
    ),


{

    if (gMin(n_) <= 0)
    {
        FatalErrorIn
            (
                "in krBazant.C"
            )
            << "Diffusion Model coefficient n equal or less than 0" 
                << exit(FatalError);
    }
 
/*    Info << "Brooks and Corey parameters for relative permeability model" << nl << "{" << endl;
    Info << "    n ";
    if (n_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(n_).value() << endl;}
    Info << "    Smax ";
    if (Smax_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(Smax_).value() << endl;}
    Info <<  "    Smin ";
    if (Smin_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(Smin_).value() << endl;}
    Info << "    kramax ";
    if (kramax_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(kramax_).value() << endl;}
    Info << "    krbmax ";
    if (krbmax_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(krbmax_).value() << endl;}
    Info << "} \n" << endl;*/
}

// ************************************************************************* //
