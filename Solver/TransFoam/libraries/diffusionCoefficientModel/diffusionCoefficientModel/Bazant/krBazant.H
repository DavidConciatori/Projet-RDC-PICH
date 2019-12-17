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

Class
    Foam::diffusionCoefficientModels::krBazant

Description
    Standard Brooks and Corey relative permeability model.

SourceFiles
    krBazant.C

\*---------------------------------------------------------------------------*/

#ifndef krBazant_H
#define krBazant_H

#include "diffusionCoefficientModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diffusionCoefficientModels
{

/*---------------------------------------------------------------------------*\
    Class krBazant Declaration
    \*---------------------------------------------------------------------------*/

class krBazant
:
    public diffusionCoefficientModel
{
    //- Private data

    dictionary krBazantCoeffs_;
    volScalarField n_;
    volScalarField Sc;
    volScalarField alpha0;
    //- effective saturation of phase b
    // volScalarField Se_;

    //- end points
    // volScalarField kramax_;
    // volScalarField krbmax_;

        
public:

    //- Runtime type information
    TypeName("Bazant");

    // Constructors

    //- Construct from components
    krBazant
    (
        const word& name,
        const dictionary& transportProperties,
        const volScalarField& S
    );

    //- Destructor
    ~krBazant()
        {}

    // Member Functions

    //- Correct the diffusion coefficient
    void correct()
        {
            kr_ = alpha0_+((scalar(1)-alpha0_)/(scalar(1)+(scalar(4)*(scalar(1)-S_)/(scalar(1)-Sc_)))); //compute the diffusion multiplizer using Bazant's Model

            dkrdS_ = scalar(4)*(scalar(1)-alpha0_)/(scalar(1)-Sc_)*(scalar(1)/pow((1-scalar(4)*(scalar(1)-S_)/(scalar(1)-Sc_)),2));
  
        }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diffusionCoefficientModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
