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

Namespace
    Foam::capillarityModel

Description
    A namespace for various capillarity Effect implementations.

Class
    Foam::capillarityModel

Description
    An abstract base class for incompressible capillarityModels.

SourceFiles
    diffusionCoefficientModel.C

\*---------------------------------------------------------------------------*/

#ifndef diffusionCoefficientModel_H
#define diffusionCoefficientModel_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class diffusionCoefficientModel Declaration
\*---------------------------------------------------------------------------*/

class diffusionCoefficientModel
{

protected:

    // Protected data
    word name_;
    dictionary transportProperties_;

    // reference to the saturation field 
    const volScalarField& S_;

    // capillary pressure, derivative and capillary capacity
    volScalarField kr_;
    volScalarField dkrdS_;

    // Private Member Functions

    //- Disallow copy construct
    diffusionCoefficientModel(const diffusionCoefficientModel&);

    //- Disallow default bitwise assignment
    void operator=(const diffusionCoefficientModel&);

public:

    //- Runtime type information
    TypeName("diffusionCoefficientModel");


    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        diffusionCoefficientModel,
        dictionary,
        (
            const word& name,
            const dictionary& transportProperties,
            const volScalarField& S
        ),
        (name, transportProperties, S)
    );


    // Selectors

    //- Return a reference to the selected capillarity model
    static autoPtr<diffusionCoefficientModel> New
    (
        const word& name,
        const dictionary& transportProperties,
        const volScalarField& S
    );


    // Constructors

    //- Construct from components
    diffusionCoefficientModel
    (
        const word& name,
        const dictionary& transportProperties,
        const volScalarField& S
    );


    //- Destructor
    virtual ~diffusionCoefficientModel() {}


    // Member Functions

    //- Return the phase transport properties dictionary
    const dictionary& transportProperties() const
        {
            return transportProperties_;
        }
      
    //- Return the Bazant coefficient as a function of S 
    tmp<volScalarField> kr() const
        {
            return kr_;
        }

    //- Return the derivative of Bazant coefficient 
    tmp<volScalarField> dkrdS() const
        {
            return dkrdS_;
        }

    //- Correct the Bazant coefficient
    virtual void correct() = 0; 


};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
