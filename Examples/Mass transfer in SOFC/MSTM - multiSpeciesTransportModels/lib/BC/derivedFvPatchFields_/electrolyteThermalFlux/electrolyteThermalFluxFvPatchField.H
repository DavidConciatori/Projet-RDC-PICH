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

Class
    Foam::electrolyteThermalFluxFvPatchField

Description
    Thermal boundary condition for electrolyte

SourceFiles
    electrolyteThermalFluxFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef electrolyteThermalFluxFvPatchField_H
#define electrolyteThermalFluxFvPatchField_H

#include "fixedValueFvPatchFields.H"
#include "electrolyte.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class electrolyteThermalFluxFvPatch Declaration
\*---------------------------------------------------------------------------*/

class electrolyteThermalFluxFvPatchField
:
    public fixedValueFvPatchScalarField
{
  
    // Private data
    
        //- Faraday constant
        static const dimensionedScalar F;
      
        //- Heat flux
        scalarField heatFlux_;

public:

    //- Runtime type information
    TypeName("electrolyteThermalFlux");


    // Constructors

        //- Construct from patch and internal field
        electrolyteThermalFluxFvPatchField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        electrolyteThermalFluxFvPatchField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given electrolyteThermalFluxFvPatchField onto a new patch
        electrolyteThermalFluxFvPatchField
        (
            const electrolyteThermalFluxFvPatchField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        electrolyteThermalFluxFvPatchField
        (
            const electrolyteThermalFluxFvPatchField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new electrolyteThermalFluxFvPatchField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        electrolyteThermalFluxFvPatchField
        (
            const electrolyteThermalFluxFvPatchField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new electrolyteThermalFluxFvPatchField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();
        
        //- Return the heat flux
        const scalarField& heatFlux() const 
        {
          return heatFlux_;
        }        
        
        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
