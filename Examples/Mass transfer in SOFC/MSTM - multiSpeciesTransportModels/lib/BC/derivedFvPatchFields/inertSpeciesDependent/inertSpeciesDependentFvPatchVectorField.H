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
    Foam::inertSpeciesDependentFvPatchVectorField

Description
    Describes a surface normal vector boundary condition by its magnitude.
    Note: The value is positive for outward-pointing vectors

SourceFiles
    inertSpeciesDependentFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef inertSpeciesDependentFvPatchVectorField_H
#define inertSpeciesDependentFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "multiSpeciesTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class inertSpeciesDependentFvPatch Declaration
\*---------------------------------------------------------------------------*/

class inertSpeciesDependentFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data
    
        //- Value of the mass/molar (depeding on the multispecies transport
        //  model) flux of the not directly solved specie
        scalarField inertFlux_;

public:

    //- Runtime type information
    TypeName("inertSpeciesDependent");


    // Constructors

        //- Construct from patch and internal field
        inertSpeciesDependentFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        inertSpeciesDependentFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given inertSpeciesDependentFvPatchVectorField onto a new patch
        inertSpeciesDependentFvPatchVectorField
        (
            const inertSpeciesDependentFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        inertSpeciesDependentFvPatchVectorField
        (
            const inertSpeciesDependentFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new inertSpeciesDependentFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        inertSpeciesDependentFvPatchVectorField
        (
            const inertSpeciesDependentFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new inertSpeciesDependentFvPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Mapping functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
