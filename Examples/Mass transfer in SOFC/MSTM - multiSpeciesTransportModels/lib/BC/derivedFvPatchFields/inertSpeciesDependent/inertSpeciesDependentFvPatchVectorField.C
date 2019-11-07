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

#include "inertSpeciesDependentFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inertSpeciesDependentFvPatchVectorField::
inertSpeciesDependentFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    inertFlux_(p.size(), 0.0)
{}


inertSpeciesDependentFvPatchVectorField::
inertSpeciesDependentFvPatchVectorField
(
    const inertSpeciesDependentFvPatchVectorField& pivpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(p, iF),
    inertFlux_(pivpvf.inertFlux_, mapper) 
{
    fvPatchVectorField::operator=(0*patch().nf());
}


inertSpeciesDependentFvPatchVectorField::
inertSpeciesDependentFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    inertFlux_("inertMolarFlux", dict, p.size())
{
    fvPatchVectorField::operator=(0*patch().nf());
}


inertSpeciesDependentFvPatchVectorField::
inertSpeciesDependentFvPatchVectorField
(
    const inertSpeciesDependentFvPatchVectorField& pivpvf
)
:
    fixedValueFvPatchVectorField(pivpvf),
    inertFlux_(pivpvf.inertFlux_)
{}


inertSpeciesDependentFvPatchVectorField::
inertSpeciesDependentFvPatchVectorField
(
    const inertSpeciesDependentFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    inertFlux_(pivpvf.inertFlux_) 
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void inertSpeciesDependentFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const multiSpeciesTransportModel& mstm = 
       this->db().objectRegistry::lookupObject
       <const multiSpeciesTransportModel>("transportProperties");

    const fvPatchField<scalar>& rho =
        patch().lookupPatchField<volScalarField, scalar>("rho");

    const scalarField& A = this->patch().magSf();

    scalarField value = rho * 0;
    label inertIndex = -1;

    forAll(mstm.mix().species(), specieI)
    {
        if(mstm.mix().species()[specieI] != mstm.inertSpecie())
        {
            const fvsPatchField<scalar>& n =
                patch().lookupPatchField<surfaceScalarField, scalar>
                (("n_"+ mstm.mix().species()[specieI]));

            value += n;
        }
        else
        {
            inertIndex = specieI;
        }
    }

    if(mstm.bases() == "mass")
    {
        value += inertFlux_;
    }
    
    if(mstm.bases() == "molar")
    {
        value += inertFlux_ * mstm.speciesThermo()[inertIndex].W();
    }    
  
    value /= (rho * A);
   
    operator== (value * patch().nf());

    fixedValueFvPatchVectorField::updateCoeffs();
}


void inertSpeciesDependentFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    inertFlux_.writeEntry("inertMolarFlux", os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    inertSpeciesDependentFvPatchVectorField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
