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

#include "fixedFluxFvPatchField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedFluxFvPatchField::fixedFluxFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchField<scalar>(p, iF)
{
    this->gradient() = pTraits<scalar>::zero;
}


fixedFluxFvPatchField::fixedFluxFvPatchField
(
    const fixedFluxFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchField<scalar>(ptf, p, iF, mapper)
{}


fixedFluxFvPatchField::fixedFluxFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchField<scalar>(p, iF, dict)
{
    // Set dummy gradient
    this->gradient() = pTraits<scalar>::zero;

    // Read the value entry from the dictionary
    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            Field<scalar>("value", dict, p.size())
        );
    }
    else
    {
        FatalIOErrorIn
        (
            "fluxFvPatchField<scalar>::fluxFvPatchField"
            "("
            "const fvPatch& p,"
            "const DimensionedField<scalar, volMesh>& iF,"
            "const dictionary& dict,"
            "const bool valueRequired"
            ")",
            dict
        )   << "Essential entry 'value' missing"
            << exit(FatalIOError);
    }
}


fixedFluxFvPatchField::fixedFluxFvPatchField
(
    const fixedFluxFvPatchField& ptf
)
:
    fixedGradientFvPatchField<scalar>(ptf)
{}



fixedFluxFvPatchField::fixedFluxFvPatchField
(
    const fixedFluxFvPatchField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchField<scalar>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void fixedFluxFvPatchField::updateCoeffs()
{

    if (this->updated())
    {
        return;
    }

    const multiSpeciesTransportModel& mstm = 
       this->db().objectRegistry::lookupObject
       <const multiSpeciesTransportModel>("transportProperties");

    const word& fieldName =  this->dimensionedInternalField().name();
    const word& speciesName =  fieldName(2,fieldName.size());
    
    const scalarField& A = this->patch().magSf();
       
    const fvPatchField<scalar>& Gamma =
        patch().lookupPatchField<volScalarField, scalar>
        ("Gamma_" + speciesName);
 
    const fvsPatchField<scalar>& phi_P =
        patch().lookupPatchField<surfaceScalarField, scalar>
        ("phi_P_" + speciesName);

    const fvsPatchField<scalar>& phi_N =
        patch().lookupPatchField<surfaceScalarField, scalar>
        ("phi_N_" + speciesName);    
    
    if(mstm.bases() == "mass")
    { 
        const fvsPatchField<scalar>& nAlpha =
            patch().lookupPatchField<surfaceScalarField, scalar>
            ("n_" + speciesName);
      
        this->gradient() = ( (phi_P+phi_N) * *this - nAlpha ) / Gamma / A;
    }
    if(mstm.bases() == "molar")
    { 
        const fvsPatchField<scalar>& NAlpha =
            patch().lookupPatchField<surfaceScalarField, scalar>
            ("N_" + speciesName);
     
        this->gradient() = ( (phi_P+phi_N) * *this - NAlpha ) / Gamma / A; 
    }

    fixedGradientFvPatchField<scalar>::updateCoeffs();
}


void fixedFluxFvPatchField::write(Ostream& os) const
{
    fixedGradientFvPatchField<scalar>::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, fixedFluxFvPatchField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
