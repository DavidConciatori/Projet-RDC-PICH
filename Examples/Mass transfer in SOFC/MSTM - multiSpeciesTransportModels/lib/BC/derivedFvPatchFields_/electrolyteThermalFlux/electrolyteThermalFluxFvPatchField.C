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

#include "electrolyteThermalFluxFvPatchField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //  
  
const  dimensionedScalar Foam::electrolyteThermalFluxFvPatchField::F
(
    "F", dimensionSet(0, 0, 1, 0, -1, 1, 0), 96485.3415 * 1000
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

electrolyteThermalFluxFvPatchField::electrolyteThermalFluxFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    heatFlux_(p.size(), 0.0)
{}


electrolyteThermalFluxFvPatchField::electrolyteThermalFluxFvPatchField
(
    const electrolyteThermalFluxFvPatchField& pivpvf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(p, iF),
    heatFlux_(pivpvf.heatFlux_, mapper)
{
    //fvPatchScalarField::operator=(refValue_*patch().nf());
}


electrolyteThermalFluxFvPatchField::electrolyteThermalFluxFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    heatFlux_("heatFlux", dict, p.size())
{
    //fvPatchScalarField::operator=(refValue_*patch().nf());
}


electrolyteThermalFluxFvPatchField::electrolyteThermalFluxFvPatchField
(
    const electrolyteThermalFluxFvPatchField& pivpvf
)
:
    fixedValueFvPatchScalarField(pivpvf),
    heatFlux_(pivpvf.heatFlux_)
{}


electrolyteThermalFluxFvPatchField::electrolyteThermalFluxFvPatchField
(
    const electrolyteThermalFluxFvPatchField& pivpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(pivpvf, iF),
    heatFlux_(pivpvf.heatFlux_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::electrolyteThermalFluxFvPatchField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    heatFlux_.autoMap(m);
}


void Foam::electrolyteThermalFluxFvPatchField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const electrolyteThermalFluxFvPatchField& hfptf =
        refCast<const electrolyteThermalFluxFvPatchField>(ptf);

    heatFlux_.rmap(hfptf.heatFlux_, addr);
}


void electrolyteThermalFluxFvPatchField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const label index = patch().index();
    
    const SOFC::electrolyte& electr = 
       this->db().objectRegistry::lookupObject<const SOFC::electrolyte>("elementProperties");
 
    const scalarField& A = this->patch().magSf();

    const scalarField& q_homic = electr.q().internalField();
    const scalarField& q_act_and_flux = electr.q().boundaryField()[index];

    heatFlux_ = (0.5 * q_homic + q_act_and_flux) * A;	
    
    fixedValueFvPatchScalarField::updateCoeffs();
}


void electrolyteThermalFluxFvPatchField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    heatFlux_.writeEntry("heatFlux", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    electrolyteThermalFluxFvPatchField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
