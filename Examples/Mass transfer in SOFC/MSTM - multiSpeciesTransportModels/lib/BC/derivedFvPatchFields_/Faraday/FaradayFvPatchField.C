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

#include "FaradayFvPatchField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //    
  
const scalar FaradayFvPatchField::F = 96485.3415 * 1000;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

FaradayFvPatchField::FaradayFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    molarFlux_(p.size(), 0.0)
{}


FaradayFvPatchField::FaradayFvPatchField
(
    const FaradayFvPatchField& pivpvf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(p, iF),
    molarFlux_(pivpvf.molarFlux_, mapper)
{
    //fvPatchScalarField::operator=(refValue_*patch().nf());
}


FaradayFvPatchField::FaradayFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    molarFlux_("molarFlux", dict, p.size())
{
    //fvPatchScalarField::operator=(refValue_*patch().nf());
}


FaradayFvPatchField::FaradayFvPatchField
(
    const FaradayFvPatchField& pivpvf
)
:
    fixedValueFvPatchScalarField(pivpvf),
    molarFlux_(pivpvf.molarFlux_)
{}


FaradayFvPatchField::FaradayFvPatchField
(
    const FaradayFvPatchField& pivpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(pivpvf, iF),
    molarFlux_(pivpvf.molarFlux_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::FaradayFvPatchField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    molarFlux_.autoMap(m);
}


void Foam::FaradayFvPatchField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const FaradayFvPatchField& hfptf =
        refCast<const FaradayFvPatchField>(ptf);

    molarFlux_.rmap(hfptf.molarFlux_, addr);
}


void FaradayFvPatchField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const SOFC::electrolyte& electr = 
       this->db().objectRegistry::lookupObject
       <const SOFC::electrolyte>("elementProperties");

    const scalarField& T =
        patch().lookupPatchField<volScalarField, scalar>("T");

    const word& fieldName =  this->dimensionedInternalField().name();
    const word& name =  fieldName(2,fieldName.size());
    
    const scalarField& A = this->patch().magSf();
     
    const scalarField& eta =
        patch().lookupPatchField<volScalarField, scalar>("eta");    

        
    forAll(electr.anodicSpecies(), speciesI)
    {
        if(electr.anodicSpecies()[speciesI] == name)
        {
            forAll(electr.etaConcAnode().etaConc_, reactionI)
            {
                const wordList& reactionSpecie =
                    electr.etaConcAnode().etaConc_[reactionI].species();
                    
                forAll(reactionSpecie, specieI)
                {
                    if (reactionSpecie[specieI]==name)
                    {
                        scalarField i =
                            electr.etaActAnode().etaAct_[reactionI].i(T, eta);
                            
                        const scalar& n_ =
                            electr.etaConcAnode().etaConc_[reactionI].n();
                            
                        const scalarList& stechCoeff_ =
                            electr.etaConcAnode()
                            .etaConc_[reactionI].stechCoeff();
                            
                        molarFlux_ = stechCoeff_[specieI]*i/(n_*F)*A;
                    }
                }
            }        
        }
    }
    
    forAll(electr.cathodicSpecies(), speciesI)
    {
        if(electr.cathodicSpecies()[speciesI] == name)
        {
            forAll(electr.etaConcCathode().etaConc_, reactionI)
            {
                const wordList& reactionSpecie =
                    electr.etaConcCathode().etaConc_[reactionI].species();
                    
                forAll(reactionSpecie, specieI)
                {
                    if (reactionSpecie[specieI]==name)
                    {
                        scalarField i =
                            electr.etaActCathode()
                            .etaAct_[reactionI].i(T, eta);
                            
                        const scalar& n_ =
                            electr.etaConcCathode().etaConc_[reactionI].n();
                            
                        const scalarList& stechCoeff_ =
                            electr.etaConcCathode()
                            .etaConc_[reactionI].stechCoeff();
                            
                        molarFlux_ = stechCoeff_[specieI]*i/(n_*F)*A;
                    }
                }
            }
        }
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void FaradayFvPatchField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    molarFlux_.writeEntry("molarFlux", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    FaradayFvPatchField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
