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

#include "electrolyte.H"

#include "specieThermo.H"
#include "NASAThermo.H"
#include "perfectGas.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace SOFC
{
  
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

    defineTypeNameAndDebug(electrolyte, 0);      

    const  scalar electrolyte::iToll_ = 1e-2;
    const  scalar electrolyte::VToll_ = 1e-8;
}

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SOFC::electrolyte::electrolyte
(
    const Foam::word name,
    Foam::regionPropertiesSOFC& regions,
    const Foam::dictionary& dic,
    Foam::SOFC::solidComponent& anode,
    Foam::SOFC::solidComponent& cathode,
    Foam::SOFC::channel& anodeChannel,
    Foam::SOFC::channel& cathodeChannel
)
:
    component(name, regions), 
    
    solidComponent(name, regions),
  
    anodeChannel_(anodeChannel),
    
    cathodeChannel_(cathodeChannel),
  
    anodeReactions_
    (
        wordList(dic.subDict("electroChemicalReactions").lookup("anode"))
    ),          
    
    cathodeReactions_
    (
        wordList(dic.subDict("electroChemicalReactions").lookup("cathode"))
    ),
  
    anodicSpecies_(wordList(dic.lookup("anodicSpecies"))),          
    
    cathodicSpecies_(wordList(dic.lookup("cathodicSpecies"))),
  
    anodePatchIndex_
    (
        mesh().boundaryMesh().findPatchID("electrolyte_to_anode")
    ),
    
    cathodePatchIndex_
    (
        mesh().boundaryMesh().findPatchID("electrolyte_to_cathode")
    ),
  
    etaActAnode_(anodeReactions_, dic),
    
    etaActCathode_(cathodeReactions_, dic),

    etaConcAnode_(anodeReactions_, dic),
    
    etaConcCathode_(cathodeReactions_, dic),
  
    etaOhm_(*this, anode, cathode),
  
    T_
    (
        IOobject
        (
            "T",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh()
    ),
  
    q_
    (
        IOobject
        (
            "q",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("q", dimensionSet(1, 0, -3, 0, 0, 0, 0), 0)
    ),
  
    eta_
    (
        IOobject
        (
            "eta",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("eta", dimensionSet(0, -2, 0, 0, 0, 1, 0), 0)
    ),
  
    report_("report.txt")  
{
    speciesName_ = anodicSpecies_; 
    speciesName_.append(cathodicSpecies_);
    
    x_.setSize(speciesName_.size());
    N_.setSize(speciesName_.size());
    
    forAll(speciesName_, speciesI)
    {  
        x_.set
        (
           speciesI, new volScalarField
            (
                IOobject
                (
                    "x_" + speciesName_[speciesI],
                    T_.time().timeName(),
                    T_.db(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                T_.mesh()
            )
        ); 
        
        N_.set
        (
            speciesI, new surfaceScalarField
            (
              IOobject
              (
                  "N_" + speciesName_[speciesI],
                  T_.time().timeName(),
                  T_.db(),
                  IOobject::MUST_READ,
                  IOobject::AUTO_WRITE
              ),
              T_.mesh()
            )
        );       
    }
    
    etaActAnode_.initX(x_);
    etaActCathode_.initX(x_);
    
    etaConcAnode_.initX(x_, anodeChannel);
    etaConcCathode_.initX(x_, cathodeChannel);
    
    scalar Tref(readScalar(dic.lookup("Tref")));
    NernstVoltage_ = etaConcAnode_.NernstVoltage(Tref) + etaConcCathode_.NernstVoltage(Tref);
    Eocv_ = readScalar(dic.lookup("Eocv"));
    iLim_ = readScalar(dic.lookup("iLim"));
    
    sy_.setSize(x_.size());  

    Vrev_ = Eocv_;
    
    V_ = Vrev_;
    
    //setS();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::volScalarField& Foam::SOFC::electrolyte::x(Foam::word name)
{ 
    forAll (x_, speciesI)
    {
        if (x_[speciesI].name() == ("x_"+name))
        {
            return (x_[speciesI]);
        }
    }
    
    Info << "\nFatal Error!!!\n" <<
    "The species " << name << " is not present\n" <<"(electrolyte)\n";
    Foam::FatalError.exit();
    
    return (x_[0]);
}


Foam::surfaceScalarField& Foam::SOFC::electrolyte::N(Foam::word name)
{ 
    forAll (N_, speciesI)
    {
        if (N_[speciesI].name() == ("N_"+name))
        {
            return (N_[speciesI]);
        }
    }
    
    Info << "\nFatal Error!!!\n" <<
    "The species " << name << " is not present\n" <<"(electrolyte)\n";
    Foam::FatalError.exit();
    
    return (N_[0]);
}
  
  
const Foam::volScalarField& Foam::SOFC::electrolyte::sy(Foam::word name) const
{
    forAll (sy_, speciesI)
    {
      if (sy_[speciesI].name() == ("sy_"+name))
      {
        return (sy_[speciesI]);
      }
    }
    
    Info << "\nFatal Error!!!\n" <<
    "The species " << name << " is not present\n" <<"(electrolyte)\n";
    Foam::FatalError.exit();
    
    return (sy_[0]);
}
  
  
void Foam::SOFC::electrolyte::updateTemperature
(
    Foam::SOFC::fluidComponent& anode,
    Foam::SOFC::fluidComponent& cathode
)
{
    mapBoundaryField<volScalarField>(T_, "electrolyte_to_anode");
    mapBoundaryField<volScalarField>(T_, "electrolyte_to_cathode");
    etaOhm_.update(T_.internalField());
}
  
  
void Foam::SOFC::electrolyte::updateMolarFraction
(
    Foam::SOFC::fluidComponent& anode,
    Foam::SOFC::fluidComponent& cathode
)
{  
    forAll(anodicSpecies_, specieI)
    {
        word& name = anodicSpecies_[specieI];
        mapBoundaryField<volScalarField>(x(name), "electrolyte_to_anode");
        scalarField& internalField = x(name).internalField();
        internalField = x(name).boundaryField()[anodePatchIndex_];
    }
    etaConcAnode_.initX(x_, anodeChannel_);
    etaActAnode_.initX(x_);
    etaActAnode_.update(T_.internalField());
    
    forAll(cathodicSpecies_, specieI)
    {
        word& name = cathodicSpecies_[specieI];
        mapBoundaryField<volScalarField>(x(name), "electrolyte_to_cathode");
        scalarField& internalField = x(name).internalField();
        internalField = x(name).boundaryField()[cathodePatchIndex_];
    }
    etaConcCathode_.initX(x_, cathodeChannel_);
    etaActCathode_.initX(x_);
    etaActCathode_.update(T_.internalField());
}


void Foam::SOFC::electrolyte::updateCurrentDensity(Foam::scalar V)
{
    scalarField& T = T_.internalField();
    scalarField etaTOT =
        Vrev_ - V - etaConcAnode_.eta(T) - etaConcCathode_.eta(T);
    
    forAll (i_.internalField(), faceI)
    {
        scalar& i = i_.internalField()[faceI];

        scalar etaTOTnew = 0, etaRes = 10;     
      
        scalar& etaA = eta_.boundaryField()[anodePatchIndex_][faceI];
        scalar& etaC = eta_.boundaryField()[cathodePatchIndex_][faceI];   
        scalar& etaHom = eta_.internalField()[faceI];      
      
        if(etaTOT[faceI] < 0)
        {
            i = 0; etaRes = 0;
            etaA = 0; etaC = 0; etaHom = 0;
        }
     
        while(mag(etaRes) > VToll_)
        {
            if(i == 0) i = 100;
            
            etaA = etaTOT[faceI]/5;
            etaC = etaA * 4;
                  
            etaA = etaActAnode_.eta(T[faceI], i, faceI,etaA);
            etaC = etaActCathode_.eta(T[faceI], i, faceI, etaC);
            etaHom = etaOhm_.eta(i, faceI);  

            etaTOTnew =  (etaA + etaC) + etaHom;
            etaRes = etaTOT[faceI]-etaTOTnew;
            
            i = etaTOT[faceI]/etaTOTnew * i;
        }
    }
}


// WARNING!!!!!!! Must be checked!!!!!
void Foam::SOFC::electrolyte::setS() 
{   
    IOdictionary thermophysicalProperties
    (
        IOobject
        (
            "thermophysicalProperties",
            T_.mesh().time().constant(),
            T_.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    forAll(x_, speciesI)
    {
        word xName = x_[speciesI].name();
        word name = "";
        for(label i=2; i < label(xName.size()); i++)
        {
          name += xName[i];
        }

        specieThermo<NASAThermo<perfectGas> >
            sT(thermophysicalProperties.lookup(name));
        
        const fvMesh& mesh = this->T_.mesh();

        tmp<volScalarField> tES
        (
            new volScalarField
            (
                IOobject
                (
                    "ES",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimEnergy/dimMass
            )
        );   volScalarField& ES = tES();

        forAll(this->T_, celli)
        {
            ES[celli] = sT.G(T_[celli]);
        }

        forAll(this->T_.boundaryField(), patchi)
        {
            const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
            fvPatchScalarField& pES = ES.boundaryField()[patchi];
            forAll(pT, facei)
            {
                pES[facei] = sT.G(pT[facei]);
            }
        }
    
        sy_.set(speciesI,tES); 
        sy_[speciesI].rename("sy_" + name);
    }
}


void Foam::SOFC::electrolyte::updateVoltage(Foam::scalar I)
{       
    Vrev_ = Eocv_ + I / iLim_ * (NernstVoltage_ - Eocv_);  
        
    //Vrev_ = NernstVoltage_ - min(X("H2")) / xH2_in * (NernstVoltage_-Eocv_);  
    //Vrev_ = NernstVoltage_;
    //Vrev_ = Eocv_;
    
    scalar iRes = 10;
    
    scalar increment = 0.001;
    
    label updateIncrement = 0;
    label positive = 0;
    label negative = 0;
        
    while ( mag(iRes) > iToll_ ) 
    {
        if(iRes>0)
        {
            V_ -= increment;
            positive = 1;
        }
        else
        {
            V_ += increment;
            negative = 1;
        }
        if(positive*negative)
        {
            updateIncrement = 1;
            positive = 0;
            negative = 0;
        }
        if(updateIncrement)
        {
            increment = increment/5;
            updateIncrement = 0;
        }
        
        updateCurrentDensity(V_);
        
        iRes = I-currentDensity().value();
    }
}


void Foam::SOFC::electrolyte::updateBoundary()
{   
    i_.correctBoundaryConditions();

    forAll(anodicSpecies_, speciesI)
    {
        word nameSpecie = anodicSpecies_[speciesI];
        x(nameSpecie).correctBoundaryConditions();
        
        const FaradayFvPatchField& molarFraction =
            refCast<FaradayFvPatchField>(x(nameSpecie)
            .boundaryField()[anodePatchIndex_]);
            
        const scalarField& molarFluxOnMolarFraction =
            molarFraction.molarFlux();
            
        fvsPatchScalarField& molarFlux =
            N(nameSpecie).boundaryField()[anodePatchIndex_];
            
        forAll(molarFluxOnMolarFraction, faceI)
        {
            molarFlux[faceI] = molarFluxOnMolarFraction[faceI];
        }
    }   
    forAll(cathodicSpecies_, speciesI)
    {
        word nameSpecie = cathodicSpecies_[speciesI];
        x(nameSpecie).correctBoundaryConditions();
        
        const FaradayFvPatchField& molarFraction =
            refCast<FaradayFvPatchField>(x(nameSpecie)
            .boundaryField()[cathodePatchIndex_]);
            
        const scalarField& molarFluxOnMolarFraction =
            molarFraction.molarFlux();
            
        fvsPatchScalarField& molarFlux =
            N(nameSpecie).boundaryField()[cathodePatchIndex_];
            
        forAll(molarFluxOnMolarFraction, faceI)
        {
            molarFlux[faceI] = molarFluxOnMolarFraction[faceI];
        }      
    }   
    
    IOdictionary thermophysicalProperties
    (
        IOobject
        (
            "thermophysicalProperties",
            mesh().time().constant(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    
    setS();
    dimensionedScalar thickness("thickness", dimLength, etaOhm_.thickness(1));
    q_ = sqr(i_) / solid()->sigma(T_) * thickness;    
    q_.boundaryField()[anodePatchIndex_] = eta_.boundaryField()[anodePatchIndex_];
    q_.boundaryField()[cathodePatchIndex_] = eta_.boundaryField()[cathodePatchIndex_];
    forAll(anodicSpecies_, speciesI)
    {
        word nameSpecie = anodicSpecies_[speciesI];
        const FaradayFvPatchField& xBoundary = refCast<FaradayFvPatchField>(x(nameSpecie).boundaryField()[anodePatchIndex_]);
        const scalarField& xBoundaryMolarFlux = xBoundary.molarFlux();
        const scalarField& A = mesh().magSf().boundaryField()[anodePatchIndex_];        
        specie sp(thermophysicalProperties.lookup(nameSpecie));
        q_.boundaryField()[anodePatchIndex_] += sp.W() *  xBoundaryMolarFlux * (-1) * sy(nameSpecie).boundaryField()[anodePatchIndex_] / A;
    }
    forAll(cathodicSpecies_, speciesI)
    {
        word nameSpecie = cathodicSpecies_[speciesI];
        const FaradayFvPatchField& xBoundary = refCast<FaradayFvPatchField>(x(nameSpecie).boundaryField()[cathodePatchIndex_]);
        const scalarField& xBoundaryMolarFlux = xBoundary.molarFlux();
        const scalarField& A = mesh().magSf().boundaryField()[cathodePatchIndex_];        
        specie sp(thermophysicalProperties.lookup(nameSpecie));
        q_.boundaryField()[cathodePatchIndex_] += sp.W() * xBoundaryMolarFlux * (-1) * sy(nameSpecie).boundaryField()[cathodePatchIndex_] / A;
    }    
    T_.correctBoundaryConditions();
}


void Foam::SOFC::electrolyte::write()
{
    forAll(anodicSpecies_, i)
    {
        x(anodicSpecies_[i]).write();
        N(anodicSpecies_[i]).write();
    }
    forAll(cathodicSpecies_, i)
    {
        x(cathodicSpecies_[i]).write();
        N(cathodicSpecies_[i]).write();
    }
  
    solidComponent::write();
    T_.write();
    
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
