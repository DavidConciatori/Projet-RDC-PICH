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

#include "MaxwellStefan.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(MaxwellStefan, 0);
    addToRunTimeSelectionTable
    (
        multiSpeciesTransportModel, MaxwellStefan, dictionary
    );
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::MaxwellStefan::updateCoefficients()
{     
    DijModel_().update();

    forAll(mix().species(), i)
    {
        if(mix().species()[i] != inertSpecie_)
        {
            volScalarField tmpGamma = 0.0 / Dij(0,0);
            forAll(mix().species(), j)
            {
                if(j != i)
                {
                    tmpGamma += mix().x(j) / Dij(i,j);
                }
            }
            Gamma_.set(i, 1 / (tmpGamma * RR * thermo().T()));
            Gamma_[i].rename("Gamma_" + mix().species()[i]);

            surfaceScalarField tmpPhi =
                0.0 * mix().N(0) / fvc::interpolate(Dij(0,0));
            forAll(mix().species(), j)
            {
                if(j != i)
                {
                     tmpPhi += mix().N(j) / fvc::interpolate(Dij(i,j));
                }
            }
            phi_N_.set
            (
                i,
                fvc::interpolate(Gamma_[i]*RR*thermo().T()/thermo().p())
                    *tmpPhi
            );
            phi_N_[i].rename("phi_N_" + mix().species()[i]);
  
            dimensionedScalar Wi
            (
                "Wi",
                dimMass/dimMoles,
                speciesThermo()[i].W()
            );
            phi_P_.set
            (
                i,
                fvc::interpolate(Gamma_[i] * Wi / mix().W() / thermo().p())
                    * fvc::snGrad(thermo().p()) * Gamma_[i].mesh().magSf()
            );  
            phi_P_[i].rename("phi_P_" + mix().species()[i]);
        }
    }
}
  
  
void Foam::MaxwellStefan::correctInertFlux()
{  
    label inertIndex = -1;
    mix().n(inertSpecie_) = mt().phi();
    forAll(mix().species(), i)
    {    
        if (mix().species()[i] != inertSpecie_)
        { 
            mix().n(inertSpecie_) -= mix().n(i);
        }
        else
        {
          inertIndex = i;
        }
    }
    dimensionedScalar WInert
    (
        "Wnn",
        dimMass/dimMoles,
        speciesThermo()[inertIndex].W()
    );
    mix().N(inertIndex) = mix().n(inertIndex) / WInert; 
}
  

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MaxwellStefan::MaxwellStefan
(
    const fvMesh& mesh
)
:

    multiSpeciesTransportModel(mesh)  
{   
    bases_ = "molar";
    
    // Construct the thermo pakage 
    word thermoTypeName = 
    "hsPsiMixtureThermoSOFC<multiComponentMixtureMolarBases<gasThermoPhysicsSOFC>>";
    thermo_ = hsCombustionThermoSOFC::New(mesh, thermoTypeName);
    
    // Construct the momentum (and continuity) solver    
    mt_ = momentumTransport::New(thermo_(), mesh);
        
    // Construct the diffusivity model here because thermo pakage
    // isn't available before
    DijModel_.set(new diffusivityModel(thermo().p(), thermo().T()));
    
    // Construct the partial pressures fields
    constructP(mesh);
    
    // The last specie in the list is the default inertspecie
    inertSpecie_ = lookupOrDefault<word>
    (
        "inertSpecie",
        mix().species()[mix().species().size()-1]
    );
    
    // Initialize the species equations coefficients
    updateCoefficients();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::MaxwellStefan::correct()
{  
    multiComponentMixtureMolarBases<gasThermoPhysicsSOFC>& thermoMix = 
    reinterpret_cast
    <
        multiComponentMixtureMolarBases<gasThermoPhysicsSOFC>&
    >(mix());   
    
    scalar nit = 0;
    bool convergenceFlag = 1;
    
    bool mustSolveMomentum = 1; 
    bool mustSolveMasstransport = 1; 

    volScalarField& p = thermo().p();
    
    while( (convergenceFlag) && (nit < nitMax_)) 
    {   
        scalar maxResidual = 0;
      
        Info << "\nMass transport n° iteration = " << nit+1 << endl;

        if(mustSolveMomentum)
        {
            scalar eqnResidual = 0;
  
            eqnResidual = mt().momentumLoop();
            maxResidual = max(eqnResidual, maxResidual);

            if(eqnResidual < convergenceCriterion_)
            {
                mustSolveMomentum = 0;
                mustSolveMasstransport = 1;
            }        
        }
        
        if(mustSolveMasstransport)
        {
            scalar eqnResidual = 0;

            eqnResidual = massTransportLoop();
            maxResidual = max(eqnResidual, maxResidual);

            if(eqnResidual < convergenceCriterion_)
            {
                mustSolveMasstransport = 0;
                mustSolveMomentum = 1;
                eqnResidual = mt().momentumLoop();
                maxResidual = max(eqnResidual, maxResidual);
            }
        }
           
        thermoMix.correct();      
        thermoMix.correctMassFluxes();
        
//        mt().rho().storePrevIter();
        mt().rho() = mix().W() * p / (RR * thermo().T());
//        mt().rho().relax();

        correctInertFlux();

        updateCoefficients();

        if (maxResidual < convergenceCriterion_)
        {
          convergenceFlag = 0;
        }
        
        nit++;
    }  
} 
    

Foam::scalar Foam::MaxwellStefan::massTransportLoop()
{        
    scalar maxResidual = 0;
    scalar eqnResidual = 1;
 
    updateCoefficients(); 
    label inertIndex = -1;
    volScalarField pt = 0.0*thermo().p();     
      
//     multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;
   
    forAll(mix().species(), i)
    {
        volScalarField& pi = p_[i];
        volScalarField& xi = mix().x(i);
        surfaceScalarField& Ni = mix().N(i);    
    
        if (mix().species()[i] != inertSpecie_)
        {    
/*            tmp<fv::convectionScheme<scalar> > mvConvection
            (
                fv::convectionScheme<scalar>::New
                (
                    Ni.mesh(),
                    fields,
                    phi_P_[i]+phi_N_[i], // Ni,
                    Ni.mesh().divScheme("div(phi,pAlpha)")
                )
            );  */      

            pi.storePrevIter();
            Ni.storePrevIter();

            tmp<fvScalarMatrix> pEqn
            (
                fvm::div(phi_N_[i], pi, "div(phi,pAlpha)")
              + fvm::div(phi_P_[i], pi, "div(phi,pAlpha)")
              - fvm::laplacian(Gamma_[i], pi, "laplacian(D,pAlpha)")      
//              - fvc::SuSp(fvc::div(Ni), xi)            
//              + mvConvection->fvmDiv(phi_N_[i], pi)
//              + mvConvection->fvmDiv(phi_P_[i], pi)     
            );
  
            pEqn() -= Sx_[i];  
              
            eqnResidual = solve(pEqn() , pi.mesh().solver("pAlpha")).initialResidual();
            maxResidual = max(eqnResidual, maxResidual);
            
            pi.max(100); 
//          pi.min(101325);
  
            Ni = pEqn().flux();
            
            pEqn.clear();
  
            pi.relax(pi.mesh().relaxationFactor("pAlpha"));

            Ni.relax(Ni.mesh().relaxationFactor("NAlpha"));
  
            xi = pi / thermo().p();
  
            pt += pi;
        }
        else
        {
            inertIndex = i;
        }
    }
    
    // Calculate inert species
    volScalarField& pInert = p_[inertIndex];
    pInert = thermo().p() - pt;
    pInert.max(100);
    mix().x(inertIndex) = p_[inertIndex] / thermo().p();
  
    return maxResidual;    
} 
    
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
