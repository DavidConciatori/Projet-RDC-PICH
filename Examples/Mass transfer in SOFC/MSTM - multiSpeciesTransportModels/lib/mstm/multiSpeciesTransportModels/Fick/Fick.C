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

#include "Fick.H"
#include "addToRunTimeSelectionTable.H"
#include "fixedFluxFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Fick, 0);
    addToRunTimeSelectionTable
    (
        multiSpeciesTransportModel, Fick, dictionary
    );
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::Fick::updateCoefficients()
{
    DijModel_().update();

    forAll(mix().species(), i)
    {
        if(mix().species()[i] != inertSpecie_)
        {
            volScalarField tmpGamma = 0 / Dij(0,0);
            forAll(mix().species(), j)
            {
                if (j != i)
                {
                    tmpGamma += mix().x()[j] / Dij(i,j);
                }
            }
            Gamma_.set(i, mt().rho() * (1-mix().x()[i]) / tmpGamma );
            Gamma_[i].rename("Gamma_" + mix().species()[i]);

            // Instantiated only for compatibility while using
            // fixedFlux boundary condition
            phi_N_.set(i, mt().phi());
            phi_N_[i].rename("phi_N_" + mix().species()[i]);

            // Instantiated only for compatibility while using
            // fixedFlux boundary condition    
            phi_P_.set(i, 0*mt().phi());
            phi_P_[i].rename("phi_P_" + mix().species()[i]);
        }
    }
} 


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Fick::Fick
(
    const fvMesh& mesh
)
:
    multiSpeciesTransportModel(mesh)
{    
    bases_ = "mass";  
  
    // Construct the therm pakage 
    word thermoTypeName = 
    "hsPsiMixtureThermoSOFC<multiComponentMixtureMassBases<gasThermoPhysicsSOFC>>";
    thermo_ = hsCombustionThermoSOFC::New(mesh, thermoTypeName);
    
    // Construct the momentum (and continuity) solver    
    mt_ = momentumTransport::New(thermo_(), mesh);
    
    // Construct the diffusivity model here because thermo pakage
    // isn't available before
    DijModel_.set(new diffusivityModel(thermo().p(), thermo().T()));
    
    // Construct the partial pressures fields for compatibility with the models
    // that use partial pressures while coupling regions
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

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
 
void Foam::Fick::correct()
{
    multiComponentMixtureMassBases<gasThermoPhysicsSOFC>& thermoMix = 
    reinterpret_cast
    <
        multiComponentMixtureMassBases<gasThermoPhysicsSOFC>&
    >(mix());
    
    scalar nit = 0;
    bool convergenceFlag = 1;
    
    bool mustSolveMomentum = 1; 
    bool mustSolveMasstransport = 1; 
                 
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
        thermoMix.correctMolarFluxes();

        mt().rho() = mix().W() * thermo().p() / (RR * thermo().T());

        forAll(mix().species(), i)
        {
            p_[i] = mix().x(i) * thermo().p();
        }

        updateCoefficients();

        if (maxResidual < convergenceCriterion_)
        {
          convergenceFlag = 0;
        }

        nit++;
    }        
}
 
 
Foam::scalar Foam::Fick::massTransportLoop()
{       
    scalar maxResidual = 0;
    scalar eqnResidual = 1;

    updateCoefficients(); 
    label inertIndex = -1;
    volScalarField yt = 0.0*mix().y()[0];
    surfaceScalarField nt = mt().phi();
                
//    multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

    forAll(mix().species(), i)
    {        
        volScalarField& yi = mix().y(i);
        surfaceScalarField& ni = mix().n(i);
      
        if (mix().species()[i] != inertSpecie_)
        {
//         tmp<fv::convectionScheme<scalar> > mvConvection
//         (
//             fv::convectionScheme<scalar>::New
//             (
//                 mt_().phi().mesh(),
//                 fields,
//                 mt_().phi(),
//                 mt_().phi().mesh().divScheme("div(phi,yAlpha)")
//             )
//         );
        
            yi.storePrevIter();

            tmp<fvScalarMatrix> yEqn
            (
                fvm::div(mt_().phi(), yi, "div(phi,yAlpha)")
              - fvm::laplacian(Gamma_[i],yi, "laplacian(D,yAlpha)")
            );

            yEqn() -= Sy_[i];          
                              
            eqnResidual = solve(yEqn() , yi.mesh().solver("yAlpha")).initialResidual();
            maxResidual = max(eqnResidual, maxResidual);

            yi.max(0.0);
//          yi.min(1.0);
          
            ni = yEqn().flux();
            nt -= ni;

            yi.relax(yi.mesh().relaxationFactor("yAlpha"));
          
            yt += yi;  
        }
        else
        {
            inertIndex = i;
        }
    }
        
    // Calculate inert species
    volScalarField& yInert = mix().y()[inertIndex];
    yInert = 1 - yt;
    yInert.max(0.0);
    mix().n()[inertIndex] = nt;
          
    return maxResidual;
}
   

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
