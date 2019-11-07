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

#include "dustyGas.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dustyGas, 0);
    addToRunTimeSelectionTable
    (
        multiSpeciesTransportModel, dustyGas, dictionary
    );
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::dustyGas::updateCoefficients()
{   
    DijModel_().update();

    forAll(mix().species(), i)
    {
        volScalarField tmpGamma = 1 / DK(i);
        forAll(mix().species(), j)
        {
            if(j != i)
            {
                tmpGamma += mix().x(j) / Dij(i,j);
            }
        }
        Gamma_.set(i, 1 / (tmpGamma*RR*thermo().T()));
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
            fvc::interpolate(Gamma_[i] * RR*thermo().T() / thermo().p())
                * tmpPhi
        );
        phi_N_[i].rename("phi_N_" + mix().species()[i]);
      
        phi_P_.set(i, fvc::interpolate(Gamma_[i] /  DK(i)) * phiDarcy_);  
        phi_P_[i].rename("phi_P_" + mix().species()[i]);
    }
  } 


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dustyGas::dustyGas
(
    const fvMesh& mesh
)
:
    multiSpeciesTransportModel(mesh),
    
    B0_(lookup("B0")),
    
    DarcySwitch_(lookupOrDefault<scalar>("DarcySwitch", 1)),
 
    phiDarcy_
    (
      IOobject
      (
        "phiDarcy",
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
      ),
      mesh,
      dimensionedScalar("phiDarcy", dimVolume/dimTime, 0)
    )
{        
    bases_ = "molar";
  
    // Construct the thermo pakage 
    word thermoTypeName = 
    "hsPsiMixtureThermoSOFC<multiComponentMixtureMolarBases<gasThermoPhysicsSOFC>>";
    thermo_ = hsCombustionThermoSOFC::New(mesh, thermoTypeName);
    
    // Construct the momentum (and continuity) solver    
    mt_ = momentumTransport::New(thermo_(), mesh, "Darcy");
        
    // Construct the diffusivity model here because thermo pakage
    // isn't available before
    DijModel_.set(new diffusivityModel(thermo().p(), thermo().T()));
    
    // Construct the Knudsen diffusivity model
    DKModel_.set(new KnudsenDiffusivityModel(thermo().T()));
    
    // Construct the partial pressures fields
    constructP(mesh);
    
    // Initialize the species equations coefficients
    updateCoefficients();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::dustyGas::correct()
{   
    multiComponentMixtureMolarBases<gasThermoPhysicsSOFC>& thermoMix = 
    reinterpret_cast
    <
        multiComponentMixtureMolarBases<gasThermoPhysicsSOFC>&
    >(mix());  
  
    scalar nit = 0;
    bool convergenceFlag = 1;
    
    volScalarField& p = thermo().p();
              
    while( (convergenceFlag) && (nit < nitMax_)) 
    {   
      Info << "\nMass transport n° iteration = " << nit+1 << endl;
      
      scalar maxResidual = massTransportLoop();
      
      thermoMix.correct();      
      thermoMix.correctMassFluxes();
      
      mt().rho() = mix().W() * p / (RR * thermo().T());

      // If grad(p) is used to calculate the Darcy term, boundary conditions
      // for U will be not considered. Alternativelly directly U can be used to
      // calculate the Darcy term. In this case it is necessary to solve the 
      // momentum equation (Darcy) at each iterations ann non olny at the end.
      // Here grad(p) is used.
      if(DarcySwitch_)
      {
          phiDarcy_ = - B0_ / fvc::interpolate(thermo().mu())
              * (fvc::interpolate(fvc::grad(p)) & p.mesh().Sf());
      }  
          
      updateCoefficients();

      if (maxResidual < convergenceCriterion_)
      {
	convergenceFlag = 0;
      }

      nit++;
    }
    
    mt_->momentumSolve(nitMax_);
} 


Foam::scalar Foam::dustyGas::massTransportLoop()
{        
    scalar maxResidual = 0;
    scalar eqnResidual = 1;
    
    updateCoefficients(); 
    thermo().p() *= 0.0;
    
//    multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;
    
    forAll(mix().species(), i)
    {
        volScalarField& pi = p_[i];
        surfaceScalarField& Ni = mix().N(i);     
    
// 	  tmp<fv::convectionScheme<scalar> > mvConvection
// 	  (
// 	      fv::convectionScheme<scalar>::New
// 	      (
// 		  Ni.mesh(),
// 		  fields,
// 		  Ni,
// 		  Ni.mesh().divScheme("div(phi,pAlpha)")
// 	      )
// 	  );	

        pi.storePrevIter();
        Ni.storePrevIter();

        tmp<fvScalarMatrix> pEqn
        (
            fvm::div(phi_N_[i], pi, "div(phi,pAlpha)")
          + fvm::div(phi_P_[i], pi, "div(phi,pAlpha)")
          - fvm::laplacian(Gamma_[i], pi, "laplacian(D,pAlpha)")
        );    

        pEqn() -= Sx_[i];

        eqnResidual = solve(pEqn() , pi.mesh().solver("pAlpha")).initialResidual();
        maxResidual = max(eqnResidual, maxResidual);

        pi.max(100); 
//	pi.min(101325);

        Ni = pEqn().flux();

        pi.relax(pi.mesh().relaxationFactor("pAlpha"));
        Ni.relax(Ni.mesh().relaxationFactor("NAlpha"));

        thermo().p() += pi;
    }
    
    forAll(mix().species(), i)
    {
        volScalarField& pi = p_[i];
        volScalarField& xi = mix().x(i);
        xi = pi / thermo().p();	
    }
    
    return maxResidual;    
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
