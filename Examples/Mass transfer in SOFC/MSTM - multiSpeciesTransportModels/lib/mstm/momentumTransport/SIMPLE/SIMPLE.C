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

#include "SIMPLE.H"
#include "addToRunTimeSelectionTable.H"
#include "adjustPhi.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(SIMPLE, 0);
    addToRunTimeSelectionTable(momentumTransport, SIMPLE, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SIMPLE::SIMPLE
(
    hsCombustionThermoSOFC& thermo,
    const fvMesh& mesh
)
:
    momentumTransport(thermo, mesh)
{
    Info<< "Creating turbulence model\n" << endl;
    
    // WARNING!!! thermo.rho() is not the updated rho...
    turbulence_ = 
    (
        compressible::RASModel::New
        (
            thermo.rho(),
            U_,
            phi_,
            thermo
        )
    );
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::SIMPLE::momentumLoop()
{   
    volVectorField& U = U_;
    volScalarField& p = thermo_.p();
    surfaceScalarField& phi = phi_;
    
    const volScalarField& rho = rho_;
      
    scalar maxResidual = 0;
    scalar eqnResidual = 1;
      
    p.storePrevIter();

    tmp<fvVectorMatrix> UEqn
    ( 
        fvm::div(phi, U, "div(phi,U)")
	+ turbulence_->divDevRhoReff(U)
//      - fvm::SuSp(fvc::div(phi), U)
//      - fvm::laplacian(mu, U, "laplacian(mu,U)")
//      - fvc::div(mu*dev2(fvc::grad(U)().T()))   
    );
    
//        - fvm::laplacian(nuEff(), U)
//       - fvc::div(nuEff()*dev(fvc::grad(U)().T()))   
  
    UEqn().relax();

    eqnResidual = solve(UEqn() == -fvc::grad(p)).initialResidual();
    maxResidual = max(eqnResidual, maxResidual);
  

    volScalarField rUA = 1.0/UEqn().A();
    surfaceScalarField rhorUAf("(rho*(1|A(U)))", fvc::interpolate(rho*rUA));
  
    U = rUA*UEqn().H();
    UEqn.clear();
  
    phi = fvc::interpolate(rho)*(fvc::interpolate(U) & U.mesh().Sf());
    adjustPhi(phi, U, p);
    
    // Non-orthogonal pressure corrector loop
    for (int nonOrth=0; nonOrth<=nNonOrthCorr_; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rho*rUA, p) == fvc::div(phi + surfMassSource_)
        );
    
        if (nonOrth == 0)
        {
            eqnResidual = pEqn.solve().initialResidual();
            maxResidual = max(eqnResidual, maxResidual);
        }
        else
        {
            pEqn.solve();
        }

        if (nonOrth == nNonOrthCorr_)
        {
          // Calculate the conservative fluxes
          phi -= pEqn.flux();
        } 
    }
 
    // Explicitly relax pressure for momentum corrector
    p.relax();
 
    U -= rUA*fvc::grad(p);
    U.correctBoundaryConditions();
//     Info << fvc::interpolate(rho);
//     Info << fvc::interpolate(thermo.rho());
    
 //   turbulence_().correct();
   
    return (maxResidual);
}


// ************************************************************************* //
