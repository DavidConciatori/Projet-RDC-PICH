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

#include "porousSIMPLE.H"
#include "addToRunTimeSelectionTable.H"
#include "adjustPhi.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(porousSIMPLE, 0);
    addToRunTimeSelectionTable(momentumTransport, porousSIMPLE, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porousSIMPLE::porousSIMPLE
(
    hsCombustionThermoSOFC& thermo,
    const fvMesh& mesh
)
:
    momentumTransport(thermo, mesh),
    
    pZones_(mesh)
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::porousSIMPLE::momentumLoop()
{   
    volVectorField& U = U_;
    volScalarField& p = thermo_.p();
    surfaceScalarField& phi = phi_;
    
    const volScalarField& rho = rho_;
    const volScalarField& mu = thermo_.mu();    
   
    scalar maxResidual = 0;
    scalar eqnResidual = 1;
    
    p.storePrevIter();  
    
    tmp<fvVectorMatrix> UEqn
    ( 
        fvm::div(phi, U, "div(phi,U)")
      - fvm::SuSp(fvc::div(phi), U)
      - fvm::laplacian(mu, U, "laplacian(mu,U)")
    );

    UEqn().relax();

    pZones_.addResistance(UEqn());

    solve(UEqn() == -fvc::grad(p)).initialResidual();
    maxResidual = max(eqnResidual, maxResidual);

    volScalarField trAU = 1.0/UEqn().A();
    trAU.rename("rAU");
    surfaceScalarField rhorUAf("(rho*(1|A(U)))", fvc::interpolate(rho*trAU));

    U = trAU*UEqn().H();
    UEqn.clear();
    
    phi = fvc::interpolate(rho*U) & U.mesh().Sf();
    adjustPhi(phi, U, p);
    
    for (int nonOrth=0; nonOrth<=nNonOrthCorr_; nonOrth++)
    {
        tmp<fvScalarMatrix> tpEqn;
        tpEqn = (fvm::laplacian(rho*trAU, p) == fvc::div(phi + surfMassSource_));

        // retain the residual from the first iteration
        if (nonOrth == 0)
        {
            eqnResidual = tpEqn().solve().initialResidual();
            maxResidual = max(eqnResidual, maxResidual);
        }
        else
        {
            tpEqn().solve();
        }

        if (nonOrth == nNonOrthCorr_)
        {
            // Calculate the conservative fluxes
            phi -= tpEqn().flux();

            // Explicitly relax pressure for momentum corrector
            p.relax();

            U -= trAU*fvc::grad(p);
            U.correctBoundaryConditions();

            //bound(p, pMin);
        }
    }
    
    return (maxResidual);
}
  
  
// ************************************************************************* //
