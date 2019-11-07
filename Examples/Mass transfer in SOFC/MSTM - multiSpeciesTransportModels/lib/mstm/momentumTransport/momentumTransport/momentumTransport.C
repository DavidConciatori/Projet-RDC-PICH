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

#include "momentumTransport.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(momentumTransport, 0);
    defineRunTimeSelectionTable(momentumTransport, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::momentumTransport::momentumTransport
(
    hsCombustionThermoSOFC& thermo,
    const fvMesh& mesh
)
:
    thermo_(thermo),
    
    U_
    (
      IOobject
      (
        "U",
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
      ),
      mesh
    ),
     
    phi_
    (
      IOobject
      (
        "phi",
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
      ),
      mesh,
      dimensionedScalar("phi", dimMass/dimTime, 0)
    ),

    rho_
    (
      IOobject
      (
        "rho",
        mesh.time().timeName(),
        mesh
      ),
      thermo.rho()
    ),
  
    surfMassSource_
    (
      IOobject
      (
        "surfMassSource",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
      ),
      mesh,
      dimensionedScalar("surfMassSource", dimMass/dimTime, 0)
    )
{
  dictionary convergenceCriterionDictionary_ = mesh.solutionDict().subDict("SIMPLE");
  convergenceCriterion_ = convergenceCriterionDictionary_.lookupOrDefault<scalar>("convergence", 0);
    
  dictionary simpleDictionary_ = mesh.solutionDict().subDict("SIMPLE");
  nNonOrthCorr_ = simpleDictionary_.lookupOrDefault<label>("nNonOrthogonalCorrectors", 0);   
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::momentumTransport::momentumSolve(label nitMax)
{       
    scalar nit = 0;
    bool convergenceFlag = 1;
   
    while( (convergenceFlag) && (nit < nitMax)) 
    {
        Info<< "\nNavier Stokes n° iteration = " << nit+1 << endl;
    
        scalar maxResidual = momentumLoop();
    
        if (maxResidual < convergenceCriterion_)
        {
            convergenceFlag = 0;
        }    
        nit++;        
    }
}


void Foam::momentumTransport::printMassFlux(List<word> boundaryName)
{
    Info << "\n\n";
    scalar sumPhi = 0;
    scalar sumBoundaryPhi = 0;
    forAll(boundaryName, boundaryI)
    {
        label index = phi_.mesh().boundaryMesh().findPatchID(boundaryName[boundaryI]);
        scalar massFlow = gSum(phi_.boundaryField()[index]);
        sumBoundaryPhi += gSum(surfMassSource_.boundaryField()[index]);
        Info << "Global mass flow at " << boundaryName[boundaryI] << ": " << massFlow << "\n";
        sumPhi += massFlow;
    }
    Info << "Global mass source: " << sumBoundaryPhi << "\n";
    Info << "Global mass flow balance: " << sumPhi + sumBoundaryPhi;
}


void Foam::momentumTransport::write()
{
    U_.write();
    phi_.write();
    rho_.write();
}


// ************************************************************************* //
