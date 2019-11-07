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

Application
    stefanTube

Description
    Solver that simulate a Stefan tube
    
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "multiSpeciesTransportModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Useful patch index and area
    label indexLiquid = mesh.boundaryMesh().findPatchID("liquid");
    label indexExit = mesh.boundaryMesh().findPatchID("exit");
    scalar A = mesh.magSf().boundaryField()[indexExit][0];
    
    
    // The multi species transport model
    autoPtr<multiSpeciesTransportModel> mstm_ = multiSpeciesTransportModel::New(mesh);
    
    
    // Solving the problem...
    mstm_->correct();
    
    
    // Write the results
    runTime++;
    mstm_->write();
    
    
    // Print on screen mass fluxes balance
    Info << "\n\n";
    forAll(mstm_().mix().species(), specieI)
    {
      surfaceScalarField phi = mstm_().mix().n(specieI);
      
      Info << "\n\n" << mstm_().mix().species()[specieI] << " mass flow [kg/s]: " << token::TAB 
           << " inlet = " << sum(phi.boundaryField()[indexLiquid]) << token::TAB 
           << " outlet = " << sum(phi.boundaryField()[indexExit]) << token::TAB << token::TAB << token::TAB
            << " balance = " << sum(phi.boundaryField()[indexLiquid] + phi.boundaryField()[indexExit]);
    }
      
    surfaceScalarField phi = mstm_().mt().phi();
 
    Info << "\n\nGlobal mass flow [kg/s]: " << token::TAB 
         << " inlet = " << sum(phi.boundaryField()[indexLiquid]) << token::TAB  
         << " outlet = " << sum(phi.boundaryField()[indexExit]) << token::TAB << token::TAB << token::TAB
         << " balance = " << sum(phi.boundaryField()[indexLiquid] + phi.boundaryField()[indexExit]);
    

    // Print on screen molar fluxes balance
    Info << "\n\n";
    forAll(mstm_().mix().species(), specieI)
    {    
      Info << "\n\n" << mstm_().mix().species()[specieI] << " mass flow = " << mstm_().mix().n(specieI).boundaryField()[indexExit][0] << " kg/s " << token::TAB 
           << "(" << mstm_().mix().n(specieI).boundaryField()[indexExit][0] / A << " kg/s/m^2)";

      Info << "\n\n" << mstm_().mix().species()[specieI] << " molar flow = " << mstm_().mix().N(specieI).boundaryField()[indexExit][0] << " kmol/s " << token::TAB 
           << "(" << mstm_().mix().N(specieI).boundaryField()[indexExit][0] / A << " kmol/s/m^2)";
    }
   
    // Print on screen velocity
    Info << "\n\n";
    Info << "\n\nU at inlet = " << mstm_().mt().U().boundaryField()[indexLiquid][0];
    Info << "\n\nU at outlet = " << mstm_().mt().U().boundaryField()[indexExit][0];


    Info << "\n\n";
    Info<< "\nExecutionTime = "
    << runTime.elapsedCpuTime()
    << " s\n\n" << endl;

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
