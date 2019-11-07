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
    sofcFoam

Description
    Solver that can simulate solid oxide fuel cells (SOFC)
    
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "OFstream.H"
#include "regionPropertiesSOFC.H"

#include "channel.H"
#include "electrode.H"
#include "electrolyte.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
  
#   include "setRootCase.H"  
#   include "createTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    regionPropertiesSOFC regions(runTime);
    
#   include "readParameters.H"    

    OFstream polarization("polarization.txt");

    runTime.setTime(startValue,0);

    SOFC::channel anodeChannel("anodeChannel", regions);
    SOFC::channel cathodeChannel("cathodeChannel", regions);
    SOFC::electrode anode("anode", regions,reactiveAnode);
    SOFC::electrode cathode("cathode", regions);
    SOFC::electrolyte electrolyte
        ("electrolyte", regions, electrolyteDictionary,
         anode, cathode, anodeChannel, cathodeChannel);

    anodeChannel.activated_ = anodeChannelSwitch;
    cathodeChannel.activated_ = cathodeChannelSwitch;
    anode.activated_ = anodeSwitch;
    cathode.activated_ = cathodeSwitch;
    electrolyte.activated_ = 1;

    runTime.setTime(0,0);
    
    
    forAll(values, index)
    {
        scalar oldValue = -10000;
        scalar newValue = 0;
        
    while ( mag(newValue - oldValue) > topLevelTollerance[index] )
    { 

        if (electricalSwitch)
        {        
            Info << "\n\nSolving electrolyte\n\n";
            electrolyte.setBoundaryConditions(anode,cathode);
            if(solverSwitch)
            {
                oldValue = electrolyte.voltage().value();
                electrolyte.solve(values[index]);
                newValue = electrolyte.voltage().value();
            }
            else
            {
                oldValue = electrolyte.currentDensity().value();
                electrolyte.solveV(values[index]);
                newValue = electrolyte.currentDensity().value();
            }
            electrolyte.updateBoundary();
            electrolyte.write();
        }
      
      
        if (massTransportSwitch)
        {
          anode.solveMassTransport
              (anodeChannel, electrolyte, "electrolyte_to_anode");
              
          cathode.solveMassTransport
              (cathodeChannel, electrolyte, "electrolyte_to_cathode");
          
          anodeChannel.solveMassTransport(anode);
          cathodeChannel.solveMassTransport(cathode);
          
          anode.solveMassTransport
              (anodeChannel, electrolyte, "electrolyte_to_anode");
              
          cathode.solveMassTransport
              (cathodeChannel, electrolyte, "electrolyte_to_cathode");
        }   


        electrolyte.write();
        anode.write();
        anodeChannel.write();
        cathode.write();
        cathodeChannel.write();

        
        if (printFlowSwitch)
        {
          anode.printMassFlux();
          anodeChannel.printMassFlux();
          cathode.printMassFlux();
          cathodeChannel.printMassFlux();
        }    
    }    
    
    
    scalar averageCurrentDensity = electrolyte.currentDensity().value();
    scalar voltage = electrolyte.voltage().value();
    scalar time = runTime.elapsedCpuTime();

    if (Pstream::master())
    {
       polarization << "Current density   " << averageCurrentDensity << tab
                    << "Cell voltage   " << voltage << tab 
                    << "Time   " << time << endl;         
    }
    
    
    runTime.setTime(values[index],0);
    anodeChannel.write();
    cathodeChannel.write();
    anode.write();
    cathode.write();
    electrolyte.write();
    runTime.setTime(0,0);
  }

  Info<< "\nExecutionTime = "
  << runTime.elapsedCpuTime()
  << " s\n\n" << endl;

  Info<< "End\n" << endl;

  return(0);
}


// ************************************************************************* //
