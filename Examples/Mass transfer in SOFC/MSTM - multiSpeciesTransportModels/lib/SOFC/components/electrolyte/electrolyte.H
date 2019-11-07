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

Class
    Foam::SOFC::electrolyte

Description
    electrolyte component for solid oxide fuel cell

SourceFiles
    electrolyte.C

\*---------------------------------------------------------------------------*/

#ifndef electrolyte_H
#define electrolyte_H

#include "fluidComponent.H"
#include "solidComponent.H"

#include "activationOverpotentials.H"
#include "concentrationOverpotentials.H"
#include "ohmicOverpotential.H"

#include "FaradayFvPatchField.H"

#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace SOFC
{
   
/*---------------------------------------------------------------------------*\
                          Class electrolyte Declaration
\*---------------------------------------------------------------------------*/  
  
class electrolyte
:

    public solidComponent
{

protected:

    // Protected data  
  
        //- Numerical tolerance for current density
        static const scalar iToll_;
        
        //- Numerical tolerance for voltage
        static const scalar VToll_;
    
        //- Anodic gas channel component
        channel& anodeChannel_;
        
        //- Cathodic gas channel component
        channel& cathodeChannel_;
    
        //- Anodic reactions
        wordList anodeReactions_;
        
        //- Cathodic reactions
        wordList cathodeReactions_;
        
        //- Species
        wordList speciesName_;
        
        //- Anodic species
        wordList anodicSpecies_;
        
        //- Cathodic species
        wordList cathodicSpecies_;      
        
        //- Index of the anodic boundary
        label anodePatchIndex_;
        
        //- Index of the cathodic boundary
        label cathodePatchIndex_;   
        
        //- Anodic activation overpotential
        activationOverpotentials etaActAnode_;
        
        //- Anodic activation overpotential
        activationOverpotentials etaActCathode_; 
      
        //- Anodic concentratio overpotential
        concentrationOverpotentials etaConcAnode_;
        
        //- Cathodic concentratio overpotential
        concentrationOverpotentials etaConcCathode_;       
        
        //- Ohmic overpotential
        ohmicOverpotential etaOhm_; 
        
        //- Temperature [K]
        volScalarField T_;
        
        //- Thermal flux
        volScalarField q_;
        
        //- Overpotential
        volScalarField eta_;

        //- Nernst voltage
        scalar NernstVoltage_;
        
        //- Open circuit voltage 
        scalar Eocv_;
        
        //- Limiting current density
        scalar iLim_;

        //- Report file
        OFstream report_;    
        
        //- Voltage
        scalar V_;

        //- Molar fractions
        PtrList<volScalarField> x_;
        
        //- Molar fluxes
        PtrList<surfaceScalarField> N_;
        
        //- Molar sources
        PtrList<volScalarField> sy_;

        //- Reversible voltage
        scalar Vrev_;
        
    // Protected Member Functions

        //- Set thermal source
        void setS();

        //- Disallow default bitwise copy construct
        electrolyte(const electrolyte&);

        //- Disallow default bitwise assignment
        void operator=(const electrolyte&);

        
public:
      
      // Runtime type information
      TypeName("electrolyte");


    // Constructors

        //- Construct from components
        electrolyte
        (
            const word name,
            regionPropertiesSOFC& regions,
            const dictionary& dic,
            solidComponent& anode,
            solidComponent& cathode,
            channel& anodeChannel,
            channel& cathodeChannel
        );

        
    //- Destructor

        ~electrolyte()
        {}


    // Member Functions
    
        //- Return a reference to the anodic species names
        const wordList& anodicSpecies() const
        {
          return anodicSpecies_;
        }
      
        //- Return a reference to the cathodic species names
        const wordList& cathodicSpecies() const
        {
          return cathodicSpecies_;
        }
        
        //- Return a reference to the anodic activation overpotential
        const activationOverpotentials& etaActAnode() const
        {
          return etaActAnode_;
        }   
      
        //- Return a reference to the cathodic activation overpotential
        const activationOverpotentials& etaActCathode() const
        {
          return etaActCathode_;
        }     
        
        //- Return a reference to the anodic concentration overpotential
        const concentrationOverpotentials& etaConcAnode() const
        {
          return etaConcAnode_;
        }   
        
        //- Return a reference to the cathodic concentration overpotential
        const concentrationOverpotentials& etaConcCathode() const
        {
          return etaConcCathode_;
        }   
      
        //- Return a reference to the temperature
        volScalarField& T()
        {
            return (T_);
        }
        
        //- Return a reference to the heat flux
        const volScalarField& q() const
        {
          return (q_);
        }
          
        //- Return a reference to molar fraction of a species selected by name
        volScalarField& x(word nome);
        
        //- Return a reference to molar flux of a species selected by name
        surfaceScalarField& N(word nome);
        
        //- Return a reference to molar source of a species selected by name
        const volScalarField& sy(word name) const;
        
        //- Return a reference to molar fraction of a species selected by name
        //- like a scalarField (internalField of 1D domain)
        scalarField& X(word name)
        {
            return (x(name).internalField());
        }
    
        //- Update molar fractions values (also for submodels)
        void updateMolarFraction
        (
            fluidComponent& anode,
            fluidComponent& cathode
        );    

        //- Update the temperature
        void updateTemperature(fluidComponent& anode, fluidComponent& cathode);

        //- Calculate the new current density for given voltage value
        void updateCurrentDensity(scalar V);

        //- calculate the new voltage for given current densit
        void updateVoltage(scalar I);
        
        //- Solve the electrochemical problem in term of current density
        void solve(scalar i)
        {
          updateVoltage(i);
        
          scalar averageCurrentDensity = currentDensity().value();
          if (Pstream::master())
          {
            report_ << "Current density   " << averageCurrentDensity << tab 
                    << "Voltage   " << V_ << endl;
          }
        }        
            
        //- Solve the electrochemical problem in term of voltage
        void solveV(scalar V)
        {
          updateCurrentDensity(V);
          V_ = V;
          report_ << "Voltage   " << V << tab            
                  << "Current density   " << currentDensity().value() << endl; 
        }
    
        //- Set the boundary conditions for molar fractions and temperature
        void setBoundaryConditions
        (
            fluidComponent& anode,
            fluidComponent& cathode
        )
        {
          updateMolarFraction(anode,cathode);
          updateTemperature(anode,cathode);
        }
    
        //- Calculate the boundary molar fluxes
        void updateBoundary();
        
        //- Return the current density (dimensioned type)
        dimensionedScalar currentDensity()
        {
          return( i_.weightedAverage(mesh().V()) );
        } 
        
        // Return the voltage (dimensioned type)
        dimensionedScalar voltage()
        {
          return
          (
              dimensionedScalar("V", dimensionSet(0, -2, 0, 0, 0, 1, 0), V_)
          );
        } 


    // I-O
        
        //- Write fields on file
        void write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace SOFC

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
