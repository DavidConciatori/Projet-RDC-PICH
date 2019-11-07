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
    Foam::multiSpeciesTransportModel

Description
    Abstract base class for multispecies mass transport model models.

SourceFiles
    multiSpeciesTransportModelI.H
    multiSpeciesTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef multiSpeciesTransportModel_H
#define multiSpeciesTransportModel_H

#include "fvc.H"
#include "fvm.H"

#include "hsCombustionThermoSOFC.H"
#include "basicMultiComponentMixtureSOFC.H"
#include "thermoPhysicsTypesSOFC.H"

#include "diffusivityModel.H"
#include "KnudsenDiffusivityModel.H"
#include "momentumTransport.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  
/*---------------------------------------------------------------------------*\
                    Class multiSpeciesTransportModel Declaration
\*---------------------------------------------------------------------------*/

class multiSpeciesTransportModel
:

    public IOdictionary
{
    
protected:  
  
    // Protected data
      
        //- Static dimensioned universal gas constant [J/(kmol K)]
        static const dimensionedScalar RR;   
    
        //- Thermo package
        autoPtr<hsCombustionThermoSOFC> thermo_;  
    
        //- Momentum transport package
        autoPtr<momentumTransport>  mt_;  
    
        //- Partial pressures [Pa]
        PtrList<volScalarField> p_;

        //- Species mass sources [kg/s/m]
        PtrList<volScalarField> Sy_;

        //- Species molar sources [kmol/s/m]
        PtrList<volScalarField> Sx_;
    
        //- Diffusivity models 
        autoPtr<diffusivityModel> DijModel_;

        //- Diffusivity models 
        autoPtr<KnudsenDiffusivityModel> DKModel_;
    
        //- Coefficients of the partial pressures gradients
        PtrList<volScalarField> Gamma_;
    
        //- Coefficients of the coupling fluxes
        PtrList<surfaceScalarField> phi_N_;
    
        //- Coefficients of the global pressures gradients
        PtrList<surfaceScalarField> phi_P_;

        //- Mass/molar bases information useful to generalize
        //  boundary conditions
        word bases_;

        //- Number of species in the mixture
        label ns_;

        //- Convergence tolerance
        scalar convergenceCriterion_;

        //- Max iterations number
        label nitMax_;
    

    // Protected member functions

        //- Construct the partial pressures
        void constructP(const fvMesh&);
      
 
private:
  
    // Private Member Functions
    
        //- Disallow copy construct
        multiSpeciesTransportModel(const multiSpeciesTransportModel&);

        //- Disallow default bitwise assignment
        void operator=(const multiSpeciesTransportModel&);


public:
 
      //- Runtime type information
      TypeName("multiSpeciesTransportModel");

      
      // Declare run-time constructor selection table
     
          declareRunTimeSelectionTable
          (
              autoPtr,
              multiSpeciesTransportModel,
              dictionary,
              (
                  const fvMesh& mesh
              ),
              (mesh)
          );
  

    // Constructors	  

        //- Construct from components
        multiSpeciesTransportModel
        (
            const fvMesh& mesh
        );  


    // Selectors	  
  
        //- Return a reference to the selected multispecies transport model

        static autoPtr<multiSpeciesTransportModel> New
        (
            const fvMesh& mesh
        );


    //- Destructor

        virtual ~multiSpeciesTransportModel()
        {}


    // Member Functions

        // Access
        
            //- Return the bases variable of the model
            inline const word& bases() const;        
           
            //- Return the const binary diffusion coefficient field
            //  for a couple of species given by indexs
            inline const volScalarField& Dij(label i, label j) const;

            //- Return the const Knudsen diffusion coefficient field
            //  given by index
            inline const volScalarField& DK(label i) const;	

            //- Return the thermo package
            inline hsCombustionThermoSOFC& thermo();    

            //- Return the const thermo package
            inline const hsCombustionThermoSOFC& thermo() const;	

            //- Return the momentumTransport package
            inline momentumTransport& mt();
      
            //- Return the const momentumTransport package
            inline const momentumTransport& mt() const;	

            //- Return a pointer to the specieMixture object
            inline basicMultiComponentMixtureSOFC& mix();
    
            //- Return const-access to the specieMixture object
            inline const basicMultiComponentMixtureSOFC& mix() const;

            //- Return const-access to the species mass sources
            inline PtrList<volScalarField>& Sy();

            //- Return const-access to the species mass sources
            inline const PtrList<volScalarField>& Sy() const;

            //- Return const-access to the species molar sources
            inline PtrList<volScalarField>& Sx();

            //- Return const-access to the species molar sources
            inline const PtrList<volScalarField>& Sx() const;
  
        //- Return const-access to the inert specie
        virtual const word inertSpecie() const;
      
        //- Return the thermodynamic species data
        virtual const PtrList<gasThermoPhysicsSOFC>& speciesThermo() const = 0;

        //- Solve the mass transport
        virtual void correct() = 0;


    // I-O

        //- Write the mixture fields      
        void write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "multiSpeciesTransportModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
