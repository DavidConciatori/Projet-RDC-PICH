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
    Foam::multiComponentMixtureMassBases

Description
    Multi component mixture mixture that considers
    mass fractions as a main variable.

SourceFiles
    multiComponentMixtureMassBases.C

\*---------------------------------------------------------------------------*/

#ifndef multiComponentMixtureMassBases_H
#define multiComponentMixtureMassBases_H

#include "basicMultiComponentMixtureSOFC.H"
#include "HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class multiComponentMixtureMassBases Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class multiComponentMixtureMassBases
:
    public basicMultiComponentMixtureSOFC
{
    // Private data

        //- Species data
        PtrList<ThermoType> speciesData_;

        //- Temporary storage for the cell/face mixture thermo data
        mutable ThermoType mixture_;


    // Private member functions

        //- Construct the species data from the given dictionary and return the
        //  data for the first specie to initialise the mixture thermo data
        const ThermoType& constructSpeciesData(const dictionary& thermoDict);

        //- Construct the mixture composition
        void constructMixtureComposition(const fvMesh&);	

        //- Construct as copy (not implemented)
        multiComponentMixtureMassBases(const multiComponentMixtureMassBases<ThermoType>&);


public:

    //- The type of thermodynamics this mixture is instantiated for
    typedef ThermoType thermoType;


    // Constructors

        //- Construct from dictionary, specie names, thermo database and mesh
        multiComponentMixtureMassBases
        (
            const dictionary&,
            const wordList& specieNames,
            const HashPtrTable<ThermoType>& specieThermoData,
            const fvMesh&
        );

        //- Construct from dictionary and mesh
        multiComponentMixtureMassBases(const dictionary&, const fvMesh&);


    //- Destructor
    virtual ~multiComponentMixtureMassBases()
    {}


    // Member functions
    
        //- Update molar fractions and the molar weight of the mixture 
        void correct();    

        //- Update molar fluxes
        void correctMolarFluxes();    

        const ThermoType& cellMixture(const label celli) const;

        const ThermoType& patchFaceMixture
        (
            const label patchi,
            const label facei
        ) const;

        //- Return the raw specie thermodynamic data
        const PtrList<ThermoType>& speciesData() const
        {
            return speciesData_;
        }

        //- Read dictionary
        void read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "multiComponentMixtureMassBases.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
