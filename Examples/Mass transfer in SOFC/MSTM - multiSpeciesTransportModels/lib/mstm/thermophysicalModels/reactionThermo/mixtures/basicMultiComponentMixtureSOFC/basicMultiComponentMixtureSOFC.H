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
    Foam::basicMultiComponentMixtureSOFC

Description
    Base class for multi component mixture mixture based on both
    molar and mass fractions. Fluxes (mass and molar) are also defined.

SourceFiles
    basicMultiComponentMixtureSOFC.C

\*---------------------------------------------------------------------------*/

#ifndef basicMultiComponentMixtureSOFC_H
#define basicMultiComponentMixtureSOFC_H

#include "volFields.H"
#include "surfaceFields.H"
#include "PtrList.H"
#include "speciesTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class basicMultiComponentMixtureSOFC Declaration
\*---------------------------------------------------------------------------*/

class basicMultiComponentMixtureSOFC
{

protected:

    // Protected data
    
        //- Table of specie names
        speciesTable species_;
	
	//- Molar weight of the mixture
	volScalarField W_;
	
        //- Species mass fractions
        PtrList<volScalarField> y_;
	
	//- Species mass fractions
	PtrList<volScalarField> x_;
	
        //- Species total (convective + diffusive) mass fluxes
        PtrList<surfaceScalarField> n_;
	
	//- Species total (convective + diffusive) molar fluxes
	PtrList<surfaceScalarField> N_;
	
	//- Construct the mixture composition
        virtual void constructMixtureComposition(const fvMesh&) = 0;

public:

    // Constructors

        //- Construct from dictionary and mesh
        basicMultiComponentMixtureSOFC
        (
            const dictionary&,
            const wordList& specieNames,
            const fvMesh&
        );


    //- Destructor
    virtual ~basicMultiComponentMixtureSOFC()
    {}


    // Member functions

        //- Return the table of species
        const speciesTable& species() const
        {
            return species_;
        }

        //- Return the mass-fraction fields
        inline PtrList<volScalarField>& y();

        //- Return the const mass-fraction fields
        inline const PtrList<volScalarField>& y() const;

        //- Return the mass-fraction field for a specie given by index
        inline volScalarField& y(const label i);

        //- Return the const mass-fraction field for a specie given by index
        inline const volScalarField& y(const label i) const;

        //- Return the mass-fraction field for a specie given by name
        inline volScalarField& y(const word& specieName);

        //- Return the const mass-fraction field for a specie given by name
        inline const volScalarField& y(const word& specieName) const;

        //- Return the molar-fraction fields
        inline PtrList<volScalarField>& x();

        //- Return the const molar-fraction fields
        inline const PtrList<volScalarField>& x() const;

        //- Return the molar-fraction field for a specie given by index
        inline volScalarField& x(const label i);

        //- Return the const molar-fraction field for a specie given by index
        inline const volScalarField& x(const label i) const;

        //- Return the molar-fraction field for a specie given by name
        inline volScalarField& x(const word& specieName);

        //- Return the const molar-fraction field for a specie given by name
        inline const volScalarField& x(const word& specieName) const;

        //- Return the total mass flux fields
        inline PtrList<surfaceScalarField>& n();

        //- Return the const total mass flux fields
        inline const PtrList<surfaceScalarField>& n() const;

        //- Return the total mass flux field for a specie given by index
        inline surfaceScalarField& n(const label i);

        //- Return the const total mass flux field for a specie given by index
        inline const surfaceScalarField& n(const label i) const;

        //- Return the total mass flux field for a specie given by name
        inline surfaceScalarField& n(const word& specieName);

        //- Return the const total mass flux field for a specie given by name
        inline const surfaceScalarField& n(const word& specieName) const;

        //- Return the total molar flux fields
        inline PtrList<surfaceScalarField>& N();

        //- Return the const total molar flux fields
        inline const PtrList<surfaceScalarField>& N() const;

        //- Return the total molar flux field for a specie given by index
        inline surfaceScalarField& N(const label i);

        //- Return the const total molar flux field for a specie given by index
        inline const surfaceScalarField& N(const label i) const;

        //- Return the total molar flux field for a specie given by name
        inline surfaceScalarField& N(const word& specieName);

        //- Return the const total molar flux field for a specie given by name
        inline const surfaceScalarField& N(const word& specieName) const;

        //- Return the molar-weight field for the mixture
        inline volScalarField& W();

        //- Return the const molar-weight field for the mixture
        inline const volScalarField& W() const;	

        //- Does the mixture include this specie?
        inline bool contains(const word& specieName) const;

        inline scalar fres(const scalar ft, const scalar stoicRatio) const;

        inline tmp<volScalarField> fres
        (
            const volScalarField& ft,
            const dimensionedScalar& stoicRatio
        ) const;
	
	void write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "basicMultiComponentMixtureSOFCI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
