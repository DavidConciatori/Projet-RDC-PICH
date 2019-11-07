/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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

#include "material.H"
#include "error.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(material, 0);
defineRunTimeSelectionTable(material,);
defineRunTimeSelectionTable(material, Istream);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

autoPtr<material> material::New(word name)

{
    if (debug)
    {
        Info<< "material::New(Istream&): "
            << "constructing material"
            << endl;
    }

    word materialName(name);

        ConstructorTable::iterator cstrIter =
            ConstructorTablePtr_->find(materialName);

        if (cstrIter == ConstructorTablePtr_->end())
        {
            FatalErrorIn("material::New(Istream&)")
                << "Unknown material type " << materialName << nl << nl
                << "Valid material types are :" << endl
                << ConstructorTablePtr_->toc()
                << exit(FatalError);
        }

	return autoPtr<material>(cstrIter()());

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


  // Return the T-interpolation for thermal conductivity
  tmp<volScalarField> material::sigma(volScalarField& T)
  {
      const fvMesh& mesh = T.mesh();

      tmp<volScalarField> tSigma
      (
	  new volScalarField
	  (
	      IOobject
	      (
		  "lambda",
		  mesh.time().timeName(),
		  mesh,
		  IOobject::NO_READ,
		  IOobject::NO_WRITE
	      ),
	      mesh,
	      dimensionSet(-1, -3, 3, 0, 0, 2, 0)
	  )
      );
	      
      volScalarField& sigma = tSigma();
      forAll(T, celli)
      {
	  sigma[celli] = this->sigma(T[celli]);
      }
	      
      forAll(T.boundaryField(), patchi)
      {
	  forAll(T.boundaryField()[patchi], celli)
	  {
	      sigma.boundaryField()[patchi][celli] = this->sigma(T.boundaryField()[patchi][celli]);
	  }
      }

      return tSigma;
  }







} // End namespace Foam

// ************************************************************************* //
