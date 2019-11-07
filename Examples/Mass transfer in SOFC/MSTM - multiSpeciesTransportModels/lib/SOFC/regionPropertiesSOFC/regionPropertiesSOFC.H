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
    regionPropertiesSOFC

Description
    Simple class to hold region information for coupled region simulations

SourceFiles
    regionPropertiesSOFC.C

\*---------------------------------------------------------------------------*/

#ifndef regionPropertiesSOFC_H
#define regionPropertiesSOFC_H

#include "IOdictionary.H"
#include "Time.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  
/*---------------------------------------------------------------------------*\
                       Class regionPropertiesSOFC Declaration
\*---------------------------------------------------------------------------*/

class regionPropertiesSOFC
:
    public IOdictionary
{
  
    // Private data

        //- List of the fluid region names
        List<word> regionNames_;

        //- List of Meshes
        PtrList<fvMesh> mesh_;

        //- Number of elements
        label nElements_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        regionPropertiesSOFC(const regionPropertiesSOFC&);

        //- Disallow default bitwise assignment
        void operator=(const regionPropertiesSOFC&);


public:

    // Constructors

        //- Construct from components
        regionPropertiesSOFC(const Time& runTime);


    // Destructor

        ~regionPropertiesSOFC()
        {}


    // Member Functions

        // Access

            //- Return const reference to the list of fluid region names
            const List<word>& regionNames() const;

            //- Return the region index
            label regionIndex(word regionName);

            //- Return the region name
            word& regionNames(label i);

            //- Return the region mesh
            fvMesh& mesh(label i);

            //- Return the region mesh
            fvMesh& mesh(word regionName);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
