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
    Foam::SOFC::ohmicOverpotential

Description
    Ohmic overpotential 
    
SourceFiles
    ohmicOverpotential.C

\*---------------------------------------------------------------------------*/

#ifndef ohmicOverpotential_H
#define ohmicOverpotential_H

#include "solidComponent.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  
namespace SOFC
{ 

/*---------------------------------------------------------------------------*\
                           Class ohmicOverpotential Declaration
\*---------------------------------------------------------------------------*/

class ohmicOverpotential

{
    // Private data
    
        //- Layers thickness
        scalarList thickness_;
        
        //- Total resistivity
        scalarField totalSigma_;
        
        //- Anode component
        solidComponent& anode_;
        
        //- Anode electrolyte
        solidComponent& electrolyte_;
        
        //- Anode cathode
        solidComponent& cathode_;

        
public:
  
      //- Runtime type information
      TypeName("ohmicOverpotential");


    // Constructors

        //- Construct from components
        ohmicOverpotential(solidComponent& elec, solidComponent& anode, solidComponent& cathode);
         

        // Destructor
        virtual ~ohmicOverpotential()
        {}


    // Member Functions
    
        //- Return the current density (point value)
        inline scalar eta(scalar& i, int index) const
        { 
            return i*totalSigma_[index];  
        }
        
        //- Return the current density (global value)
        scalarField eta(scalarField& i) const
        { 
            return i*totalSigma_;  
        }
        
        //- Update the overpotential
        void update(scalarField& T);
        
        //- Return the thickness of a layer specified by index
        //  [0] anode, [1] electrolyte, [2] cathode
        scalar thickness(int i)
        { 
            return thickness_[i];  
        }
        
        //- Return the thickness of a layer specified by name
        scalar thickness(solidComponent& component);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace SOFC

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
