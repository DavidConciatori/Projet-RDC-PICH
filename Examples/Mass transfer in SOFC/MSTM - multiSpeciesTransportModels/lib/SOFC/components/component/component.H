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
    Foam::SOFC::component

Description
    A base class for all component (channel, electrode, electrolite) of
    solid oxide fuel cell

SourceFiles
    component.C

\*---------------------------------------------------------------------------*/

#ifndef component_H
#define component_H

#include "fvCFD.H"

#include "regionPropertiesSOFC.H"

#include "fixedFluxFvPatchField.H"

#include "directMappedPatchBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace SOFC
{
  
/*---------------------------------------------------------------------------*\
                           Class component Declaration
\*---------------------------------------------------------------------------*/

class component
:

    public IOdictionary
{
  
    // Private data
      
        //- Mesh
        const fvMesh& mesh_;        
      
        //- Component name
        word name_;


    // Private member functions

        //- Check if in the domain is present the specified boundary
        void check
        (
            label& fromIndex,
            label& toIntex,
            word& fromPatchName,
            word& toPatchName
        );

      
public:

    // Public data
        
        //- Flag that activate the component
        bool activated_;

        
    //- Constructor
        
        //- Construct from components
        component
        (
            const word elementName,
            regionPropertiesSOFC& regions
        );

        
    //- Destructor
        
        virtual ~component() {}


    // Member Functions

        //- Return the name of the component
        word name() const
        {
            return (name_);
        }

        //- Return const acces to the mesh
        const fvMesh& mesh() const
        {
            return (mesh_);
        }
                
        //- Map the boundary values from values of neighbour component
        template<typename fieldType>
        void mapBoundaryField(fieldType& field, word patchName)
        {                   
            int coeff = 1;;
            if(field.typeName == "surfaceScalarField")
            {
                coeff = -1;
            }
          
                  label index = mesh().boundaryMesh().findPatchID(patchName);
          
            // Get the coupling information from the directMappedPatchBase
            const directMappedPatchBase& mpp =
            refCast<const directMappedPatchBase>
            (
                mesh().boundaryMesh()[index]
            );

            // Force recalculation of mapping and schedule
            const mapDistribute& distMap = mpp.map();
            (void)distMap.schedule();
            
            const fvMesh& nbrMesh = refCast<const fvMesh>(mpp.sampleMesh());
            const word& fldName = field.name();

            // Result of obtaining remote values
            scalarField newValues;

            switch (mpp.mode())
            {
                case directMappedPatchBase::NEARESTPATCHFACE:
                {
                    const label nbrPatchID = nbrMesh.boundaryMesh().findPatchID
                    (
                        mpp.samplePatch()
                    );
                    if (nbrPatchID < 0)
                    {
                        FatalErrorIn
                        (
                            "void directMappedFixedValueFvPatchField<Type>::"
                            "updateCoeffs()"
                        )<< "Unable to find sample patch " << mpp.samplePatch()
                        << " in region " << mpp.sampleRegion()
                        << " for patch " << mesh().boundaryMesh()[index].name() 
                        << nl << abort(FatalError);
                    }

                    const fieldType& nbrField = nbrMesh.lookupObject<fieldType>
                    (
                        fldName
                    );
                    
                    if(coeff)
                    {
                      newValues = nbrField.boundaryField()[nbrPatchID];
                    }
                    else
                    {
                      newValues = nbrField.boundaryField()[nbrPatchID]
                          / nbrMesh.magSf().boundaryField()[nbrPatchID];
                    }
                  
                    mapDistribute::distribute
                    (
                        Pstream::defaultCommsType,
                        distMap.schedule(),
                        distMap.constructSize(),
                        distMap.subMap(),
                        distMap.constructMap(),
                        newValues
                    );

                    break;
                }
                default:
                {
                    FatalErrorIn
                    (
                        "SOFC::component::mapField(fieldType&, word&)"
                    )<< "Only sampling mode NEARESTPATCHFACE is implemented"
                    << nl << abort(FatalError);
                }
            }          
          
            const scalarField& area = mesh().magSf().boundaryField()[index];
            scalarField& boundaryValues = field.boundaryField()[index];

            if(coeff)
            {
              forAll(boundaryValues, faceI)
              {
                boundaryValues[faceI] = coeff * newValues[faceI];
              }
            }
            else
            {
              forAll(boundaryValues, faceI)
              {
                boundaryValues[faceI] = coeff * newValues[faceI] * area[faceI];
              }
            }
        }   
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace SOFC

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
