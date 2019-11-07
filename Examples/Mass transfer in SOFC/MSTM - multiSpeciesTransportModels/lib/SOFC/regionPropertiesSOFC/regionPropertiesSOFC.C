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

\*---------------------------------------------------------------------------*/

#include "regionPropertiesSOFC.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionPropertiesSOFC::regionPropertiesSOFC(const Time& runTime)
:
    IOdictionary
    (
        IOobject
        (
            "regionPropertiesSOFC",
            "modelInput",
            runTime.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    regionNames_(lookup("regionNames"))
{
    nElements_ = regionNames_.size();
    mesh_.setSize(nElements_);


    for(label i=0; i<nElements_; i++)
    {
        mesh_.set 
        (
            i,
            new fvMesh
            (
                IOobject
                (
                    regionNames_[i],
                    runTime.timeName(),
                    runTime,
                    IOobject::MUST_READ
                )
            )
        ); 
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::List<Foam::word>& Foam::regionPropertiesSOFC::regionNames() const
{
    return regionNames_;
}


Foam::label Foam::regionPropertiesSOFC::regionIndex(word regionName)
{
  forAll (regionNames_, regionI)
  {
      if (regionNames_[regionI] == regionName)
      {
          return (regionI);
      }
  }
  
  Info << "\nFatal Error!!!\n" <<
  "The region " << regionName << " doesen't exist\n";
  Foam::FatalError.exit();
  
  return (0);
}


Foam::word& Foam::regionPropertiesSOFC::regionNames(label i)
{
    if (i > nElements_ || i < 1)
    {
        Info << "\nFatal Error!!!\n" <<
        "The region number " << i << " doesen't exist\n";
        Foam::FatalError.exit();
        return (regionNames_[0]);
    }
    
    i--;
    
    return regionNames_[i];
}


Foam::fvMesh& Foam::regionPropertiesSOFC::mesh(label i)
{
    if (i > nElements_ || i < 1)
    {
        Info << "\nFatal Error!!!\n" <<
        "The region number " << i << " doesen't exist\n";
        Foam::FatalError.exit();
        return (mesh_[0]);
    }
    
    i--;
    
    return(mesh_[i]);
}


Foam::fvMesh& Foam::regionPropertiesSOFC::mesh(word regionName)
{
    forAll (regionNames_, regionI)
    {
        if (regionNames_[regionI] == regionName)
        {
            return (mesh_[regionI]);
        }
    }
    
    Info << "\nFatal Error!!!\n" <<
    "The region " << regionName << " doesen't exist\n";
    Foam::FatalError.exit();
    
    return (mesh_[0]);
}


// ************************************************************************* //
