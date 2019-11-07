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

#include "diffusivityModel.H"
#include "volFields.H"

namespace Foam
{
    defineTypeNameAndDebug(diffusivityModel, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diffusivityModel::diffusivityModel
(
    const volScalarField& p,
    const volScalarField& T
)
:
    transportDic_
    (
        IOobject
        (
            "transportProperties",
            T.mesh().time().constant(),
            T.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    
    thermoDic_
    (
        IOobject
        (
            "thermophysicalProperties",
            T.mesh().time().constant(),
            T.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    p_(p),
    
    T_(T),
    
    eps_(transportDic_.lookupOrDefault<scalar>("porosity", 1)),

    tau_(transportDic_.lookupOrDefault<scalar>("tortuosity", 1)),
    
    species_(thermoDic_.lookup("species"))
{
    DijModels_.setSize(0.5*species_.size()*(species_.size()+1));
    Dij_.setSize(0.5*species_.size()*(species_.size()+1));
        
    for(label i=0; i < species_.size(); i++)
    {
        for(label j=i; j < species_.size(); j++)
        {
            label k = species_.size()*i+j-0.5*i*(i+1);
            
            DijModels_.set
            (
                k,
                binaryDiffusivityModel::New
                (
                    species_[i],
                    species_[j],
                    transportDic_,
                    thermoDic_,
                    p,
                    T
                )
            );

            Dij_.set
            (
                k,
                new volScalarField
                (
                   eps_/tau_*DijModels_[k].D() 
                )
            );   
        }
    } 
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::diffusivityModel::update()
{    
    for(label i=0; i < species_.size(); i++)
    {
        for(label j=i; j < species_.size(); j++)
        {
            label k = species_.size()*i+j-0.5*i*(i+1);
            Dij_[k] = eps_/tau_*DijModels_[k].D();
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
