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

inline Foam::PtrList<Foam::volScalarField>&
Foam::basicMultiComponentMixtureSOFC::y()
{
    return y_;
}


inline const Foam::PtrList<Foam::volScalarField>&
Foam::basicMultiComponentMixtureSOFC::y() const
{
    return y_;
}


inline Foam::volScalarField& Foam::basicMultiComponentMixtureSOFC::y
(
    const label i
)
{
    return y_[i];
}


inline const Foam::volScalarField&
Foam::basicMultiComponentMixtureSOFC::y
(
    const label i
) const
{
    return y_[i];
}


inline Foam::volScalarField& Foam::basicMultiComponentMixtureSOFC::y
(
    const word& specieName
)
{
    return y_[species_[specieName]];
}


inline const Foam::volScalarField&
Foam::basicMultiComponentMixtureSOFC::y
(
    const word& specieName
) const
{
    return y_[species_[specieName]];
}


inline Foam::PtrList<Foam::volScalarField>&
Foam::basicMultiComponentMixtureSOFC::x()
{
    return x_;
}


inline const Foam::PtrList<Foam::volScalarField>&
Foam::basicMultiComponentMixtureSOFC::x() const
{
    return x_;
}


inline Foam::volScalarField& Foam::basicMultiComponentMixtureSOFC::x
(
    const label i
)
{
    return x_[i];
}


inline const Foam::volScalarField& Foam::basicMultiComponentMixtureSOFC::x
(
    const label i
) const
{
    return x_[i];
}


inline Foam::volScalarField& Foam::basicMultiComponentMixtureSOFC::x
(
    const word& specieName
)
{
    return x_[species_[specieName]];
}


inline const Foam::volScalarField& Foam::basicMultiComponentMixtureSOFC::x
(
    const word& specieName
) const
{
    return x_[species_[specieName]];
}


inline Foam::PtrList<Foam::surfaceScalarField>&
Foam::basicMultiComponentMixtureSOFC::n()
{
    return n_;
}


inline const Foam::PtrList<Foam::surfaceScalarField>&
Foam::basicMultiComponentMixtureSOFC::n() const
{
    return n_;
}


inline Foam::surfaceScalarField& Foam::basicMultiComponentMixtureSOFC::n
(
    const label i
)
{
    return n_[i];
}


inline const Foam::surfaceScalarField& Foam::basicMultiComponentMixtureSOFC::n
(
    const label i
) const
{
    return n_[i];
}


inline Foam::surfaceScalarField& Foam::basicMultiComponentMixtureSOFC::n
(
    const word& specieName
)
{
    return n_[species_[specieName]];
}


inline const Foam::surfaceScalarField& Foam::basicMultiComponentMixtureSOFC::n
(
    const word& specieName
) const
{
    return n_[species_[specieName]];
}


inline Foam::PtrList<Foam::surfaceScalarField>&
Foam::basicMultiComponentMixtureSOFC::N()
{
    return N_;
}


inline const Foam::PtrList<Foam::surfaceScalarField>&
Foam::basicMultiComponentMixtureSOFC::N() const
{
    return N_;
}


inline Foam::surfaceScalarField& Foam::basicMultiComponentMixtureSOFC::N
(
    const label i
)
{
    return N_[i];
}


inline const Foam::surfaceScalarField& Foam::basicMultiComponentMixtureSOFC::N
(
    const label i
) const
{
    return N_[i];
}


inline Foam::surfaceScalarField& Foam::basicMultiComponentMixtureSOFC::N
(
    const word& specieName
)
{
    return N_[species_[specieName]];
}


inline const Foam::surfaceScalarField& Foam::basicMultiComponentMixtureSOFC::N
(
    const word& specieName
) const
{
    return N_[species_[specieName]];
}


inline Foam::volScalarField& Foam::basicMultiComponentMixtureSOFC::W()
{
    return W_;
}


inline const Foam::volScalarField&
Foam::basicMultiComponentMixtureSOFC::W() const
{
    return W_;
}


inline bool Foam::basicMultiComponentMixtureSOFC::contains
(
    const word& specieName
) const
{
    return species_.contains(specieName);
}


inline Foam::scalar Foam::basicMultiComponentMixtureSOFC::fres
(
    const scalar ft,
    const scalar stoicRatio
) const
{
    return max(ft - (1.0 - ft)/stoicRatio, 0.0);
}


inline Foam::tmp<Foam::volScalarField>
Foam::basicMultiComponentMixtureSOFC::fres
(
    const volScalarField& ft,
    const dimensionedScalar& stoicRatio
) const
{
    return max(ft - (scalar(1) - ft)/stoicRatio.value(), scalar(0));
}


// ************************************************************************* //
