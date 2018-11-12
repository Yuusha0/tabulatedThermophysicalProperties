/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 Yuusha and tilasoldo
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of tilasoldo and Yuusha contribution to OpenFOAM.
    It is based on chriss85 contribution for OpenFOAM 2.3.x.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "tabularEOS.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie>
inline Foam::tabularEOS<Specie>::tabularEOS
(
    const Specie& sp,
    const extrapolation2DTable<scalar>& densityTable
)
:
    Specie(sp),
    densityTable(densityTable)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::tabularEOS<Specie>::tabularEOS
(
    const word& name,
    const tabularEOS<Specie>& ico
)
:
    Specie(name, ico),
    densityTable(ico.densityTable)
{}


template<class Specie>
inline Foam::autoPtr<Foam::tabularEOS<Specie> >
Foam::tabularEOS<Specie>::clone() const
{
    return autoPtr<tabularEOS<Specie>>
    (
       new tabularEOS<Specie>(*this)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::tabularEOS<Specie> >
Foam::tabularEOS<Specie>::New(const dictionary& dict)
{
    return autoPtr<tabularEOS<Specie>>
    (
        new tabularEOS<Specie>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::tabularEOS<Specie>::rho
(
    scalar p,
    scalar T
) const
{
    return densityTable(T, p);
}

// Totally not sure of what I am doing here...
template<class Specie>
inline Foam::scalar Foam::tabularEOS<Specie>::H
(
    scalar p,
    scalar T
) const
{
    return 0.;
}

template<class Specie>
inline Foam::scalar Foam::tabularEOS<Specie>::Cp
(
    scalar p,
    scalar T
) const
{
    return 0.;
}

template<class Specie>
inline Foam::scalar Foam::tabularEOS<Specie>::S
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::tabularEOS<Specie>::psi
(
    scalar p,
    scalar T
) const
{
    return rho(p, T) / (p + SMALL);
}


template<class Specie>
inline Foam::scalar Foam::tabularEOS<Specie>::Z
(
    scalar p,
    scalar T
) const
{
    return 0.0;
}


template<class Specie>
inline Foam::scalar Foam::tabularEOS<Specie>::CpMCv
(
    scalar p,
    scalar T
) const
{
    return -(p/sqr(densityTable(T,p)))*densityTable.Tderivative(T,p);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::tabularEOS<Specie>::operator+=(const tabularEOS<Specie>& ico)
{
    Specie::operator+=(ico);

    //scalar molr1 = this->nMoles();

    //Specie::operator+=(ico);
/*
    molr1 /= this->nMoles();
    scalar molr2 = ico.nMoles()/this->nMoles();

    rho_ = molr1*rho_ + molr2*ico.rho_;
*/
}


template<class Specie>
inline void Foam::tabularEOS<Specie>::operator-=(const tabularEOS<Specie>& ico)
{
    Specie::operator-=(ico);

/*
    scalar molr1 = this->nMoles();
*/
    //Specie::operator-=(ico);
/*
    molr1 /= this->nMoles();
    scalar molr2 = ico.nMoles()/this->nMoles();

    rho_ = molr1*rho_ - molr2*ico.rho_;
*/
}


template<class Specie>
inline void Foam::tabularEOS<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);

    //Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
Foam::tabularEOS<Specie> Foam::operator+
(
    const tabularEOS<Specie>& ico1,
    const tabularEOS<Specie>& ico2
)
{
    Specie sp
    (
        static_cast<const Specie&>(ico1)
      + static_cast<const Specie&>(ico2)
    );

    if (mag(sp.Y()) < SMALL)
    {
        return tabularEOS<Specie>
        (
            sp,
            ico1.densityTable
        );
    }
    else
    {
        const scalar Y1 = ico1.Y()/sp.Y();
        const scalar Y2 = ico2.Y()/sp.Y();

        return tabularEOS<Specie>
        (
            sp,
            Y1*ico1.densityTable + Y2*ico2.densityTable
        );
    }
}


template<class Specie>
Foam::tabularEOS<Specie> Foam::operator-
(
    const tabularEOS<Specie>& ico1,
    const tabularEOS<Specie>& ico2
)
{
    Specie sp
    (
        static_cast<const Specie&>(ico1)
      - static_cast<const Specie&>(ico2)
    )
;

    if (mag(sp.Y()) < SMALL)
    {
        return tabularEOS<Specie>
        (
            sp,
            ico1.densityTable
        );
    }
    else
    {
        const scalar Y1 = ico1.Y()/sp.Y();
        const scalar Y2 = ico2.Y()/sp.Y();

        return tabularEOS<Specie>
        (
            sp,
            Y1*ico1.densityTable - Y2*ico2.densityTable
        );
    }

}


template<class Specie>
Foam::tabularEOS<Specie> Foam::operator*
(
    const scalar s,
    const tabularEOS<Specie>& ico
)
{
    return tabularEOS<Specie>
    (
        s*static_cast<const Specie&>(ico),
        ico.densityTable
    );
}


template<class Specie>
Foam::tabularEOS<Specie> Foam::operator==
(
    const tabularEOS<Specie>& ico1,
    const tabularEOS<Specie>& ico2
)
{
    Specie sp
    (
        static_cast<const Specie&>(ico1)
     == static_cast<const Specie&>(ico2)
    );

    const scalar Y1 = ico1.Y()/sp.Y();
    const scalar Y2 = ico2.Y()/sp.Y();

    return tabularEOS<Specie>
    (
        sp,
        Y2*ico2.densityTable - Y1*ico1.densityTable
    );
}


// ************************************************************************* //
