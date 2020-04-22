/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018-2020 Yuusha and tilasoldo
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of tabulatedThermophysicalProperties on OpenFOAM.
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

Class
    Foam::tabularEOS

Description
    Tabular form of the equation of state

SourceFiles
    tabularEOSI.H
    tabularEOS.C

\*---------------------------------------------------------------------------*/

#ifndef tabularEOS_H
#define tabularEOS_H

#include "autoPtr.H"
#include "extrapolation2DTable.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie>
class tabularEOS;

template<class Specie>
tabularEOS<Specie> operator+
(
    const tabularEOS<Specie>&,
    const tabularEOS<Specie>&
);

template<class Specie>
tabularEOS<Specie> operator-
(
    const tabularEOS<Specie>&,
    const tabularEOS<Specie>&
);

template<class Specie>
tabularEOS<Specie> operator*
(
    const scalar,
    const tabularEOS<Specie>&
);

template<class Specie>
tabularEOS<Specie> operator==
(
    const tabularEOS<Specie>&,
    const tabularEOS<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const tabularEOS<Specie>&
);


/*---------------------------------------------------------------------------*\
                        Class tabularEOS Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class tabularEOS
:
    public Specie
{
    // Private data

        //- "Perfect gaz" model for Cp - Cv
        Switch perfectGas_;

        //- Density table
        extrapolation2DTable<scalar> densityTable;


public:

    // Constructors

        //- Construct from components
        inline tabularEOS
	(
	    const Specie& sp,
	    const Switch& perfectGas,
	    const extrapolation2DTable<scalar>& densityTable
	);

        //- Construct from Istream
        tabularEOS(Istream&);

        //- Construct from dictionary
        tabularEOS(const dictionary& dict);

        //- Construct as named copy
        inline tabularEOS(const word& name, const tabularEOS&);

        //- Construct and return a clone
        inline autoPtr<tabularEOS> clone() const;

        // Selector from Istream
        inline static autoPtr<tabularEOS> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "tabularEOS<" + word(Specie::typeName_()) + '>';
        }


        // Fundamental properties

            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = false;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

            //- Return density [kg/m^3]
            inline scalar rho(scalar p, scalar T) const;

            //- Return enthalpy departure [J/kg]
            inline scalar H(scalar p, scalar T) const;

            //- Return Cp departure [J/(kg K]
            inline scalar Cp(scalar p, scalar T) const;

            //- Return entropy [J/(kg K)]
            inline scalar S(const scalar p, const scalar T) const;

            //- Return compressibility rho/p [s^2/m^2]
            inline scalar psi(scalar p, scalar T) const;

            //- Return compression factor []
            inline scalar Z(scalar p, scalar T) const;

            //- Return (cp - cv) [J/(kmol K]
            inline scalar CpMCv(scalar p, scalar T) const;


        // IO

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator+=(const tabularEOS&);
        inline void operator-=(const tabularEOS&);
        inline void operator*=(const scalar);


    // Friend operators

        friend tabularEOS operator+ <Specie>
        (
            const tabularEOS&,
            const tabularEOS&
        );

        friend tabularEOS operator- <Specie>
        (
            const tabularEOS&,
            const tabularEOS&
        );

        friend tabularEOS operator* <Specie>
        (
            const scalar s,
            const tabularEOS&
        );

        friend tabularEOS operator== <Specie>
        (
            const tabularEOS&,
            const tabularEOS&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const tabularEOS&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tabularEOSI.H"

#ifdef NoRepository
#   include "tabularEOS.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
