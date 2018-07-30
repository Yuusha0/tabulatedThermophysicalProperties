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

#include "hTabularThermo.H"
#include "IOstreams.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::hTabularThermo<EquationOfState>::hTabularThermo
(
    Istream& is
)
:
    EquationOfState(is),
    mode_(hTabularThermo::tabulated),
    Hf_(0),
    Sf_(0)
{

    cpTable = extrapolation2DTable<scalar>("constant/cpTable");
    hTable = extrapolation2DTable<scalar>("constant/hTable");
    hfTable = extrapolation2DTable<scalar>("constant/hfTable");
    cpTable.outOfBounds(extrapolation2DTable<scalar>::EXTRAPOLATE);
    hTable.outOfBounds(extrapolation2DTable<scalar>::EXTRAPOLATE);
    hfTable.outOfBounds(extrapolation2DTable<scalar>::EXTRAPOLATE);
}


template<class EquationOfState>
Foam::hTabularThermo<EquationOfState>::hTabularThermo
(
    const dictionary& dict
)
:
    EquationOfState(dict),
    mode_
    (
	wordToHfMode
	(
	    dict.subDict("thermodynamics").subDict("hf").lookupOrDefault<word>
	    (
		"mode","constant"
	    )
	)
    ),
    cpTable(dict.subDict("thermodynamics").subDict("Cp")),
    hTable(dict.subDict("thermodynamics").subDict("h"))
{

    switch(mode_)
    {
    case hTabularThermo::constant:
	{
	    // Create constant enthalpy of formation
	    Hf_ = readScalar
		(
		    dict.subDict("thermodynamics").subDict("hf").lookup("Hf")
		);

            // Create empty table for enthalpy of formation
	    hfTable = extrapolation2DTable<scalar>();

	    break;
	}

    case hTabularThermo::tabulated:
       {
	   hfTable =
	       extrapolation2DTable<scalar>
	       (
		   dict.subDict("thermodynamics").subDict("hf")
	       );

	   Hf_ = 0;

	   break;
       }
    }

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
void Foam::hTabularThermo<EquationOfState>::write
(
    Ostream& os
) const
{
    EquationOfState::write(os);

    dictionary dict("thermodynamics");
    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const hTabularThermo<EquationOfState>& pt
)
{
    os  << static_cast<const EquationOfState&>(pt) << tab;

    os.check
    (
        "operator<<"
        "("
            "Ostream&, "
            "const hTabularThermo<EquationOfState>&"
        ")"
    );

    return os;
}


// ************************************************************************* //
