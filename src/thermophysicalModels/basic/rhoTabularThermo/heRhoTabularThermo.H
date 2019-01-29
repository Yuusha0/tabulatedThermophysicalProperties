/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018-2019 Yuusha, cbunge
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
    Foam::heRhoTabularThermo

Description
    Energy for a mixture based on density

SourceFiles
    heRhoTabularThermo.C

\*---------------------------------------------------------------------------*/

#ifndef heRhoTabularThermo_H
#define heRhoTabularThermo_H

#include "rhoTabularThermo.H"
#include "extrapolation2DTable.H"
#include "heThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class heRhoTabularThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTabularThermo, class MixtureType>
class heRhoTabularThermo
:
    public heThermo<BasicTabularThermo, MixtureType>
{
    // Private Member Functions

        //- Calculate the thermo variables
        void calculate();

        //- Construct as copy (not implemented)
        heRhoTabularThermo
	(
	    const heRhoTabularThermo<BasicTabularThermo, MixtureType>&
	);

public:

    //- Runtime type information
    TypeName("heRhoTabularThermo");


    // Constructors

        //- Construct from mesh and phase name
        heRhoTabularThermo
        (
            const fvMesh&,
            const word& phaseName
        );

    extrapolation2DTable<scalar> TTable;

    //- Destructor
    virtual ~heRhoTabularThermo();


    // Member functions

        //- Update properties
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
#   include "heRhoTabularThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
