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

Typedefs
    Foam::tabularThermoPhysicsTypes

Description
    Type definitions for tabular thermo-physics models

\*---------------------------------------------------------------------------*/

#ifndef tabularThermoPhysicsTypes_H
#define tabularThermoPhysicsTypes_H

#include "specie.H"

#include "sensibleEnthalpy.H"
#include "sensibleInternalEnergy.H"
#include "thermo.H"

#include "tabularEOS.H"
#include "hTabularThermo.H"
#include "tabularTransport.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // thermo physics types based on sensibleEnthalpy
    typedef
    tabularTransport
    <
        species::thermo
        <
            hTabularThermo
            <
                tabularEOS<specie>
            >,
            sensibleEnthalpy
        >
    > tabularHThermoPhysics;


    // thermo physics types based on sensibleInternalEnergy
    typedef
    tabularTransport
    <
        species::thermo
        <
            hTabularThermo
            <
                tabularEOS<specie>
            >,
            sensibleInternalEnergy
        >
    > tabularEThermoPhysics;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
