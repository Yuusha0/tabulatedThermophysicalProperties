/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018-2019 Yuusha and cbunge
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

\*---------------------------------------------------------------------------*/

#include "rhoThermo.H"
#include "makeThermo.H"

#include "specie.H"
#include "perfectGas.H"
#include "tabularEOS.H"
#include "incompressiblePerfectGas.H"
#include "Boussinesq.H"
#include "rhoConst.H"
#include "perfectFluid.H"
#include "PengRobinsonGas.H"
#include "adiabaticPerfectFluid.H"

#include "hConstThermo.H"
#include "janafThermo.H"
#include "hTabularThermo.H"
#include "sensibleEnthalpy.H"
#include "sensibleInternalEnergy.H"
#include "thermo.H"

#include "constTransport.H"
#include "sutherlandTransport.H"
#include "tabularTransport.H"

#include "icoPolynomial.H"
#include "hPolynomialThermo.H"
#include "polynomialTransport.H"

#include "hePsiThermo.H"
#include "heRhoThermo.H"
#include "heRhoTabularThermo.H"
#include "pureMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * Internal-energy-based * * * * * * * * * * * * * */
// Some simple combinations for testing purposes
makeThermos
(
    rhoThermo,
    heRhoTabularThermo,
    pureMixture,
    constTransport,
    sensibleInternalEnergy,
    hTabularThermo,
    perfectGas,
    specie
);

makeThermos
(
    rhoThermo,
    heRhoTabularThermo,
    pureMixture,
    constTransport,
    sensibleInternalEnergy,
    hConstThermo,
    tabularEOS,
    specie
);

makeThermos
(
    rhoThermo,
    heRhoTabularThermo,
    pureMixture,
    tabularTransport,
    sensibleInternalEnergy,
    hConstThermo,
    tabularEOS,
    specie
);

// Usually used for plasmas, using tabulated data
// These use two tables for h(p, T) and T(p, h)
makeThermos
(
    rhoThermo,
    heRhoTabularThermo,
    pureMixture,
    tabularTransport,
    sensibleInternalEnergy,
    hTabularThermo,
    tabularEOS,
    specie
);

makeThermos
(
    rhoThermo,
    heRhoTabularThermo,
    pureMixture,
    tabularTransport,
    sensibleEnthalpy,
    hTabularThermo,
    tabularEOS,
    specie
);

// This one only uses one table for h(p,T) and calculates T(p, h) manually.
makeThermos
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    tabularTransport,
    sensibleInternalEnergy,
    hTabularThermo,
    tabularEOS,
    specie
);

    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
