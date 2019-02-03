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

#include "rhoTabularThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(rhoTabularThermo, 0);
    defineRunTimeSelectionTable(rhoTabularThermo, fvMesh);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rhoTabularThermo::rhoTabularThermo(const fvMesh& mesh, const word& phaseName)
:
    fluidThermo(mesh, phaseName),

    rho_
     (
         IOobject
         (
             phasePropertyName("thermo:rho"),
             mesh.time().timeName(),
             mesh,
             IOobject::NO_READ,
             IOobject::AUTO_WRITE
         ),
         mesh,
         dimDensity
     ),

    psi_
    (
        IOobject
        (
            phasePropertyName("thermo:psi"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(0, -2, 2, 0, 0)
    ),

    mu_
    (
        IOobject
        (
            phasePropertyName("thermo:mu"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(1, -1, -1, 0, 0)
    ),
    densityTable("constant/densityTable")
{
    densityTable.outOfBounds(extrapolation2DTable<scalar>::CLAMP);
}


Foam::rhoTabularThermo::rhoTabularThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    fluidThermo(mesh, dict, phaseName),
    rho_
    (
        IOobject
        (
            phasePropertyName("thermo:rho"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimDensity
    ),

    psi_
    (
        IOobject
        (
            phasePropertyName("thermo:psi"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(0, -2, 2, 0, 0)
    ),

    mu_
    (
        IOobject
        (
            phasePropertyName("thermo:mu"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(1, -1, -1, 0, 0)
    )
{}
// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::rhoTabularThermo> Foam::rhoTabularThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return basicThermo::New<rhoTabularThermo>(mesh, phaseName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::rhoTabularThermo::~rhoTabularThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::volScalarField& Foam::rhoTabularThermo::lookupOrConstruct2
(
    const fvMesh& mesh,
    const char* name,
    dimensionSet units
) const
{
    if (!mesh.objectRegistry::foundObject<volScalarField>(name))
    {
        volScalarField* fPtr
        (
            new volScalarField
            (
                IOobject
                (
                    name,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
		dimensionedScalar("tmp", units, 0)
            )
        );

        // Transfer ownership of this object to the objectRegistry
        fPtr->store(fPtr);
    }

    return const_cast<volScalarField&>
    (
        mesh.objectRegistry::lookupObject<volScalarField>(name)
    );
}

Foam::tmp<Foam::volScalarField> Foam::rhoTabularThermo::rho() const
{
    volScalarField rho_ =
	lookupOrConstruct2
	(
	    T_.mesh(),
	    phasePropertyName("thermo:rho_").c_str(),
	    dimensionSet(1, -1, -2, 0, 0)
	 );
    volScalarField T_back =
	lookupOrConstruct2
	(
	    T_.mesh(),
	    phasePropertyName("thermo:T_back").c_str(),
	    dimensionSet(0, 0, 0, 1, 0)
	 );

    volScalarField p_back =
	lookupOrConstruct2
	(
	    T_.mesh(),
	    phasePropertyName("thermo:p_back").c_str(),
	    dimensionSet(1, -3, 0, 0, 0)
	 );
    /* Check if the temperature or pressure fields
       have changed since the last iteration. */
    forAll(T_, faceI)
    {
	if(T_[faceI] != T_back[faceI] || p_[faceI] != p_back[faceI])
	{
		//Info << "difference" << endl;
    		forAll(T_, faceI2)
		{
			rho_[faceI2] = densityTable(T_[faceI2], p_[faceI2]);
			T_back[faceI2] = T_[faceI2];
			p_back[faceI2] = p_[faceI2];
		}
		return rho_;
	}
    }
    return rho_;
}

void Foam::rhoTabularThermo::correctRho(const Foam::volScalarField& deltaRho)
{
    rho_ += deltaRho;
}

const Foam::volScalarField& Foam::rhoTabularThermo::psi() const
{
    return rho()/(p_ + dimensionedScalar("tmpstabs", p_.dimensions(), SMALL));
}

Foam::tmp<Foam::volScalarField> Foam::rhoTabularThermo::mu() const
{
    return mu_;
}


Foam::tmp<Foam::scalarField> Foam::rhoTabularThermo::mu(const label patchi) const
{
    return mu_.boundaryField()[patchi];
}


// ************************************************************************* //
