/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018-2020 Yuusha and cbunge
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
    Foam::rhoTabularThermo

Description
    Basic thermodynamic properties based on density

SourceFiles
    rhoTabularThermo.C

\*---------------------------------------------------------------------------*/

#ifndef rhoTabularThermo_H
#define rhoTabularThermo_H

#include "fluidThermo.H"
#include "rhoThermo.H"
#include "runTimeSelectionTables.H"
#include "extrapolation2DTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class rhoTabularThermo Declaration
\*---------------------------------------------------------------------------*/

class rhoTabularThermo
:
    public fluidThermo
{

protected:

    // Protected data
        //- Density field [kg/m^3]
        volScalarField rho_;

        //- Compressibility [s^2/m^2]
       volScalarField psi_;

        //- Dynamic viscosity [kg/m/s]
        volScalarField mu_;

        //- Density [kg/m^3]
        extrapolation2DTable<scalar> densityTable;

    // Protected Member Functions

        //- Construct as copy (not implemented)
        rhoTabularThermo(const rhoTabularThermo&);


public:
    //- Runtime type information
    TypeName("rhoTabularThermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        rhoTabularThermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );


    // Constructors
        //- Construct from mesh and phase name
        rhoTabularThermo
        (
	    const fvMesh&,
            const word& phaseName
	);
        rhoTabularThermo
        (
            const fvMesh&,
            const dictionary&,
            const word& phaseName
        );

    //- Selector
    static autoPtr<rhoTabularThermo> New
    (
        const fvMesh& mesh,
        const word& phaseName=word::null
    );


    //- Destructor
    virtual ~rhoTabularThermo();


    // Member functions

        volScalarField& lookupOrConstruct2
        (
            const fvMesh& mesh,
            const char* name,
	    dimensionSet units
        ) const;

        // Fields derived from thermodynamic state variables

            //- Density [kg/m^3] - uses current value of pressure
            virtual tmp<volScalarField> rho() const;

            //- Add the given density correction to the density field.
            //  Used to update the density field following pressure solution
            virtual void correctRho(const volScalarField& deltaRho);

            //- Compressibility [s^2/m^2]
            virtual const volScalarField& psi() const;


        // Access to transport state variables

            //- Dynamic viscosity of mixture [kg/m/s]
            virtual tmp<volScalarField> mu() const;

            //- Dynamic viscosity of mixture for patch [kg/m/s]
            virtual tmp<scalarField> mu(const label patchi) const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
