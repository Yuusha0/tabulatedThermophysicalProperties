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

Class
    Foam::extrapolation2DTable

Description
    2D table extrapolation. The data must be in ascending order in both
    dimensions x and y.

SourceFiles
    extrapolation2DTable.C

\*---------------------------------------------------------------------------*/

#ifndef extrapolation2DTable_H
#define extrapolation2DTable_H

#include "List.H"
#include "Tuple2.H"
#include "tableReader.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Type> class extrapolation2DTable;

template<class Type>
inline extrapolation2DTable<Type> operator+
(
    const extrapolation2DTable<Type>&,
    const extrapolation2DTable<Type>&
);


template<class Type>
inline extrapolation2DTable<Type> operator-
(
    const extrapolation2DTable<Type>&,
    const extrapolation2DTable<Type>&
);


template<class Type>
inline extrapolation2DTable<Type> operator*
(
    const scalar,
    const extrapolation2DTable<Type>&
);


/*---------------------------------------------------------------------------*\
                    Class extrapolation2DTable Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class extrapolation2DTable:
	public List<Tuple2<scalar, List<Tuple2<scalar, Type> > > >
{
public:

    // Public data types

        //- Enumeration for handling out-of-bound values
        enum boundsHandling
        {
            ERROR,         /*!< Exit with a FatalError */
            WARN,          /*!< Issue warning and extrapolate value (default) */
            EXTRAPOLATE    /*!< Extrapolate value to the start/end value */
        };

        //- Enumeration for search in the tables
        enum searchMethod
	{
	    simple,  /*!< The whole table's covered */
	    uniform, /*!< Using constant spacing of the table */
	    bisect   /*!< Bisection search */
	};

        //- Convenience typedef
        typedef List<Tuple2<scalar, List<Tuple2<scalar, Type> > > > table;


private:

    // Private data

        //- Enumeration for handling out-of-bound values
        boundsHandling boundsHandling_;

        //- Enumeration for search in the table
        searchMethod searchMethod_;

        //- File name
        fileName fileName_;

        //- the actual reader
        autoPtr<tableReader<Type> > reader_;

        //- true or yes if all the values are 0
        Switch isNull_;


    // Private Member Functions

        //- Read the table of data from file
        void readTable();

        //- Return extrapolated value in List
        Type extrapolateValue
        (
            const List<Tuple2<scalar, Type> >& data,
            const scalar
        ) const;

        //- Return an X index from the matrix
        template<class L>
        label Xi
        (
	    const List<Tuple2<scalar, L> >& t,
            const scalar valueX,
            const bool reverse
        ) const;


public:

    // Constructors

        //- Construct null
        extrapolation2DTable();

        //- Construct from components
        extrapolation2DTable
        (
            const List<Tuple2<scalar, List<Tuple2<scalar, Type> > > >& values,
            const boundsHandling bounds,
	    const searchMethod method,
            const fileName& fName,
	    const Switch isNull
        );

        //- Construct given the name of the file containing the table of data
        extrapolation2DTable(const fileName& fName);

        //- Construct by reading the fileName and boundsHandling from dictionary
        extrapolation2DTable(const dictionary& dict);

        //- Construct copy
        extrapolation2DTable(const extrapolation2DTable& extrapTable);


    // Member Functions

        //- Return the out-of-bounds handling as a word
        word boundsHandlingToWord(const boundsHandling& bound) const;

        //- Return the out-of-bounds handling as an enumeration
        boundsHandling wordToBoundsHandling(const word& bound) const;

        //- Set the out-of-bounds handling from enum, return previous setting
        boundsHandling outOfBounds(const boundsHandling& bound);

        //- Return the search method as a word
        word searchMethodToWord(const searchMethod& searchMethod) const;

        //- Return the search method as an enumeration
        searchMethod wordToSearchMethod(const word& searchMethod) const;

        //- Check that list is monotonically increasing
        //  Exit with a FatalError if there is a problem
        void checkOrder() const;

        //- Derivate a tabulated value.
        //  Simple derivation for testing purpose. Must be improved.
        Type Tderivative(const scalar valueX,const scalar valueY) const;

        //- Write
        void write(Ostream& os) const;


    // Member Operators

        //- Return an element of constant Tuple2<scalar, Type>
        const List<Tuple2<scalar, Type> >& operator[](const label) const;

        //- Return an extrapolated value
        Type operator()(const scalar, const scalar) const;

        inline extrapolation2DTable& operator=(const extrapolation2DTable&);


   // Friend operators

        friend extrapolation2DTable operator+ <Type>
        (
            const extrapolation2DTable&,
            const extrapolation2DTable&
        );

        friend extrapolation2DTable operator- <Type>
        (
            const extrapolation2DTable&,
            const extrapolation2DTable&
        );

        friend extrapolation2DTable operator* <Type>
	(
            const scalar,
            const extrapolation2DTable&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "extrapolation2DTable.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
