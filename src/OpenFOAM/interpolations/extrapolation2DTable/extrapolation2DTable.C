/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Xavier Lamboley
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of my contributions of OpenFOAM.

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

#include "IFstream.H"
#include "openFoamTableReader.H"
#include "Vector.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class Type>
void Foam::extrapolation2DTable<Type>::readTable()
{
    fileName fName(fileName_);
    fName.expand();

    // Read data from file
    reader_()(fName, *this);

    if (this->empty())
    {
        FatalErrorIn
        (
            "Foam::extrapolation2DTable<Type>::readTable()"
        )   << "table read from " << fName << " is empty" << nl
            << exit(FatalError);
    }

    // Check that the data are in ascending order
    checkOrder();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::extrapolation2DTable<Type>::extrapolation2DTable()
:
    List<Tuple2<scalar, List<Tuple2<scalar, Type> > > >(),
    boundsHandling_(extrapolation2DTable::WARN),
    fileName_("fileNameIsUndefined"),
    reader_(NULL),
    isNull_(true)
{}


template<class Type>
Foam::extrapolation2DTable<Type>::extrapolation2DTable
(
    const List<Tuple2<scalar, List<Tuple2<scalar, Type> > > >& values,
    const boundsHandling bounds,
    const fileName& fName,
    const Switch isNull
)
:
    List<Tuple2<scalar, List<Tuple2<scalar, Type> > > >(values),
    boundsHandling_(bounds),
    fileName_(fName),
    reader_(NULL),
    isNull_(isNull)
{}


template<class Type>
Foam::extrapolation2DTable<Type>::extrapolation2DTable(const fileName& fName)
:
    List<Tuple2<scalar, List<Tuple2<scalar, Type> > > >(),
    boundsHandling_(extrapolation2DTable::WARN),
    fileName_(fName),
    reader_(new openFoamTableReader<Type>(dictionary())),
    isNull_(false)
{
    readTable();
}


template<class Type>
Foam::extrapolation2DTable<Type>::extrapolation2DTable(const dictionary& dict)
:
    List<Tuple2<scalar, List<Tuple2<scalar, Type> > > >(),
    boundsHandling_(wordToBoundsHandling(dict.lookup("outOfBounds"))),
    fileName_(dict.lookup("fileName")),
    reader_(tableReader<Type>::New(dict)),
    isNull_(false)
{
    readTable();
}


template<class Type>
Foam::extrapolation2DTable<Type>::extrapolation2DTable
(
     const extrapolation2DTable& extrapTable
)
:
    List<Tuple2<scalar, List<Tuple2<scalar, Type> > > >(extrapTable),
    boundsHandling_(extrapTable.boundsHandling_),
    fileName_(extrapTable.fileName_),
    reader_(extrapTable.reader_),    // note: steals reader. Used in write().
    isNull_(extrapTable.isNull_)
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::extrapolation2DTable<Type>::extrapolateValue
(
    const List<Tuple2<scalar, Type> >& data,
    const scalar lookupValue
) const
{
    label n = data.size();

    scalar minLimit = data.first().first();
    scalar maxLimit = data.last().first();
    if (n == 1)
      return data.first().second();
    if (lookupValue < minLimit)
    {
        switch (boundsHandling_)
        {
            case extrapolation2DTable::ERROR:
            {
                FatalErrorIn
                (
                    "Foam::extrapolation2DTable<Type>::extrapolateValue"
                    "("
                        "List<Tuple2<scalar, Type> >&, "
                        "const scalar"
                    ")"
                )   << "value (" << lookupValue << ") less than lower "
                    << "bound (" << minLimit << ")" << nl
                    << exit(FatalError);
                break;
            }
            case extrapolation2DTable::WARN:
            {
                WarningIn
                (
                    "Foam::extrapolation2DTable<Type>::extrapolateValue"
                    "("
                        "List<Tuple2<scalar, Type> >&, "
                        "const scalar"
                    ")"
                )   << "value (" << lookupValue << ") less than lower "
                    << "bound (" << minLimit << ")" << nl
                    << "    extrapolating the first entry"
                    << endl;
                // fall-through to 'EXTRAPOLATE'
            }
            case extrapolation2DTable::EXTRAPOLATE:
            {
		scalar x1 = minLimit;
		scalar x2 = data[1].first();
		Type y1 = data[0].second();
		Type y2 = data[1].second();
                return y1 + (lookupValue - x1)/(x2-x1) * (y2 - y1); //extrapolation
                break;
            }
        }
    }
    else if (lookupValue > maxLimit)
    {
        switch (boundsHandling_)
        {
            case extrapolation2DTable::ERROR:
            {
                FatalErrorIn
                (
                    "Foam::extrapolation2DTable<Type>::extrapolateValue"
                    "("
                        "List<Tuple2<scalar, Type> >&, "
                        "const scalar"
                    ")"
                )   << "value (" << lookupValue << ") greater than upper "
                    << "bound (" << maxLimit << ")" << nl
                    << exit(FatalError);
                break;
            }
            case extrapolation2DTable::WARN:
            {
                WarningIn
                (
                    "Foam::extrapolation2DTable<Type>::extrapolateValue"
                    "("
                        "List<Tuple2<scalar, Type> >&, "
                        "const scalar"
                    ")"
                )   << "value (" << lookupValue << ") greater than upper "
                    << "bound (" << maxLimit << ")" << nl
                    << "    Continuing with the last entry"
                    << endl;
                // fall-through to 'EXTRAPOLATE'
            }
            case extrapolation2DTable::EXTRAPOLATE:
            {
		 scalar x1 = maxLimit;
		 scalar x2 = data[data.size()-2].first();
		 Type y1 = data[data.size()-1].second();
		 Type y2 = data[data.size()-2].second();
                return y1 + (lookupValue - x1)/(x2-x1) * (y2 - y1); //extrapolation
                return data.last().second();
                break;
            }
        }
    }

    // look for the correct range in X
    label lo = 0;
    label hi = 0;

    for (label i = 0; i < n; ++i)
    {
        if (lookupValue >= data[i].first())
        {
            lo = hi = i;
        }
        else
        {
            hi = i;
            break;
        }
    }

    if (lo == hi)
    {
        return data[lo].second();
    }
    else
    {
        Type m =
            (data[hi].second() - data[lo].second())
           /(data[hi].first() - data[lo].first());

        // normal interpolation
        return data[lo].second() + m*(lookupValue - data[lo].first());
    }
}


template<class Type>
template<class L>
Foam::label Foam::extrapolation2DTable<Type>::Xi
(  
    const List<Tuple2<scalar, L> >& t,
    const scalar valueX,
    const bool reverse
) const
{
    /*
    label limitI = 0;
    if (reverse)
    {
        limitI = t.size() - 1;
    }

    if (bop(valueX, t[limitI].first()))
    {
        switch (boundsHandling_)
        {
            case extrapolation2DTable::ERROR:
            {
                FatalErrorIn
                (
                    "Foam::label Foam::extrapolation2DTable<Type>::Xi"
                    "("
                        "const BinaryOp&, "
                        "const scalar, "
                        "const bool"
                    ") const"
                )   << "value (" << valueX << ") out of bounds"
                    << exit(FatalError);
                break;
            }
            case extrapolation2DTable::WARN:
            {
                WarningIn
                (
                    "Foam::label Foam::extrapolation2DTable<Type>::Xi"
                    "("
                        "const BinaryOp&, "
                        "const scalar, "
                        "const bool"
                )   << "value (" << valueX << ") out of bounds"
                    << endl;
                // fall-through to 'EXTRAPOLATE'
            }
            case extrapolation2DTable::EXTRAPOLATE:
            {
                return limitI;
            }
            default:
            {
                FatalErrorIn
                (
                    "Foam::label Foam::extrapolation2DTable<Type>::Xi"
                    "("
                        "const BinaryOp&, "
                        "const scalar, "
                        "const bool"
                    ") const"
                )
                    << "Un-handled enumeration " << boundsHandling_
                    << abort(FatalError);
            }
        }
    }
*/
    label nX = t.size();
    if (reverse)
    {
	for(int i = nX - 1; i > 0; i--)
	  if(t[i].first() <= valueX)
	  {
	    if(i == nX - 1) //point lies outwards, return last index
	      return nX - 1;
	    else
	      return i + 1;
	  }
	return 1; //point lies outwards, return 2nd index
    }
    else
    {
	for(int i = 0; i < nX - 1; i++)
	  if(t[i].first() > valueX)
	  {
	    if(i == 0)
	      return 0; //point lies outwards, return first index
	    else
	      return i-1;
	  }
	return nX - 2; //point lies outwards, return second last index
    }
}

// I'm sure it's not very OpenFOAM/C++ style. It has to be improved. But it works.
template<class Type>
Type Foam::extrapolation2DTable<Type>::Tderivative
(
    const scalar valueX,
    const scalar valueY
) const
{
    const table& t = *this;

    if (t.size() <= 1)
    {
	WarningIn
        (
            "Type Foam::extrapolation2DTable<Type>::Tderivative"
            "("
	        "const scalar, "
  	        "const scalar"
            ") const"
        )
            << "cannot derivate a zero- or one-sized table - returning zero" << endl;
	return pTraits<Type>::zero;
    }

    // have 2-D data, extrapolate
    List<Tuple2<scalar, Type> > row0;
    List<Tuple2<scalar, Type> > row1;
    // find low and high indices in the X range that bound valueX
    label x0i = Xi(t, valueX, false);
    label x1i = Xi(t, valueX, true);
    row0 = t[x0i].second();
    row1 = t[x1i].second();
    //factor for interpolating between both rows
    scalar factor = (valueX - t[x0i].first()) / (t[x1i].first() - t[x0i].first());
    if(t.first().second().size() == 1)
	return (row1.first().second() - row0.first().second())/(t[x1i].first() - t[x0i].first());

    // find low and high indices in the Y range that bound valueY
    label y00 = Xi(row0, valueY, false);
    label y01 = Xi(row0, valueY, true);
    label y10 = Xi(row1, valueY, false);
    label y11 = Xi(row1, valueY, true);
	
    label ymin = min(y00, y10);
    label ymax = max(y01, y11);
    List<vector> points;
    label index = 0;
    for(label i = ymin; i <= ymax; i++)
    {
	points.append(vector(
			  row0[i].second(),
			  row0[i].first()  + factor * (row1[i].first()  - row0[i].first() ),
			  row1[i].second()
			  ));
	if(points[i - ymin].y() < valueY)
	    index = i - ymin;
    }
    return (points[index + 1].z() - points[index].x()) / (t[x1i].first() - t[x0i].first());
}
    

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
Type Foam::extrapolation2DTable<Type>::operator()
(
    const scalar valueX,
    const scalar valueY
) const
{
    
    // Considers all of the list in Y being equal
    label nX = this->size();
    label nY = this->first().second().size();
    //if(nY == 1)
      //Info << "operator(" << valueX << ", " << valueY << ")" << endl;
    //Info << "size: " << nX << " / " << nY << endl;
    const table& t = *this;
    if (nX == 0)
    {
        WarningIn
        (
            "Type Foam::extrapolation2DTable<Type>::operator()"
            "("
                "const scalar, "
                "const scalar"
            ") const"
        )
            << "cannot extrapolate a zero-sized table - returning zero" << endl;

        return pTraits<Type>::zero;
    }
    else if (nX == 1)
    {
        // only 1 row (in X) - extrapolate to find Y value
        return extrapolateValue(t.first().second(), valueY);
    }
    else
    {
        // have 2-D data, extrapolate
	//double minRowLimit = t.first().first();
	//double maxRowLimit = t.last().first();
	//idea: use plane-based inter/extrapolation, find the 4 points which surround the wanted point first
	/*vector p00;
	vector p01;
	vector p10;
	vector p11;
	*/
	List<Tuple2<scalar, Type> > row0;
	List<Tuple2<scalar, Type> > row1;
	// find low and high indices in the X range that bound valueX
	label x0i = Xi(t, valueX, false);
	label x1i = Xi(t, valueX, true);
	row0 = t[x0i].second();
	row1 = t[x1i].second();
	//factor for interpolating between both rows
	scalar factor = (valueX - t[x0i].first()) / (t[x1i].first() - t[x0i].first());
	if(nY == 1)
	  return row0.first().second() + factor * (row1.first().second() - row0.first().second());
	/*}
	else if(valueX < minRowLimit && boundsHandling_ != extrapolation2DTable::ERROR)
	{
	  if(boundsHandling_ == extrapolation2DTable::WARN)
	  {
		WarningIn
		(
			"Foam::label Foam::extrapolation2DTable<Type>::operator()"
			"("
			"const scalar valueX,"
			"const scalar valueY"
			") const"
		)   << "value (" << valueX << ") out of bounds"
		<< endl;
	  }
	  
	  p00.x() = t[0].first();
	  p01.x() = t[0].first();
	  p10.x() = t[1].first();
	  p11.x() = t[1].first();
	  row0 = t[0].second();
	  row1 = t[1].second();
	}
	else if(valueX > maxRowLimit && boundsHandling_ != extrapolation2DTable::ERROR)
	{
	  if(boundsHandling_ == extrapolation2DTable::WARN)
	  {
		WarningIn
		(
			"Foam::label Foam::extrapolation2DTable<Type>::operator()"
			"("
			"const scalar valueX,"
			"const scalar valueY"
			") const"
		)   << "value (" << valueX << ") out of bounds"
		<< endl;
	  }
	  
	  p00.x() = t[t.size()-2].first();
	  p01.x() = t[t.size()-2].first();
	  p10.x() = t[t.size()-1].first();
	  p11.x() = t[t.size()-1].first();
	  row0 = t[t.size()-2].second();
	  row1 = t[t.size()-1].second();
	}
	else //extrapolation needed but ERROR boundshandling
	{
	  FatalErrorIn
	  (
	  "Foam::label Foam::extrapolation2DTable<Type>::operator()"
	  "("
		  "const scalar valueX,"
		  "const scalar valueY"
	  ") const"
	  )   << "value (" << valueX << ") out of bounds"
	  << exit(FatalError);
	  return 0;
	}
	*/
	//Now we have the data for the first axis.
	
	//Info << "first axis" << endl;
	// find low and high indices in the Y range that bound valueY
	label y00 = Xi(row0, valueY, false);
	label y01 = Xi(row0, valueY, true);
	label y10 = Xi(row1, valueY, false);
	label y11 = Xi(row1, valueY, true);
	
	//Info << "y00:" << y00 << endl;
	//Info << "y01:" << y01 << endl;
	//Info << "y10:" << y10 << endl;
	//Info << "y11:" << y11 << endl;
	/*
	//if valueY is above the limit of the first axis
	if(y00 == y01 && row0[y00].first() < valueY)
	  y00--;
	//if valueY is below the limit of the first axis
	else if(y00 == y01 && row0[y01].first() > valueY)
	  y01++;
	
	//if valueY is above the limit of the first axis
	if(y10 == y11 && row1[y10].first() < valueY)
	  y10--;
	//if valueY is below the limit of the first axis
	else if(y10 == y11 && row1[y11].first() > valueY)
	  y11++;
	*/
	label ymin = min(y00, y10);
	label ymax = max(y01, y11);
	//Info << "ymin: " << ymin << ", ymax: " << ymax << endl;
	List<vector> points;
	//Info << "x values of rows: " << t[x0i].first() << ", " << t[x1i].first() << ", factor x: " << factor << endl;
	label index = 0;
	for(label i = ymin; i <= ymax; i++)
	{
	  points.append(vector(
	    valueX,
	    row0[i].first()  + factor * (row1[i].first()  - row0[i].first() ),
	    row0[i].second() + factor * (row1[i].second() - row0[i].second())
	  ));
	  //Info << "point " << i << ": " << points[i - ymin] << endl;
	  if(points[i - ymin].y() < valueY)
	    index = i - ymin;
	}
	//Info << "size: " << points.size() << ", index: " << index << ", index in original list: " << index + ymin << endl;
	return (points[index] + (valueY - points[index].y()) / (points[index + 1].y() - points[index].y()) * (points[index + 1] - points[index])).z();
	/*
	//Now we have the indices of 4 points in the table that lie around the searched point
	//Next, get the y and z coordinate of these points:
	p00.y() = row0[y00].first();
	p01.y() = row0[y01].first();
	p00.z() = row0[y00].second();
	p01.z() = row0[y01].second();
	
	p10.y() = row1[y10].first();
	p11.y() = row1[y11].first();
	p10.z() = row1[y10].second();
	p11.z() = row1[y11].second();
	
	vector p0i = p00 + (valueX - p00.x()) / (p10.x() - p00.x()) * (p10 - p00);
	vector p1i = p01 + (valueX - p01.x()) / (p11.x() - p01.x()) * (p11 - p01);
	vector final = p0i + (valueY - p0i.y()) / (p1i.y() - p0i.y()) * (p1i - p0i);
	*/
	/*
	vector pi0 = p00+(valueY-p00.y())/(p01.y()-p00.y()) * (p01 - p00);
	vector pi1 = p10 + (valueY - p10.y()) / (p11.y() - p10.y()) * (p11 - p10);
	vector final = pi0 + (valueX - pi0.x()) / (pi1.x() - pi0.x()) * (pi1 - pi0);
	*/
	
	
	//return final.z();
	/*Info << "p00:" << p00 << endl;
	Info << "p01:" << p01 << endl;
	Info << "p10:" << p10 << endl;
	Info << "p11:" << p11 << endl;*/
	//Now that we know these points, we can construct all triangles between them and figure out in
	//which of the triangles the point lies. Then we can interpolate the point on the plane
	//spanned up from those the triangle.
	
	/*
	//The wanted point lies in one of the triangles spanned up from the points in these lists
	//Unfortunately, these lists can have an arbitrary size for ill-conditioned sampling points
	List<vector> lPoints0;
	List<vector> lPoints1;
	label ymin = y00 < y10 ? y00 : y10;
	label ymax = y01 > y11 ? y01 : y11;
	
	for(int i = ymin; i <= ymax; i++)
	{
	  lPoints0.append(vector(p00.x(), row0[i].first(), row0[i].second()));
	  lPoints1.append(vector(p10.x(), row1[i].first(), row1[i].second()));
	}
	
	//Info << "lists" << endl;
	//Strategy: construct a triangle, check if the point lies in it, if yes, return the interpolated value.
	//How to extrapolate? If we are extrapolating, we already know the triangle and don't need to do bounds checking!
	//Triangle structure:
	// |/|
	vector p1 = lPoints0[0];
	vector p2 = lPoints1[0];
	vector p3 = lPoints1[1];
	
	//Info << "P1:" << p1 << endl;
	//Info << "P2:" << p2 << endl;
	//Info << "P3:" << p3 << endl;
	int i = 0, j = 0;
	while(true)
	{
	  //For interpolation we use barycentric coordinates, see: 
	  //http://en.wikipedia.org/wiki/Barycentric_coordinate_system
	  //They can also be used for checking if a point lies inside a triangle.
	  scalar detT = ((p2.y()-p3.y())*(p1.x()-p3.x())+(p3.x()-p2.x())*(p1.y()-p3.y()));
	  
	  //Info << "detT: " << detT << endl;
	  scalar lambda1 = ((p2.y()-p3.y())*(valueX-p3.x())+(p3.x()-p2.x())*(valueY-p3.y())) / detT;
	  scalar lambda2 = ((p3.y()-p1.y())*(valueX-p3.x())+(p1.x()-p3.x())*(valueY-p3.y())) / detT;
	  scalar lambda3 = 1 - lambda1 -lambda2;
	  //Info << "lambda1: " << lambda1 << endl;
	  //Info << "lambda2: " << lambda2 << endl;
	  //Info << "lambda3: " << lambda3 << endl;
	  scalar epsilon = 1e-10;
	  if (lambda1 >= -epsilon && lambda1 <= 1 + epsilon && lambda2 >= -epsilon && lambda2 <= 1 + epsilon && lambda3 >= -epsilon && lambda3 <= 1 + epsilon)
	    return lambda1 * p1.z() + lambda2 * p2.z() + lambda3 * p3.z();
	  
	  //if the above condition is not fulfilled, the point lies outside of the triangle and we need to check the next one
	  if((i + j) % 2 == 0)
	  {
	    //Info << "even step" << endl;
	    j++;
	    if(lPoints0.size() > i+1)
	      p2 = lPoints0[i+1];
	    else
	    {
	      FatalErrorIn
	      (
		"Foam::label Foam::extrapolation2DTable<Type>::operator()"
		"("
			"const scalar valueX,"
			"const scalar valueY"
		") const"
	      ) << "Couldn't find matching triangle for " << valueX << "/" << valueY << ", filename: " << this->fileName_ << exit(FatalError);
	      return 0;
	    }
	  }
	  else
	  {
	    //Info << "odd step" << endl;
	    i++;
	    p1 = p2;
	    p2 = p3;
	    if(lPoints1.size() > j+1)
	      p3 = lPoints1[j+1];
	    else
	    {
	      FatalErrorIn
	      (
		"Foam::label Foam::extrapolation2DTable<Type>::operator()"
		"("
			"const scalar valueX,"
			"const scalar valueY"
		") const"
	      ) << "Couldn't find matching triangle for " << valueX << "/" << valueY << ", filename: " << this->fileName_ << exit(FatalError);
	      return 0;
	    }
	  }
	}
    */
    }
}


template<class Type>
Foam::word Foam::extrapolation2DTable<Type>::boundsHandlingToWord
(
     const boundsHandling& bound
) const
{
    word enumName("warn");

    switch (bound)
    {
        case extrapolation2DTable::ERROR:
        {
            enumName = "error";
            break;
        }
        case extrapolation2DTable::WARN:
        {
            enumName = "warn";
            break;
        }
        case extrapolation2DTable::EXTRAPOLATE:
        {
            enumName = "extrapolate";
            break;
        }
    }

    return enumName;
}


template<class Type>
typename Foam::extrapolation2DTable<Type>::boundsHandling
Foam::extrapolation2DTable<Type>::wordToBoundsHandling
(
    const word& bound
) const
{
    if (bound == "error")
    {
        return extrapolation2DTable::ERROR;
    }
    else if (bound == "warn")
    {
        return extrapolation2DTable::WARN;
    }
    else if (bound == "extrapolate")
    {
        return extrapolation2DTable::EXTRAPOLATE;
    }
    else
    {
        WarningIn
        (
            "Foam::extrapolation2DTable<Type>::wordToBoundsHandling"
            "("
            "    const word&"
            ")"
        )   << "bad outOfBounds specifier " << bound << " using 'warn'" << endl;

        return extrapolation2DTable::WARN;
    }
}

template<class Type>
inline Foam::extrapolation2DTable<Type>& Foam::extrapolation2DTable<Type>::operator=
(
    const extrapolation2DTable<Type>& et
)
{
    boundsHandling_ = et.boundsHandling_;
    fileName_ = et.fileName_;
    reader_ = et.reader_;
    isNull_ = et.isNull_;

    return *this;
}

template<class Type>
inline Foam::extrapolation2DTable<Type> Foam::operator+
(
    const extrapolation2DTable<Type>& et1,
    const extrapolation2DTable<Type>& et2
)
{
    List<Tuple2<scalar, List<Tuple2<scalar, Type> > > > etn = et1, ett = et2;

    if (et1.size() != et2.size())
    {
	FatalErrorInFunction
            << "attempt to sum list with different sizes."
            << abort(FatalError);
    }

    if (et1.isNull_)
    {
    	return extrapolation2DTable<Type>
        (
            ett,
    	    et2.boundsHandling_,	
    	    et2.fileName_,
    	    et2.isNull_
        );
    }
    else if (et2.isNull_)
    {
    	return extrapolation2DTable<Type>
        (
            etn,
    	    et1.boundsHandling_,	
    	    et1.fileName_,
    	    et1.isNull_
        );
    }
    
    for (int i = 0 ; i < etn.size() ; i ++)
    {
	for (int j = 0 ; j < etn[i].second().size() ; j ++)
	{
	    etn[i].second()[j].second() += ett[i].second()[j].second();
	}
    }	    
    
    return extrapolation2DTable<Type>
    (
	etn,
	et1.boundsHandling_,	
	et1.fileName_,
	false
    );
}

template<class Type>
inline Foam::extrapolation2DTable<Type> Foam::operator-
(
    const extrapolation2DTable<Type>& et1,
    const extrapolation2DTable<Type>& et2
)
{
    List<Tuple2<scalar, List<Tuple2<scalar, Type> > > > etn = et1, ett = et2;

    if (et1.size() != et2.size())
    {
	FatalErrorInFunction
            << "attempt to sum list with different sizes."
            << abort(FatalError);
    }
    
    for (int i = 0 ; i < etn.size() ; i ++)
    {
	for (int j = 0 ; j < etn[i].second().size() ; j ++)
	{
	    etn[i].second()[j].second() -= ett[i].second()[j].second();
	}
    }	    
    
        return extrapolation2DTable<Type>
        (
            etn,
	    et1.boundsHandling_,	
	    et1.fileName_,
	    et1.isNull_
        );
}

template<class Type>
inline Foam::extrapolation2DTable<Type> Foam::operator*
(
    const scalar s,
    const extrapolation2DTable<Type>& et
)
{
    List<Tuple2<scalar, List<Tuple2<scalar, Type> > > > etn = et;

    // I'm not sure that it saves time but I try it.
    if (s == 1 || et.isNull_)
    {
    	return extrapolation2DTable<Type>
        (
    	    etn,
    	    et.boundsHandling_,	
    	    et.fileName_,
	    et.isNull_
        );
    }
    else
    {
    	for (int i = 0 ; i < etn.size() ; i ++)
    	{
    	    for (int j = 0 ; j < etn[i].second().size() ; j ++)
    	    {
    		etn[i].second()[j].second() *= s;
    	    }
    	}
    }

    if (s == 0)
    {
    	return extrapolation2DTable<Type>
    	(
    	    etn,
    	    et.boundsHandling_,	
    	    et.fileName_,
    	    true
    	);
    }
    else
    {
	return extrapolation2DTable<Type>
	(
	    etn,
	    et.boundsHandling_,	
	    et.fileName_,
	    et.isNull_
	);
    }
}

template<class Type>
typename Foam::extrapolation2DTable<Type>::boundsHandling
Foam::extrapolation2DTable<Type>::outOfBounds
(
    const boundsHandling& bound
)
{
    boundsHandling prev = boundsHandling_;
    boundsHandling_ = bound;
    return prev;
}


template<class Type>
void Foam::extrapolation2DTable<Type>::checkOrder() const
{
    label n = this->size();
    const table& t = *this;

    scalar prevValue = t[0].first();

    for (label i=1; i<n; ++i)
    {
        const scalar currValue = t[i].first();

        // avoid duplicate values (divide-by-zero error)
        if (currValue <= prevValue)
        {
            FatalErrorIn
            (
                "Foam::extrapolation2DTable<Type>::checkOrder() const"
            )   << "out-of-order value: "
                << currValue << " at index " << i << nl
                << exit(FatalError);
        }
        prevValue = currValue;
    }
}


template<class Type>
void Foam::extrapolation2DTable<Type>::write(Ostream& os) const
{
    os.writeKeyword("fileName")
        << fileName_ << token::END_STATEMENT << nl;
    os.writeKeyword("outOfBounds")
        << boundsHandlingToWord(boundsHandling_) << token::END_STATEMENT << nl;

    *this >> os;
}


// ************************************************************************* //
