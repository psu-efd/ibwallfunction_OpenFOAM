/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    ibmRBFInterpolation

Description
    Radial basis function interpolation class
    Modified based on the original version in OF for immersed boundary method.

Description
    Interpolation class which uses Radial Basis Functions to interpolate the
    fluid property (velocity, pressure, etc.) to an image point 
    based on given cloud (neighbor) values.

    The coefficient vectors, alpha and beta are determined by solving
    the system:

	| db | = | Mbb Pb | | alpha |
	| 0  |   | Pb  0  | |  beta |

	where db are the given cloud values,
	Mbb the RBF correlation matrix (NbxNb), containing RBF evaluations
    at the cloud centers, and Pb some linear polynomial matrix (Nbx4).

	Those coefficients are calculated every timestep, with the current
    cloud values db, with the inverse of Mbb. Using those
    coefficients, the RBF is evaluated at image points every
    timestep.

Original Author
    Frank Bos, TU Delft.  All rights reserved.
    Dubravko Matijasevic, FSB Zagreb.

Modification by 
    Xiaofeng Liu, UT San Antonio, Spring 2013

SourceFiles
    ibmRBFInterpolation.C
    ibmRBFInterpolationTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef ibmRBFInterpolation_H
#define ibmRBFInterpolation_H

#include "dictionary.H"
#include "RBFFunction.H"
#include "point.H"
#include "Switch.H"

#include "LUinvert.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class ibmRBFInterpolation Declaration
\*---------------------------------------------------------------------------*/

class ibmRBFInterpolation
{
    // Private data

        //- Dictionary
        const dictionary& dict_;

        //- Reference to control points (cloud of neighbors)
        const vectorField& controlPoints_;

        //- Rerefence to image points (could be more than one)
        const vectorField& ipPoints_;

        //- RBF function
        autoPtr<RBFFunction> RBF_;

        //- Interpolation matrix
        mutable scalarSquareMatrix* BPtr_;

        //- Add polynomials to RBF matrix
        Switch polynomials_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const ibmRBFInterpolation&);

        //- Calculate interpolation matrix
        void calcB() const;

        //- Clear out
        void clearOut();


public:

    // Constructors

        //- Construct from components
        ibmRBFInterpolation
        (
            const dictionary& dict,
            const vectorField& controlPoints,
            const vectorField& ipPoints
        );

        //- Construct as copy
        ibmRBFInterpolation(const ibmRBFInterpolation&);


    // Destructor

        ~ibmRBFInterpolation();


    // Member Functions
        //- Return interpolation matrix
        const scalarSquareMatrix& B() const;

        //- Interpolate
        template<class Type>
        tmp<Field<Type> > interpolate(const Field<Type>& ctrlField) const;

        //- Move points
        void movePoints();

        //- Calculate the Phi vector given the image point
        scalarField calcPhiVector(const point& ip) const;
       
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ibmRBFInterpolationTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
