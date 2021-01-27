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

Original Author
    Frank Bos, TU Delft.  All rights reserved.
    Dubravko Matijasevic, FSB Zagreb.

Modified by 
    Xiaofeng Liu, UT San Antonio, Spring 2013

\*---------------------------------------------------------------------------*/

#include "ibmRBFInterpolation.H"
#include "demandDrivenData.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
Foam::scalarField Foam::ibmRBFInterpolation::calcPhiVector(const point& ip) const
{
  return RBF_->weights(controlPoints_, ip);
}


const Foam::scalarSquareMatrix& Foam::ibmRBFInterpolation::B() const
{
    if (!BPtr_)
    {
        calcB();
    }

    return *BPtr_;
}


void Foam::ibmRBFInterpolation::calcB() const
{
    // Determine inverse of connectivity matrix
    label polySize(4);

    if (!polynomials_)
    {
        polySize = 0;
    }

    // Fill Nb x Nb matrix
//    simpleMatrix<scalar> A(controlPoints_.size()+polySize);
    scalarSquareMatrix A(controlPoints_.size()+polySize);
    
    const label nControlPoints = controlPoints_.size();
    for (label i = 0; i < nControlPoints; i++)
    {
        scalarField weights = RBF_->weights(controlPoints_, controlPoints_[i]);

        for (label col = 0; col < nControlPoints; col++)
        {
            A[i][col] = weights[col];
        }
    }

    if (polynomials_)
    {
        for
        (
            label row = nControlPoints;
            row < nControlPoints + 1;
            row++
        )
        {
            for (label col = 0; col < nControlPoints; col++)
            {
                A[col][row] = 1.0;
                A[row][col] = 1.0;
            }
        }

        // Fill in X components of polynomial part of matrix
        for
        (
            label row = nControlPoints + 1;
            row < nControlPoints + 2;
            row++
        )
        {
            for (label col = 0; col < nControlPoints; col++)
            {
                A[col][row] = controlPoints_[col].x();
                A[row][col] = controlPoints_[col].x();
            }
        }

        // Fill in Y components of polynomial part of matrix
        for
        (
            label row = nControlPoints + 2;
            row < nControlPoints + 3;
            row++
        )
        {
            for (label col = 0; col < nControlPoints; col++)
            {
                A[col][row] = controlPoints_[col].y();
                A[row][col] = controlPoints_[col].y();
            }
        }
        // Fill in Z components of polynomial part of matrix
        for
        (
            label row = nControlPoints + 3;
            row < nControlPoints + 4;
            row++
        )
        {
            for (label col = 0; col < nControlPoints; col++)
            {
                A[col][row] = controlPoints_[col].z();
                A[row][col] = controlPoints_[col].z();
            }
        }

        // Fill 4x4 zero part of matrix
        for
        (
            label row = nControlPoints;
            row < nControlPoints + 4;
            row++
        )
        {
            for
            (
                label col = nControlPoints;
                col < nControlPoints + 4;
                col++
            )
            {
                A[row][col] = 0.0;
            }
        }
    }

    // HJ and FB (05 Jan 2009)
    // Collect ALL control points from ALL CPUs
    // Create an identical inverse for all CPUs

    Info<< "Inverting RBF motion matrix" << endl;

    BPtr_ = new scalarSquareMatrix(LUinvert(A));
}


void Foam::ibmRBFInterpolation::clearOut()
{
    deleteDemandDrivenData(BPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ibmRBFInterpolation::ibmRBFInterpolation
(
    const dictionary& dict,
    const vectorField& controlPoints,
    const vectorField& ipPoints
)
:
    dict_(dict),
    controlPoints_(controlPoints),
    ipPoints_(ipPoints),
    RBF_(RBFFunction::New(word(dict.lookup("RBF")), dict)),
    BPtr_(NULL),
    polynomials_(dict.lookup("polynomials"))
{}


Foam::ibmRBFInterpolation::ibmRBFInterpolation
(
    const ibmRBFInterpolation& rbf
)
:
    dict_(rbf.dict_),
    controlPoints_(rbf.controlPoints_),
    ipPoints_(rbf.ipPoints_),
    RBF_(rbf.RBF_->clone()),
    BPtr_(NULL),
    polynomials_(rbf.polynomials_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ibmRBFInterpolation::~ibmRBFInterpolation()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ibmRBFInterpolation::movePoints()
{
    clearOut();
}


// ************************************************************************* //
