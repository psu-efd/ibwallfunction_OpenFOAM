/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "LUinvert.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Do a LU inversion of a square matrix
Foam::SquareMatrix<Foam::scalar> Foam::LUinvert
(
    const scalarSquareMatrix& matrix
)
{
    scalarSquareMatrix luMatrix = matrix;

    scalarSquareMatrix luInvert(luMatrix.n());
    scalarField column(luMatrix.n());

    labelList pivotIndices(luMatrix.n());

    LUDecompose(luMatrix, pivotIndices);

    for (label j = 0; j < luMatrix.n(); j++)
    {
        for (label i = 0; i < luMatrix.n(); i++)
        {
            column[i] = 0.0;
        }

        column[j] = 1.0;

        LUBacksubstitute(luMatrix, pivotIndices, column);

        for (label i = 0; i < luMatrix.n(); i++)
        {
            luInvert[i][j] = column[i];
        }
    }

//    Pout << "matrix = " << matrix << endl;
//    Pout << "inverse = " << luInvert << endl;

    return luInvert;

}


// ************************************************************************* //
