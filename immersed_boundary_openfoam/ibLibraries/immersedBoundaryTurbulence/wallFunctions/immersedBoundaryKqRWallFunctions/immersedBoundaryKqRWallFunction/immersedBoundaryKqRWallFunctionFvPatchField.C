/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "immersedBoundaryKqRWallFunctionFvPatchField.H"

#include "fvPatchFieldMapper.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void immersedBoundaryKqRWallFunctionFvPatchField<Type>::checkType()
{
    if (!isA<immersedBoundaryFvPatch>(this->patch()))
    {
        FatalErrorIn("immersedBoundaryKqRWallFunctionFvPatchField::checkType()")
            << "Invalid wall function specification" << nl
            << "    Patch type for patch " << this->patch().name()
            << " must be immersedBoundary" << nl
            << "    Current patch type is " << this->patch().type()
            << nl << endl
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
immersedBoundaryKqRWallFunctionFvPatchField<Type>::immersedBoundaryKqRWallFunctionFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    immersedBoundaryWallFunctionFvPatchField<Type>(p, iF)
{
    checkType();
}


template<class Type>
immersedBoundaryKqRWallFunctionFvPatchField<Type>::immersedBoundaryKqRWallFunctionFvPatchField
(
    const immersedBoundaryKqRWallFunctionFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    immersedBoundaryWallFunctionFvPatchField<Type>(ptf, p, iF, mapper)
{
    checkType();
}


template<class Type>
immersedBoundaryKqRWallFunctionFvPatchField<Type>::immersedBoundaryKqRWallFunctionFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    immersedBoundaryWallFunctionFvPatchField<Type>(p, iF, dict)
{
    checkType();
}


template<class Type>
immersedBoundaryKqRWallFunctionFvPatchField<Type>::immersedBoundaryKqRWallFunctionFvPatchField
(
    const immersedBoundaryKqRWallFunctionFvPatchField& tkqrwfpf
)
:
    immersedBoundaryWallFunctionFvPatchField<Type>(tkqrwfpf)
{
    checkType();
}


template<class Type>
immersedBoundaryKqRWallFunctionFvPatchField<Type>::immersedBoundaryKqRWallFunctionFvPatchField
(
    const immersedBoundaryKqRWallFunctionFvPatchField& tkqrwfpf,
    const DimensionedField<Type, volMesh>& iF
)
:
    immersedBoundaryWallFunctionFvPatchField<Type>(tkqrwfpf, iF)
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class Type>
void immersedBoundaryKqRWallFunctionFvPatchField<Type>::updateCoeffs()
{
    //Pout<<"void immersedBoundaryKqRWallFunctionFvPatchField<Type>::updateCoeffs()"<<endl;
    immersedBoundaryWallFunctionFvPatchField<Type>::updateCoeffs();
    //this->setIbCellValues(this->ibCellValue());//zeroGradient
    fvPatchField<Type>::updateCoeffs();
    // Insert nut values into the internal field


}

template<class Type>
void immersedBoundaryKqRWallFunctionFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes commsType
)
{



   // Pout<<"immersedBoundaryKqRWallFunctionFvPatchField<Type>::evaluate"<<endl;
    immersedBoundaryWallFunctionFvPatchField<Type>::evaluate();
 
}


template<class Type>
void immersedBoundaryKqRWallFunctionFvPatchField<Type>::write(Ostream& os) const
{
    immersedBoundaryWallFunctionFvPatchField<Type>::write(os);
    //this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
