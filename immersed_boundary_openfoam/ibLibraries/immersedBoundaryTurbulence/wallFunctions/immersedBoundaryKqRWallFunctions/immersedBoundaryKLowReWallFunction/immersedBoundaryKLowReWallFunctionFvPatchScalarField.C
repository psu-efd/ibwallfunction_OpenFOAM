/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2015 OpenFOAM Foundation
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

#include "immersedBoundaryKLowReWallFunctionFvPatchScalarField.H"
#include "immersedBoundaryKqRWallFunctionFvPatchFields.H"
//#include "turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void immersedBoundaryKLowReWallFunctionFvPatchScalarField::checkType()
{
    if (!isA<immersedBoundaryFvPatch>(patch()))
    {
        FatalErrorIn("immersedBoundaryKLowReWallFunctionFvPatchScalarField::checkType()")
            << "Invalid immersedBoundary wall function specification" << nl
            << "    Patch type for patch " << patch().name()
            << " must be immersedBoundary" << nl
            << "    Current patch type is " << patch().type() << nl << endl
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

scalar immersedBoundaryKLowReWallFunctionFvPatchScalarField::yPlusLam
(
    const scalar kappa,
    const scalar E
)
{
    scalar ypl = 11.0;

    for (int i=0; i<10; i++)
    {
        ypl = log(max(E*ypl, 1))/kappa;
    }

    return ypl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

immersedBoundaryKLowReWallFunctionFvPatchScalarField::immersedBoundaryKLowReWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    immersedBoundaryWallFunctionFvPatchField<scalar>(p, iF),
    Cmu_(0.09),
    kappa_(0.41),
    E_(9.8),
    Ceps2_(1.9),
    yPlusLam_(yPlusLam(kappa_, E_))
{
    checkType();
}


immersedBoundaryKLowReWallFunctionFvPatchScalarField::immersedBoundaryKLowReWallFunctionFvPatchScalarField
(
    const immersedBoundaryKLowReWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    immersedBoundaryWallFunctionFvPatchField<scalar>(ptf, p, iF, mapper),
    Cmu_(ptf.Cmu_),
    kappa_(ptf.kappa_),
    E_(ptf.E_),
    Ceps2_(ptf.Ceps2_),
    yPlusLam_(ptf.yPlusLam_)
{
    checkType();
}


immersedBoundaryKLowReWallFunctionFvPatchScalarField::immersedBoundaryKLowReWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    immersedBoundaryWallFunctionFvPatchField<scalar>(p, iF, dict),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09)),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    E_(dict.lookupOrDefault<scalar>("E", 9.8)),
    Ceps2_(dict.lookupOrDefault<scalar>("Ceps2", 1.9)),
    yPlusLam_(yPlusLam(kappa_, E_))
{
    checkType();
}


immersedBoundaryKLowReWallFunctionFvPatchScalarField::immersedBoundaryKLowReWallFunctionFvPatchScalarField
(
    const immersedBoundaryKLowReWallFunctionFvPatchScalarField& kwfpsf
)
:
    immersedBoundaryWallFunctionFvPatchField<scalar>(kwfpsf),
    Cmu_(kwfpsf.Cmu_),
    kappa_(kwfpsf.kappa_),
    E_(kwfpsf.E_),
    Ceps2_(kwfpsf.Ceps2_),
    yPlusLam_(kwfpsf.yPlusLam_)
{
    checkType();
}


immersedBoundaryKLowReWallFunctionFvPatchScalarField::immersedBoundaryKLowReWallFunctionFvPatchScalarField
(
    const immersedBoundaryKLowReWallFunctionFvPatchScalarField& kwfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    immersedBoundaryWallFunctionFvPatchField<scalar>(kwfpsf, iF),
    Cmu_(kwfpsf.Cmu_),
    kappa_(kwfpsf.kappa_),
    E_(kwfpsf.E_),
    Ceps2_(kwfpsf.Ceps2_),
    yPlusLam_(kwfpsf.yPlusLam_)
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void immersedBoundaryKLowReWallFunctionFvPatchScalarField::updateCoeffs()
{
    immersedBoundaryWallFunctionFvPatchScalarField::updateCoeffs();
 
}


void immersedBoundaryKLowReWallFunctionFvPatchScalarField::evaluate
(
    const Pstream::commsTypes commsType
)
{

 
 
    fvPatchScalarField::evaluate(commsType);
}


void immersedBoundaryKLowReWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    immersedBoundaryFvPatchScalarField::write(os);
    os.writeKeyword("Cmu") << Cmu_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappa") << kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("E") << E_ << token::END_STATEMENT << nl;
    os.writeKeyword("Ceps2") << Ceps2_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    immersedBoundaryKLowReWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
