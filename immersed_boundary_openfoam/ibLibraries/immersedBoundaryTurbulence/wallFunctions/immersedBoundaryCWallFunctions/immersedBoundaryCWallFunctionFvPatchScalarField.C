/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "immersedBoundaryCWallFunctionFvPatchScalarField.H"
#include "immersedBoundaryVelocityWallFunctionFvPatchVectorField.H"
#include "turbulenceModel.H"

#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
 
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

immersedBoundaryCWallFunctionFvPatchScalarField::
immersedBoundaryCWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    immersedBoundaryWallFunctionFvPatchScalarField(p, iF),
    UName_("U"),
    kName_("k"),
    GName_("kEpsilon:G"),
    nuName_("nu"),
    nutName_("nut"),
    Cmu_(0.09),
    kappa_(0.41),
    E_(9.8)
{}


immersedBoundaryCWallFunctionFvPatchScalarField::
immersedBoundaryCWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    immersedBoundaryWallFunctionFvPatchScalarField(p, iF, dict),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    kName_(dict.lookupOrDefault<word>("k", "k")),
    GName_(dict.lookupOrDefault<word>("G", "kEpsilon:G")),
    nuName_(dict.lookupOrDefault<word>("nu", "nu")),
    nutName_(dict.lookupOrDefault<word>("nut", "nut")),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09)),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    E_(dict.lookupOrDefault<scalar>("E", 9.8))
{}


immersedBoundaryCWallFunctionFvPatchScalarField::
immersedBoundaryCWallFunctionFvPatchScalarField
(
    const immersedBoundaryCWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    immersedBoundaryWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    kName_(ptf.kName_),
    GName_(ptf.GName_),
    nuName_(ptf.nuName_),
    nutName_(ptf.nutName_),
    Cmu_(ptf.Cmu_),
    kappa_(ptf.kappa_),
    E_(ptf.E_)
{}


immersedBoundaryCWallFunctionFvPatchScalarField::
immersedBoundaryCWallFunctionFvPatchScalarField
(
    const immersedBoundaryCWallFunctionFvPatchScalarField& ewfpsf
)
:
    immersedBoundaryWallFunctionFvPatchScalarField(ewfpsf),
    UName_(ewfpsf.UName_),
    kName_(ewfpsf.kName_),
    GName_(ewfpsf.GName_),
    nuName_(ewfpsf.nuName_),
    nutName_(ewfpsf.nutName_),
    Cmu_(ewfpsf.Cmu_),
    kappa_(ewfpsf.kappa_),
    E_(ewfpsf.E_)
{}


immersedBoundaryCWallFunctionFvPatchScalarField::
immersedBoundaryCWallFunctionFvPatchScalarField
(
    const immersedBoundaryCWallFunctionFvPatchScalarField& ewfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    immersedBoundaryWallFunctionFvPatchScalarField(ewfpsf, iF),
    UName_(ewfpsf.UName_),
    kName_(ewfpsf.kName_),
    GName_(ewfpsf.GName_),
    nuName_(ewfpsf.nuName_),
    nutName_(ewfpsf.nutName_),
    Cmu_(ewfpsf.Cmu_),
    kappa_(ewfpsf.kappa_),
    E_(ewfpsf.E_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void immersedBoundaryCWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }


	if(!this->fixesValue())
	{
		// gradient at IB 
		this->ibGradRef()=this->wallGrad();
		this->setIbCellValues(this->wallValue());
	}
	else
	{
		this->ibValueRef()=this->wallValue();
	}
	
    this->wallMask() = true;

    // Insert epsilon values into the internal field
    immersedBoundaryWallFunctionFvPatchScalarField::updateCoeffs();
	//fvPatchScalarField::updateCoeffs();
}

void immersedBoundaryCWallFunctionFvPatchScalarField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    // Insert epsilon values into the internal field
    //this->setIbCellValues(this->wallValue());

    fvPatchScalarField::evaluate(commsType);
}


void immersedBoundaryCWallFunctionFvPatchScalarField::
write(Ostream& os) const
{
    immersedBoundaryWallFunctionFvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "k", "k", kName_);
    writeEntryIfDifferent<word>(os, "G", GName_, GName_);
    writeEntryIfDifferent<word>(os, "nu", "nu", nuName_);
    writeEntryIfDifferent<word>(os, "nut", "nut", nutName_);
    os.writeKeyword("Cmu") << Cmu_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappa") << kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("E") << E_ << token::END_STATEMENT << nl;
}

scalar immersedBoundaryCWallFunctionFvPatchScalarField::yPlusLam
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

tmp<scalarField> immersedBoundaryCWallFunctionFvPatchScalarField::calcUTau
(
    const scalarField& magGradU,
    const scalarField& y,
    const scalarField& magUp,
    const scalarField& nuw,
    const scalarField& nutw
) const
{
    //const label patchi = patch().index();

    //const scalarField& y = turbModel.y()[patchi];

    //const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    //const scalarField magUp(mag(Uw.patchInternalField() - Uw));

    //const tmp<scalarField> tnuw = turbModel.nu(patchi);
    //const scalarField& nuw = tnuw();

    //const scalarField& nutw = *this;

    tmp<scalarField> tuTau(new scalarField(y.size(), 0.0));
    scalarField& uTau = tuTau();

    forAll(uTau, faceI)
    {
        scalar ut = sqrt((nutw[faceI] + nuw[faceI])*magGradU[faceI]);

        if (ut > ROOTVSMALL)
        {
            int iter = 0;
            scalar err = GREAT;

            do
            {
                scalar kUu = min(kappa_*magUp[faceI]/ut, 50);
                scalar fkUu = exp(kUu) - 1 - kUu*(1 + 0.5*kUu);

                scalar f =
                    - ut*y[faceI]/nuw[faceI]
                    + magUp[faceI]/ut
                    + 1/E_*(fkUu - 1.0/6.0*kUu*sqr(kUu));

                scalar df =
                    y[faceI]/nuw[faceI]
                  + magUp[faceI]/sqr(ut)
                  + 1/E_*kUu*fkUu/ut;

                scalar uTauNew = ut + f/df;
                err = mag((ut - uTauNew)/ut);
                ut = uTauNew;

            } while (ut > ROOTVSMALL && err > 0.01 && ++iter < 10);

            uTau[faceI] = max(0.0, ut);
        }
    }

    return tuTau;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    immersedBoundaryCWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

 
} // End namespace Foam

// ************************************************************************* //
