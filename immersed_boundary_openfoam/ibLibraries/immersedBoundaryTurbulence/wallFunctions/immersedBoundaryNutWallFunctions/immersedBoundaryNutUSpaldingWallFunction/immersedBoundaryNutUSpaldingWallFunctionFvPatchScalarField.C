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

#include "immersedBoundaryNutUSpaldingWallFunctionFvPatchScalarField.H"
//#include "turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "immersedBoundaryVelocityWallFunctionFvPatchVectorField.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

tmp<scalarField> immersedBoundaryNutUSpaldingWallFunctionFvPatchScalarField::calcNut() const
{

    word UName_("U");
    word nuName_("nu");

    // Velocity
    const fvPatchVectorField& Uwg =
      patch().lookupPatchField<volVectorField, vector>(UName_);

    const immersedBoundaryVelocityWallFunctionFvPatchVectorField& Uw =
        refCast<const immersedBoundaryVelocityWallFunctionFvPatchVectorField>
        (
            Uwg
        );

    const scalarField& y = ibPatch().ibDelta();
    const vectorField& n = ibPatch().ibNormals();

    const scalarField magUp(mag((I - sqr(n)) & (Uw.ibSamplingPointValue()-Uw.ibValue())));

    //const scalarField magGradUw = mag(Uw.ibGrad());
    const scalarField magGradUw = magUp/y;
    // Laminar viscosity
    const fvPatchScalarField& nuwg = 
        patch().lookupPatchField<volScalarField, scalar>(nuName_);
 
    const immersedBoundaryFvPatchScalarField& nuw =
        refCast<const immersedBoundaryFvPatchScalarField>(nuwg);
    scalarField nu = nuw.ibValue();

    return max
    (
        scalar(0),
        sqr(calcUTau(magGradUw))/(magGradUw + ROOTVSMALL) - nu
    );
}


tmp<scalarField> immersedBoundaryNutUSpaldingWallFunctionFvPatchScalarField::calcUTau
(
    const scalarField& magGradU
) const
{
    word UName_("U");
    word nuName_("nu");

    // Distance from wall to IB point
    const scalarField& y = ibPatch().ibDelta();
    const vectorField& n = ibPatch().ibNormals();

    // Velocity
    const fvPatchVectorField& Uwg =
      patch().lookupPatchField<volVectorField, vector>(UName_);

    const immersedBoundaryVelocityWallFunctionFvPatchVectorField& Uw =
        refCast<const immersedBoundaryVelocityWallFunctionFvPatchVectorField>
        (
            Uwg
        );
    
    const scalarField magUp(mag((I - sqr(n)) & (Uw.ibSamplingPointValue()-Uw.ibValue())));

    // Laminar viscosity
    const fvPatchScalarField& nuwg = 
        patch().lookupPatchField<volScalarField, scalar>(nuName_);
 
    const immersedBoundaryFvPatchScalarField& nu =
        refCast<const immersedBoundaryFvPatchScalarField>(nuwg);
    scalarField nuw = nu.ibValue();


    const scalarField& nutw = this->wallValue();

    tmp<scalarField> tuTau(new scalarField(nuw.size(), 0.0));
    scalarField& uTau = tuTau();

    forAll(uTau, cellI)
    {
        scalar ut = sqrt((nutw[cellI] + nuw[cellI])*magGradU[cellI]);

        if (ut > ROOTVSMALL)
        {
            int iter = 0;
            scalar err = GREAT;

            do
            {

                scalar kUu = min(kappa_*magUp[cellI]/ut, 50);
                scalar fkUu = exp(kUu) - 1 - kUu*(1 + 0.5*kUu);
                scalar f =
                    - ut*y[cellI]/nuw[cellI]
                    + magUp[cellI]/ut
                    + 1/E_*(fkUu - 1.0/6.0*kUu*sqr(kUu));

                scalar df =
                    y[cellI]/nuw[cellI]
                  + magUp[cellI]/sqr(ut)
                  + 1/E_*kUu*fkUu/ut;

                scalar uTauNew = ut + f/df;
                err = mag((ut - uTauNew)/ut);
                ut = uTauNew;

            } while (ut > ROOTVSMALL && err > 0.01 && ++iter < 10);

            uTau[cellI] = max(0.0, ut);
        }
    }

    return tuTau;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

immersedBoundaryNutUSpaldingWallFunctionFvPatchScalarField::
immersedBoundaryNutUSpaldingWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    immersedBoundaryWallFunctionFvPatchScalarField(p, iF),
    kappa_(0.41),
    E_(9.8)

{}


immersedBoundaryNutUSpaldingWallFunctionFvPatchScalarField::
immersedBoundaryNutUSpaldingWallFunctionFvPatchScalarField
(
    const immersedBoundaryNutUSpaldingWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    immersedBoundaryWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    kappa_(ptf.kappa_),
    E_(ptf.E_)
{}


immersedBoundaryNutUSpaldingWallFunctionFvPatchScalarField::
immersedBoundaryNutUSpaldingWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    immersedBoundaryWallFunctionFvPatchScalarField(p, iF, dict),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    E_(dict.lookupOrDefault<scalar>("E", 9.8))
{}


immersedBoundaryNutUSpaldingWallFunctionFvPatchScalarField::
immersedBoundaryNutUSpaldingWallFunctionFvPatchScalarField
(
    const immersedBoundaryNutUSpaldingWallFunctionFvPatchScalarField& wfpsf
)
:
    immersedBoundaryWallFunctionFvPatchScalarField(wfpsf),
    kappa_(wfpsf.kappa_),
    E_(wfpsf.E_)
{}


immersedBoundaryNutUSpaldingWallFunctionFvPatchScalarField::
immersedBoundaryNutUSpaldingWallFunctionFvPatchScalarField
(
    const immersedBoundaryNutUSpaldingWallFunctionFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    immersedBoundaryWallFunctionFvPatchScalarField(wfpsf, iF),
    kappa_(wfpsf.kappa_),
    E_(wfpsf.E_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void immersedBoundaryNutUSpaldingWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    scalarField& nutw = this->wallValue();
    nutw=this->calcNut();
     Info<< "nutw " << min(nutw) << " " << max(nutw) << endl;
    this->wallMask() = true;
/*
    volScalarField& nut = const_cast<volScalarField&>
        (db().lookupObject<volScalarField>("nut"));

    const labelList& ibc = ibPatch().ibCells();

    const labelList& aibc = ibPatch().adjacentIbCells();

    const scalarField& aid = ibPatch().adjacentIbDelta();

    const scalarField& id = ibPatch().ibDelta();

    forAll (aibc, ibCellI)
    {

        nut[aibc[ibCellI]]=nut[aibc[ibCellI]]/(id[ibCellI]+SMALL)*aid[ibCellI]*0.01;
    }
*/
    immersedBoundaryWallFunctionFvPatchScalarField::updateCoeffs();
}


/*tmp<scalarField> immersedBoundaryNutUSpaldingWallFunctionFvPatchScalarField::yPlus() const
{
    const label patchi = patch().index();

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            dimensionedInternalField().group()
        )
    );
    const scalarField& y = turbModel.y()[patchi];
    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    return y*calcUTau(mag(Uw.snGrad()))/nuw;
}*/

void immersedBoundaryNutUSpaldingWallFunctionFvPatchScalarField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    // Insert omega values into the internal field
    this->setIbCellValues(this->wallValue());
Info <<"nut evaluate"<<endl;
    fvPatchScalarField::evaluate(commsType);
}
void immersedBoundaryNutUSpaldingWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    immersedBoundaryWallFunctionFvPatchScalarField::write(os);
    //writeLocalEntries(os);
    //writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    immersedBoundaryNutUSpaldingWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
