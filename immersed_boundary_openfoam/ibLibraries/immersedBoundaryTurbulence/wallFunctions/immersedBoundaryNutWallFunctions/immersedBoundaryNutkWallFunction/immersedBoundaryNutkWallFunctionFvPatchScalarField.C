/*---------------------------------------------------------------------------* \
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

#include "immersedBoundaryNutkWallFunctionFvPatchScalarField.H"
//#include "turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField> immersedBoundaryNutkWallFunctionFvPatchScalarField::calcNut() const
{
 

    word UName_("U");
    word nuName_("nu");
    word kName_("k");

    // Distance from wall to IB point
    const scalarField& y = ibPatch().ibSamplingPointDelta();

    const scalar yPlusLam_ = yPlusLam(kappa_, E_);

    // Turbulence kinetic energy
    const fvPatchScalarField& kg =
        patch().lookupPatchField<volScalarField, scalar>(kName_);
    const immersedBoundaryWallFunctionFvPatchScalarField& kw =
        refCast<const immersedBoundaryWallFunctionFvPatchScalarField>(kg);

    // Current and new values of k at sampling point
    scalarField k = kw.ibSamplingPointValue(); 


    // Laminar viscosity
    const fvPatchScalarField& nuwg = 
        patch().lookupPatchField<volScalarField, scalar>(nuName_);

    const immersedBoundaryFvPatchScalarField& nu =
        refCast<const immersedBoundaryFvPatchScalarField>(nuwg);

    //scalarField nuw = nu.ibGhostCellValues();
    const scalarField& nuw = nu.ibValue();

    const scalar Cmu25 = pow025(Cmu_);

    tmp<scalarField> tnutw(new scalarField(nuw.size(), 0.0));

    scalarField& nutw = tnutw();

    forAll(nuw, cellI)
    {

        scalar yPlus = Cmu25*y[cellI]*sqrt(k[cellI])/nuw[cellI];
 
        if (yPlus > yPlusLam_)
        {
            nutw[cellI] = nuw[cellI]*(yPlus*kappa_/log(E_*yPlus) - 1.0);
        }

 
    }  
 

    tnutw = max(tnutw, SMALL);

     Info<< " nutw " << min(nutw) << " " << max(nutw) << " " << Foam::average(nutw) << endl;
     Info<< " nuEff " << min(nutw+nuw) << " " << max(nutw+nuw) << " " << Foam::average(nutw+nuw) << endl;
    return tnutw;

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

immersedBoundaryNutkWallFunctionFvPatchScalarField::immersedBoundaryNutkWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    immersedBoundaryWallFunctionFvPatchScalarField(p, iF),
    Cmu_(0.09),
    kappa_(0.41),
    E_(9.8)
{}


immersedBoundaryNutkWallFunctionFvPatchScalarField::immersedBoundaryNutkWallFunctionFvPatchScalarField
(
    const immersedBoundaryNutkWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    immersedBoundaryWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    Cmu_(ptf.Cmu_),
    kappa_(ptf.kappa_),
    E_(ptf.E_)
{}


immersedBoundaryNutkWallFunctionFvPatchScalarField::immersedBoundaryNutkWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    immersedBoundaryWallFunctionFvPatchScalarField(p, iF, dict),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09)),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    E_(dict.lookupOrDefault<scalar>("E", 9.8))
{}


immersedBoundaryNutkWallFunctionFvPatchScalarField::immersedBoundaryNutkWallFunctionFvPatchScalarField
(
    const immersedBoundaryNutkWallFunctionFvPatchScalarField& wfpsf
)
:
    immersedBoundaryWallFunctionFvPatchScalarField(wfpsf),
    Cmu_(wfpsf.Cmu_),
    kappa_(wfpsf.kappa_),
    E_(wfpsf.E_)
{}


immersedBoundaryNutkWallFunctionFvPatchScalarField::immersedBoundaryNutkWallFunctionFvPatchScalarField
(
    const immersedBoundaryNutkWallFunctionFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    immersedBoundaryWallFunctionFvPatchScalarField(wfpsf, iF),
    Cmu_(wfpsf.Cmu_),
    kappa_(wfpsf.kappa_),
    E_(wfpsf.E_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
scalar immersedBoundaryNutkWallFunctionFvPatchScalarField::yPlusLam
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

void immersedBoundaryNutkWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    scalarField& nutw = this->wallValue();

    nutw = calcNut();
    
 
    this->wallMask() = true;
 
    //immersedBoundaryWallFunctionFvPatchScalarField::updateCoeffs();
    this->setIbCellValues(this->wallValue());
    
    
    fvPatchScalarField::updateCoeffs();
    // Insert nut values into the internal field


}


void immersedBoundaryNutkWallFunctionFvPatchScalarField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    const label patchi = patch().index();
    word omegaName_("omega");
    word nutName_("nut");
    word nuName_("nu");
    word kName_("k");
    //this->setIbCellValues(this->ibSamplingPointValue());
 

    volScalarField& nutwg =
        const_cast<volScalarField&>
        (
            db().lookupObject<volScalarField>(nutName_)
        );
 

    volScalarField& nuwg =
        const_cast<volScalarField&>
        (
            db().lookupObject<volScalarField>(nuName_)
        );
 
    volScalarField& kg =
        const_cast<volScalarField&>
        (
            db().lookupObject<volScalarField>(kName_)
        );
 
    const scalar Cmu25 = pow(Cmu_, 0.25);

    const vectorField& an = ibPatch().adjacentIbNormals();
    const scalarField& aid = ibPatch().adjacentIbDelta();
    const labelList& aibc = ibPatch().adjacentIbCells();
    forAll (aibc, ibCellI)
    {
        scalar yplus=Cmu25*aid[ibCellI]*sqrt(kg[aibc[ibCellI]])/nuwg[aibc[ibCellI]];
        nutwg[aibc[ibCellI]]=0*nuwg[aibc[ibCellI]]*(yplus*kappa_/log(E_*yplus) - 1);
        
    }
    Info<<"evaluate nut"<<endl;
    fvPatchScalarField::evaluate(commsType);
}


void immersedBoundaryNutkWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    immersedBoundaryWallFunctionFvPatchScalarField::write(os);
    //writeLocalEntries(os);
    //writeEntry("value", os);
}

/*
tmp<scalarField> immersedBoundaryNutkWallFunctionFvPatchScalarField::yPlus() const
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

    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();
    tmp<scalarField> kwc = k.boundaryField()[patchi].patchInternalField();
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    return pow025(Cmu_)*y*sqrt(kwc)/nuw;
}

*/
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    immersedBoundaryNutkWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
