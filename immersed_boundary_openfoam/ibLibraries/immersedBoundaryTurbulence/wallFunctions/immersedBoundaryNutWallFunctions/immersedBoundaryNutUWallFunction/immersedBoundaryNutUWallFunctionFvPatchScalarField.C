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

#include "immersedBoundaryNutUWallFunctionFvPatchScalarField.H"
//#include "turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField> immersedBoundaryNutUWallFunctionFvPatchScalarField::calcNut() const
{
 

    word UName_("U");
    word nuName_("nu");
    word kName_("k");
    // Velocity
    const fvPatchVectorField& Uwg =
      patch().lookupPatchField<volVectorField, vector>(UName_);

    const immersedBoundaryFvPatchVectorField& Uw =
        refCast<const immersedBoundaryFvPatchVectorField>
        (
            Uwg
        );
    const vectorField& n = ibPatch().ibNormals();

    const scalarField magUp(mag((I - sqr(n)) & (Uw.ibSamplingPointValue()-Uw.ibValue())));


    // Distance from wall to IB point
    const scalarField& y = ibPatch().ibDelta();

    const scalar yPlusLam_ = yPlusLam(kappa_, E_);

    // Laminar viscosity
    const fvPatchScalarField& nuwg = 
        patch().lookupPatchField<volScalarField, scalar>(nuName_);

    const immersedBoundaryFvPatchScalarField& nu =
        refCast<const immersedBoundaryFvPatchScalarField>(nuwg);

    const scalarField& nuw = nu.ibValue();


    // Calculate yPlus using U
    tmp<scalarField> tyPlus(new scalarField(y.size(), 0.0));
    scalarField& yPlus = tyPlus();

    forAll(yPlus, cellI)
    {

        scalar kappaRe = kappa_*magUp[cellI]*y[cellI]/nuw[cellI];

        scalar yp = yPlusLam_;
        scalar ryPlusLam = 1.0/yp;

        int iter = 0;
        scalar yPlusLast = 0.0;

        do
        {
            yPlusLast = yp; 
            yp = (kappaRe + yp)/(1.0 + log(E_*yp));

        } while (mag(ryPlusLam*(yp - yPlusLast)) > 0.01 && ++iter < 10 );
        
        yPlus[cellI] = max(0.0, yp);
    }


    // Turbulence kinetic energy
    const fvPatchScalarField& kg =
        patch().lookupPatchField<volScalarField, scalar>(kName_);
    const immersedBoundaryWallFunctionFvPatchScalarField& kw =
        refCast<const immersedBoundaryWallFunctionFvPatchScalarField>(kg);

    // Current and new values of k at sampling point
    scalarField k = kw.wallValue(); 




    const scalar Cmu25 = pow025(Cmu_);

    tmp<scalarField> tnutw(new scalarField(nuw.size(), 0.0));

    scalarField& nutw = tnutw();

    forAll(nuw, cellI)
    {

        //scalar yPlus = Cmu25*y[cellI]*sqrt(k[cellI])/nuw[cellI];
//Info<<cellI<<" old nutw[cellI] "<<nutw[cellI]<<" nuw[cellI] "<<nuw[cellI]<< " yPlus "<<yPlus<<" k[cellI] "<<k[cellI]<< endl;
        if (yPlus[cellI] > yPlusLam_)
        {
            nutw[cellI] = nuw[cellI]*(yPlus[cellI]*kappa_/log(E_*yPlus[cellI]) - 1.0);
        }

//Info<<cellI<<" nutw[cellI] "<<nutw[cellI]<<" nuw[cellI] "<<nuw[cellI]<<endl;
    } //Info<<tnutw<<endl;
    //tnutw=nutw;

     scalarField magGradUw=magUp/y;
     Info<< "yPlus " << min(yPlus) << " " << max(yPlus) << " " << Foam::average(yPlus) << endl;
     Info<< "nutw " << min(nutw) << " " << max(nutw) << " " << Foam::average(nutw) << endl;
     Info<< "magGradUw " << min(magGradUw) << " " << max(magGradUw) << " " << Foam::average(magGradUw) << endl;
     Info<< "shieldsNumber " << min(magGradUw*(nutw+nuw)/9.8/0.00048) << " " << max(magGradUw*(nutw+nuw)/9.8/0.00048) << " " << Foam::average(magGradUw*(nutw+nuw)/9.8/0.00048) << endl;
    return tnutw;

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

immersedBoundaryNutUWallFunctionFvPatchScalarField::immersedBoundaryNutUWallFunctionFvPatchScalarField
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


immersedBoundaryNutUWallFunctionFvPatchScalarField::immersedBoundaryNutUWallFunctionFvPatchScalarField
(
    const immersedBoundaryNutUWallFunctionFvPatchScalarField& ptf,
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


immersedBoundaryNutUWallFunctionFvPatchScalarField::immersedBoundaryNutUWallFunctionFvPatchScalarField
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


immersedBoundaryNutUWallFunctionFvPatchScalarField::immersedBoundaryNutUWallFunctionFvPatchScalarField
(
    const immersedBoundaryNutUWallFunctionFvPatchScalarField& wfpsf
)
:
    immersedBoundaryWallFunctionFvPatchScalarField(wfpsf),
    Cmu_(wfpsf.Cmu_),
    kappa_(wfpsf.kappa_),
    E_(wfpsf.E_)
{}


immersedBoundaryNutUWallFunctionFvPatchScalarField::immersedBoundaryNutUWallFunctionFvPatchScalarField
(
    const immersedBoundaryNutUWallFunctionFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    immersedBoundaryWallFunctionFvPatchScalarField(wfpsf, iF),
    Cmu_(wfpsf.Cmu_),
    kappa_(wfpsf.kappa_),
    E_(wfpsf.E_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
scalar immersedBoundaryNutUWallFunctionFvPatchScalarField::yPlusLam
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

void immersedBoundaryNutUWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    scalarField& nutw = this->wallValue();

    nutw = calcNut();

 
    this->wallMask() = true;

    this->updateGhostCellValues();
    Info<<"updateCoeffs nut"<<endl;
    immersedBoundaryWallFunctionFvPatchScalarField::updateCoeffs();
 
}


void immersedBoundaryNutUWallFunctionFvPatchScalarField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    // Insert omega values into the internal field
    //this->setIbCellValues(this->wallValue());
Info<<"evaluate nut"<<endl;
    fvPatchScalarField::evaluate(commsType);
}


void immersedBoundaryNutUWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    immersedBoundaryWallFunctionFvPatchScalarField::write(os);
    //writeLocalEntries(os);
    //writeEntry("value", os);
}

/*
tmp<scalarField> immersedBoundaryNutUWallFunctionFvPatchScalarField::yPlus() const
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
    immersedBoundaryNutUWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
