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

#include "immersedBoundaryOmegaWallFunctionFvPatchScalarField.H"
#include "immersedBoundaryVelocityWallFunctionFvPatchVectorField.H"
#include "turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

immersedBoundaryOmegaWallFunctionFvPatchScalarField::
immersedBoundaryOmegaWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    immersedBoundaryWallFunctionFvPatchScalarField(p, iF),
    UName_("U"),
    kName_("k"),
    GName_("kOmega:G"),
    nuName_("nu"),
    nutName_("nut"),
    Cmu_(0.09),
    kappa_(0.41),
    E_(9.8),
    beta1_(0.075)
{}


immersedBoundaryOmegaWallFunctionFvPatchScalarField::
immersedBoundaryOmegaWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    immersedBoundaryWallFunctionFvPatchScalarField(p, iF, dict),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    kName_(dict.lookupOrDefault<word>("k", "k")),
    GName_(dict.lookupOrDefault<word>("G", "kOmega:G")),
    nuName_(dict.lookupOrDefault<word>("nu", "nu")),
    nutName_(dict.lookupOrDefault<word>("nut", "nut")),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09)),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    E_(dict.lookupOrDefault<scalar>("E", 9.8)),
    beta1_(dict.lookupOrDefault<scalar>("beta1", 0.075))
{}


immersedBoundaryOmegaWallFunctionFvPatchScalarField::
immersedBoundaryOmegaWallFunctionFvPatchScalarField
(
    const immersedBoundaryOmegaWallFunctionFvPatchScalarField& ptf,
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
    E_(ptf.E_),
    beta1_(ptf.beta1_)
{}


immersedBoundaryOmegaWallFunctionFvPatchScalarField::
immersedBoundaryOmegaWallFunctionFvPatchScalarField
(
    const immersedBoundaryOmegaWallFunctionFvPatchScalarField& owfpsf
)
:
    immersedBoundaryWallFunctionFvPatchScalarField(owfpsf),
    UName_(owfpsf.UName_),
    kName_(owfpsf.kName_),
    GName_(owfpsf.GName_),
    nuName_(owfpsf.nuName_),
    nutName_(owfpsf.nutName_),
    Cmu_(owfpsf.Cmu_),
    kappa_(owfpsf.kappa_),
    E_(owfpsf.E_),
    beta1_(owfpsf.beta1_)
{}


immersedBoundaryOmegaWallFunctionFvPatchScalarField::
immersedBoundaryOmegaWallFunctionFvPatchScalarField
(
    const immersedBoundaryOmegaWallFunctionFvPatchScalarField& owfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    immersedBoundaryWallFunctionFvPatchScalarField(owfpsf, iF),
    UName_(owfpsf.UName_),
    kName_(owfpsf.kName_),
    GName_(owfpsf.GName_),
    nuName_(owfpsf.nuName_),
    nutName_(owfpsf.nutName_),
    Cmu_(owfpsf.Cmu_),
    kappa_(owfpsf.kappa_),
    E_(owfpsf.E_),
    beta1_(owfpsf.beta1_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void immersedBoundaryOmegaWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            dimensionedInternalField().group()
        )
    );

    GName_=turbModel.GName();

    // If G field is not present, execute zero gradient evaluation
    // HJ, 20/Mar/2011
    if (!db().foundObject<volScalarField>(GName_))
    {
        InfoIn
        (
            "void immersedBoundaryOmegaWallFunctionFvPatchScalarField::"
            "updateCoeffs()"
        )   << "Cannot access " << GName_ << " field.  for patch "
            << patch().name() << ".  Evaluating as regular immersed boundary"
            << endl;

        immersedBoundaryWallFunctionFvPatchScalarField::evaluate();

        return;
    }

    const vectorField& n = ibPatch().ibNormals();

    const scalar yPlusLam_ = yPlusLam(kappa_, E_);

    const scalar Cmu25 = pow(Cmu_, 0.25);
    const scalar Cmu50 = sqrt(Cmu_);

    volScalarField& G = const_cast<volScalarField&>
        (db().lookupObject<volScalarField>(GName_));


    // Grab values of other fields required for wall functions

    // Velocity
   const fvPatchVectorField& Uwg =
      patch().lookupPatchField<volVectorField, vector>(UName_);
    const immersedBoundaryVelocityWallFunctionFvPatchVectorField& Uw =
        refCast<const immersedBoundaryVelocityWallFunctionFvPatchVectorField>
        (
            Uwg
        );
 
    // Calculate tangential component, taking into account wall velocity
    const vectorField UtanOld =
        (I - sqr(n)) & (Uw.ibSamplingPointValue() - Uw.ibValue());
    const scalarField magUtanOld = mag(UtanOld);

    // Tangential velocity component
    scalarField& UTangentialNew = Uw.wallTangentialValue();

    // Wall shear stress
    vectorField& tauWall = Uw.tauWall();

    // Turbulence kinetic energy
    const fvPatchScalarField& kg =
        patch().lookupPatchField<volScalarField, scalar>(kName_);
    const immersedBoundaryWallFunctionFvPatchScalarField& kw =
        refCast<const immersedBoundaryWallFunctionFvPatchScalarField>(kg);

    // Current and new values of k at sampling point
    scalarField k = kw.ibSamplingPointValue();
    scalarField& kNew = kw.wallValue();
    //volScalarField& kPsi = const_cast<volScalarField&>
            //(G.mesh().lookupObject<volScalarField>(kName_));
 

    // Laminar viscosity
    const fvPatchScalarField& nuwg =
        patch().lookupPatchField<volScalarField, scalar>(nuName_);
    const immersedBoundaryFvPatchScalarField& nuw =
        refCast<const immersedBoundaryFvPatchScalarField>(nuwg);
    scalarField nu = nuw.ibCellValue();

    // Turbulent viscosity
    const fvPatchScalarField& nutwg =
       patch().lookupPatchField<volScalarField, scalar>(nutName_);
    const immersedBoundaryWallFunctionFvPatchScalarField& nutw =
        refCast<const immersedBoundaryWallFunctionFvPatchScalarField>(nutwg);

    // New values of nut
    //scalarField nutOld = nutw.ibCellValue();
    scalarField nutOld = nutw.ibSamplingPointValue();
    scalarField& nutNew = nutw.wallValue();

    const scalarField magGradUw = mag(Uw.ibGrad());

    // Get the IB addressing and distance
    const labelList& ibc = ibPatch().ibCells();

    // Distance to sampling point
    const scalarField& ySample = ibPatch().ibSamplingPointDelta();

    // Distance from wall to IB point
    const scalarField& y = ibPatch().ibDelta();

    // Omega: store IB cell values for direct insertion
    scalarField omegaSample = this->ibSamplingPointValue();

    scalarField& omegaNew = this->wallValue();

    // Mark values to be fixed
    boolList wf(ibc.size(), false);

    // Calculate yPlus for sample points
    scalarField ypd = Cmu25*ySample*sqrt(k)/nu;

     // Calculate wall function conditions
    forAll (ibc, ibCellI)
    {

        const scalar nuLam = nu[ibCellI];

        // Calculate yPlus from k and laminar viscosity for the IB point
 
        const scalar yPlusSample = ypd[ibCellI];

        scalar uTau;

        if (yPlusSample > yPlusLam_)
        {
            // Calculate uTau from log-law, knowing sampled k and U
            //uTau = magUtanOld[ibCellI]*kappa_/log(E_*yPlusSample);
            uTau = Cmu25*sqrt(k[ibCellI]);
            //uTau = yPlusSample*nuLam/ySample[ibCellI];
        }
        else
        {
            // Sampling point is in laminar sublayer
            // Bug fix: Xu, 21/Aug/2017
            uTau = magUtanOld[ibCellI]/yPlusSample;  

        }

        // Set wall shear stress
        tauWall[ibCellI] = sqr(uTau)*UtanOld[ibCellI]/(magUtanOld[ibCellI] + SMALL);

        // Calculate yPlus for IB point
        scalar yPlusIB = uTau*y[ibCellI]/nuLam;
        //scalar yPlusIB = Cmu25*y[ibCellI]*sqrt(kPsi[ibc[ibCellI]])/nuLam;
        //scalar yPlusIB = yPlusSample*y[ibCellI]/ySample[ibCellI];

        

        // Calculate wall function data in the immersed boundary point
        if (yPlusIB > yPlusLam_)
        {

            // Logarithmic region
            wf[ibCellI] = true;

            // turbulent viscosity at IB cell and at wall
            scalar nutwb = nuLam*(yPlusIB*kappa_/log(E_*yPlusIB) - 1);

            // Fix generation even though it if is not used
            G[ibc[ibCellI]] =
                sqr((nutwb + nuLam)*magGradUw[ibCellI])/
                (Cmu25*sqrt(k[ibCellI])*kappa_*y[ibCellI]);

            // Log-Law for tangential velocity
            UTangentialNew[ibCellI] =
                min
                (
                    magUtanOld[ibCellI],
                    //uTau/kappa_*log(E_*yPlusIB)
                    uTau*(1.0/kappa_*log(E_*yPlusIB))
                );

            // Calculate k in the IB cell from G = omega
            //kNew[ibCellI] = (nutwb + nuLam)*magGradUw[ibCellI]/Cmu50;
            kNew[ibCellI] = (nutwb + nuLam)*mag(UTangentialNew[ibCellI])/y[ibCellI]/Cmu50;
            //kNew[ibCellI] = 0.001;

            //kg[ibc[ibCellI]] = kNew[ibCellI];
            // Fix generation even though it if is not used
           /* G[ibc[ibCellI]] = 
                (nutwb + nuLam)*magGradUw[ibCellI]*Cmu25*sqrt(kNew[ibCellI])/(kappa_*y[ibCellI]);*/
 
            // Calculate turbulent viscosity
            nutNew[ibCellI] = nutwb;

            // Compute omega at the IB cell
            scalar omegaVis = 6.0*nutNew[ibCellI]/(beta1_*sqr(y[ibCellI]));
            scalar omegaLog = sqrt(k[ibCellI])/(Cmu25*kappa_*y[ibCellI]);

            omegaNew[ibCellI] = sqrt(sqr(omegaVis) + sqr(omegaLog));

        }
        else
        {
            // Laminar sub-layer
            wf[ibCellI] = false;

            // G is zero  - immaterial!
            // G[ibc[ibCellI]] = 0;

            // quadratic fit
            kNew[ibCellI] = k[ibCellI]*sqr(yPlusIB/yPlusLam_);

            // Laminar sub-layer for tangential velocity: uPlus = yPlus
            UTangentialNew[ibCellI] = uTau*yPlusIB;

            // Turbulent viscosity is zero
            nutNew[ibCellI] = SMALL;

            // Compute omega at the IB cell
            scalar omegaVis = 6.0*nutNew[ibCellI]/(beta1_*sqr(y[ibCellI]));
            scalar omegaLog = sqrt(k[ibCellI])/(Cmu25*kappa_*y[ibCellI]);

            omegaNew[ibCellI] = sqrt(sqr(omegaVis) + sqr(omegaLog));

            // Bugfix - set zeroGradient bc for large omega values at ib boundary
            // to avoid k unboundedness (IG 30/OCT/2015), not
            // sure if this is a good criteria
            if(omegaNew[ibCellI] > 10.0)
            {
                wf[ibCellI] = true;
            }
        }

    }
    scalarList UPout(5);
    scalarList nutPout(5);
    scalarList kPout(5);
    scalarList omegaPout(5);
    scalarList yPout(5);
    vectorList tauWallPout(4);

    Info<< "G " << gMin(G) << " " << gMax(G) << " " << gAverage(G) << endl;
 
    UPout[0] = gMin(UTangentialNew);
    nutPout[0] = gMin(nutNew);
    kPout[0] = gMin(kNew);
    omegaPout[0] = gMin(omegaNew);
    yPout[0] = gMin(ypd);
    tauWallPout[0] = gMin(tauWall);

    UPout[1] = gMax(UTangentialNew);
    nutPout[1] = gMax(nutNew);
    kPout[1] = gMax(kNew);
    omegaPout[1] = gMax(omegaNew);
    yPout[1] = gMax(ypd);
    tauWallPout[1] = gMax(tauWall);

    UPout[2] = gAverage(UTangentialNew);
    nutPout[2] = gAverage(nutNew);
    kPout[2] = gAverage(kNew);
    omegaPout[2] = gAverage(omegaNew);
    yPout[2] = gAverage(ypd);
    tauWallPout[2] = gAverage(tauWall);

    if (Pstream::parRun())//Start of mpi run
    {
        Info<< "UTangentialNew nutNew kNew omegaNew yPlus tauWall" << endl;
 
        for (label I = 0; I < 3; I++)
        {
            if(I==0){Info<<"min ";}
            if(I==1){Info<<"max ";}
            if(I==2){Info<<"ave ";}
            Info<<UPout[I]<<" "
                <<nutPout[I]<<" "
                <<kPout[I]<<" "
                <<omegaPout[I]<<" "
                <<yPout[I]<<" "
                <<tauWallPout[I]<<" "<<endl;
        }

    }//End of mpi run
    else//Start of serial run
    {
        Info<< "UTangentialNew nutNew kNew omegaNew yPlus tauWall" << endl;
        for (label I = 0; I < 3; I++)
        {
            if(I==0){Info<<"min ";}
            if(I==1){Info<<"max ";}
            if(I==2){Info<<"ave ";}
            Info<<UPout[I]<<" "
                <<nutPout[I]<<" "
                <<kPout[I]<<" "
                <<omegaPout[I]<<" "
                <<yPout[I]<<" "
                <<tauWallPout[I]<<" "<<endl;
        }
    }//End of serial run

    // Set the fields to calculated wall function values
    Uw.wallMask() = true;
    kw.wallMask() = wf;
    nutw.wallMask() = true;
    this->wallMask() = true;

    // Insert omega values into the internal field
    immersedBoundaryWallFunctionFvPatchScalarField::updateCoeffs();


}


void immersedBoundaryOmegaWallFunctionFvPatchScalarField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    // Insert omega values into the internal field
    this->setIbCellValues(this->wallValue());
 
    fvPatchScalarField::evaluate(commsType);
}


void immersedBoundaryOmegaWallFunctionFvPatchScalarField::
write(Ostream& os) const
{
    immersedBoundaryWallFunctionFvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "k", "k", kName_);
    writeEntryIfDifferent<word>(os, "G", "kOmega:G", GName_);
    writeEntryIfDifferent<word>(os, "nu", "nu", nuName_);
    writeEntryIfDifferent<word>(os, "nut", "nut", nutName_);
    os.writeKeyword("Cmu") << Cmu_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappa") << kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("E") << E_ << token::END_STATEMENT << nl;
    os.writeKeyword("beta1") << beta1_ << token::END_STATEMENT << nl;
}

scalar immersedBoundaryOmegaWallFunctionFvPatchScalarField::yPlusLam
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
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    immersedBoundaryOmegaWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
