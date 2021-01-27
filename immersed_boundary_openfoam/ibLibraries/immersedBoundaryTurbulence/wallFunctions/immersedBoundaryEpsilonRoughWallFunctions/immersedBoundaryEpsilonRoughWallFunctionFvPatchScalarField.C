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

#include "immersedBoundaryEpsilonRoughWallFunctionFvPatchScalarField.H"
#include "immersedBoundaryVelocityWallFunctionFvPatchVectorField.H"
 
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
 
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

immersedBoundaryEpsilonRoughWallFunctionFvPatchScalarField::
immersedBoundaryEpsilonRoughWallFunctionFvPatchScalarField
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
    Ks_(0),
    Cs_(0.5),
    E_(9.8)
{}


immersedBoundaryEpsilonRoughWallFunctionFvPatchScalarField::
immersedBoundaryEpsilonRoughWallFunctionFvPatchScalarField
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
    Ks_(dict.lookupOrDefault<scalar>("Ks",0.0)),
    Cs_(dict.lookupOrDefault<scalar>("Cs",0.5)),
    E_(dict.lookupOrDefault<scalar>("E", 9.8))

{}


immersedBoundaryEpsilonRoughWallFunctionFvPatchScalarField::
immersedBoundaryEpsilonRoughWallFunctionFvPatchScalarField
(
    const immersedBoundaryEpsilonRoughWallFunctionFvPatchScalarField& ptf,
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
    Ks_(ptf.Ks_),
    Cs_(ptf.Cs_),
    E_(ptf.E_)
{}


immersedBoundaryEpsilonRoughWallFunctionFvPatchScalarField::
immersedBoundaryEpsilonRoughWallFunctionFvPatchScalarField
(
    const immersedBoundaryEpsilonRoughWallFunctionFvPatchScalarField& ewfpsf
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
    Ks_(ewfpsf.Ks_),
    Cs_(ewfpsf.Cs_),
    E_(ewfpsf.E_)
{}


immersedBoundaryEpsilonRoughWallFunctionFvPatchScalarField::
immersedBoundaryEpsilonRoughWallFunctionFvPatchScalarField
(
    const immersedBoundaryEpsilonRoughWallFunctionFvPatchScalarField& ewfpsf,
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
    Ks_(ewfpsf.Ks_),
    Cs_(ewfpsf.Cs_),
    E_(ewfpsf.E_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void immersedBoundaryEpsilonRoughWallFunctionFvPatchScalarField::updateCoeffs()
{

    if (updated())
    {
        return;
    }
    GName_="ibmKEpsilon:G";
    // If G field is not present, execute zero gradient evaluation
    // HJ, 20/Mar/2011
    if (!db().foundObject<volScalarField>(GName_))
    {GName_="kEpsilon:G";}

    if (!db().foundObject<volScalarField>(GName_))
    {
        InfoIn
        (
            "void immersedBoundaryEpsilonRoughWallFunctionFvPatchScalarField::"
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
    const scalar Cmu75 = pow(Cmu_, 0.75);

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
    //        (G.mesh().lookupObject<volScalarField>(kName_));

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
    volScalarField& nutVol = const_cast<volScalarField&>
        (db().lookupObject<volScalarField>(nutName_));

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

    // Epsilon: store IB cell values for direct insertion
    scalarField epsilonSample = this->ibSamplingPointValue();
 
    scalarField& epsilonNew = this->wallValue();

    // Mark values to be fixed
    boolList wf(ibc.size(), false);

    // Calculate yPlus for sample points
    scalarField ypd = Cmu25*ySample*sqrt(k)/nu;
    scalarField ypIB(ypd.size(),0.0);

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
            //uTau = magUtanOld[ibCellI]/yPlusSample;  
            uTau = sqrt(nuLam*magUtanOld[ibCellI]/ySample[ibCellI]);
            //uTau=UTau[ibCellI];
        }
        // Calculate KsPlus
        scalar KsPlus = uTau*Ks_/nuLam;

        // Calculate Edash
        scalar Edash = E_;

        if (KsPlus > 2.25)
        {
            scalar fnRough=1.0 + Cs_*KsPlus;
            if(KsPlus < 90.0)
            {
                fnRough = pow
                (
                  (KsPlus - 2.25)/87.75 + Cs_*KsPlus,
                  sin(0.4258*(log(KsPlus) - 0.811))
                );
            }
            Edash /= fnRough;
        }

        // Set wall shear stress
        tauWall[ibCellI] = sqr(uTau)*UtanOld[ibCellI]/(magUtanOld[ibCellI] + SMALL);

        // Calculate yPlus for IB point
        scalar yPlusIB = uTau*y[ibCellI]/nuLam;
        //scalar yPlusIB = Cmu25*y[ibCellI]*sqrt(kPsi[ibc[ibCellI]])/nuLam;
        //scalar yPlusIB = yPlusSample*y[ibCellI]/ySample[ibCellI];
        ypIB[ibCellI] = yPlusIB;
 
        // Calculate wall function data in the immersed boundary point
        if (yPlusIB > yPlusLam_)
        {
            // Logarithmic region
            wf[ibCellI] = true;

           // To avoid oscillations limit the change in the wall viscosity
           // which is particularly important if it temporarily becomes zero
            scalar limitingNutw = max(nutVol[ibc[ibCellI]], nuLam);

            scalar nutwb = nuLam*(yPlusIB*kappa_/log(Edash*yPlusIB) - 1.0);

            nutwb = max ( min ( nutwb, 2*limitingNutw), 0.5*limitingNutw);

            // Fix generation even though it if is not used
            G[ibc[ibCellI]] =
                sqr((nutwb + nuLam)*magGradUw[ibCellI])/
                (Cmu25*sqrt(k[ibCellI])*kappa_*y[ibCellI]);

            // Log-Law for tangential velocity
            UTangentialNew[ibCellI] =
                min
                (
                    magUtanOld[ibCellI],
                    //uTau/kappa_*log(Edash*yPlusIB)
                    uTau*(1.0/kappa_*log(Edash*yPlusIB))
                );
            //UTangentialNew[ibCellI] = 0.2976;
            // Calculate turbulent viscosity
            nutNew[ibCellI] = nutwb;
            //nutNew[ibCellI] = 0.000002083;
 
            // Calculate k in the IB cell from G = epsilon
            //kNew[ibCellI] = (nutwb + nuLam)*magGradUw[ibCellI]/Cmu50;
            kNew[ibCellI] = (nutwb + nuLam)*mag(UTangentialNew[ibCellI])/y[ibCellI]/Cmu50;
            //kNew[ibCellI] = 0.001;

             /*G[ibc[ibCellI]] = 
                (nutwb + nuLam)*magGradUw[ibCellI]*
                Cmu25*sqrt(kNew[ibCellI])/(kappa_*y[ibCellI]);*/

            // Calculate epsilon from yPlus and set it
            epsilonNew[ibCellI] =
                //Cmu75*pow(k[ibCellI], 1.5)/(kappa_*y[ibCellI]);
                Cmu75*pow(kNew[ibCellI], 1.5)/(kappa_*y[ibCellI]);
            //epsilonNew[ibCellI] = 0.001;
        }
        else
        {
            // Laminar sub-layer
            wf[ibCellI] = false;

            // G is zero
            G[ibc[ibCellI]] = 0;
 
            // Laminar sub-layer for tangential velocity: uPlus = yPlus
            UTangentialNew[ibCellI] = min(magUtanOld[ibCellI], uTau*yPlusIB);
 
            // Turbulent viscosity is zero
            nutNew[ibCellI] = SMALL;

            // k is zero gradient: use the sampled value
            kNew[ibCellI] = k[ibCellI];

            // Set wall shear stress
            tauWall[ibCellI] = sqr(uTau)*UtanOld[ibCellI]/(magUtanOld[ibCellI] + SMALL);
 
            //scalar grad = (k[ibCellI]-kPsi[ibc[ibCellI]])/(ySample[ibCellI]-y[ibCellI]);
            //scalar grad = 0;
            //kNew[ibCellI] = k[ibCellI]-grad*ySample[ibCellI];
            // Calculate epsilon from yPlus and set it.
            // Note:  calculating equilibrium epsilon in the sub-layer creates
            //        an unrealistic oscillation: use sampled value
            // HJ, 27/Jul/2012
            epsilonNew[ibCellI] = epsilonSample[ibCellI];
            //epsilonNew[ibCellI] = 2.0*kNew[ibCellI]*nuLam/sqr(ySample[ibCellI]);
            //epsilonNew[ibCellI] =Cmu75*pow(k[ibCellI], 1.5)/(kappa_*y[ibCellI]);
             

        }
    }



    scalarList UPout(5);
    scalarList nutPout(5);
    scalarList kPout(5);
    scalarList epsilonPout(5);
    scalarList yPout(5);
    scalarList yIBPout(5);
    vectorList tauWallPout(4);

    Info<< "G " << gMin(G) << " " << gMax(G) << " " << gAverage(G) << endl;
 
    UPout[0] = gMin(UTangentialNew);
    nutPout[0] = gMin(nutNew);
    kPout[0] = gMin(kNew);
    epsilonPout[0] = gMin(epsilonNew);
    yPout[0] = gMin(ypd);
    yIBPout[0] = gMin(ypIB);
    tauWallPout[0] = gMin(tauWall);

    UPout[1] = gMax(UTangentialNew);
    nutPout[1] = gMax(nutNew);
    kPout[1] = gMax(kNew);
    epsilonPout[1] = gMax(epsilonNew);
    yPout[1] = gMax(ypd);
    yIBPout[1] = gMax(ypIB);
    tauWallPout[1] = gMax(tauWall);

    UPout[2] = gAverage(UTangentialNew);
    nutPout[2] = gAverage(nutNew);
    kPout[2] = gAverage(kNew);
    epsilonPout[2] = gAverage(epsilonNew);
    yPout[2] = gAverage(ypd);
    yIBPout[2] = gAverage(ypIB);
    tauWallPout[2] = gAverage(tauWall);

    if (Pstream::parRun())//Start of mpi run
    {
        Info<< "UTangentialNew nutNew kNew epsilonNew yPlus yPlusIB tauWall" << endl;
 
        for (label I = 0; I < 3; I++)
        {
            if(I==0){Info<<"min ";}
            if(I==1){Info<<"max ";}
            if(I==2){Info<<"ave ";}
            Info<<UPout[I]<<" "
                <<nutPout[I]<<" "
                <<kPout[I]<<" "
                <<epsilonPout[I]<<" "
                <<yPout[I]<<" "
                <<yIBPout[I]<<" "
                <<tauWallPout[I]<<" "<<endl;
        }

    }//End of mpi run
    else//Start of serial run
    {
        Info<< "UTangentialNew nutNew kNew epsilonNew yPlus yPlusIB tauWall" << endl;
        for (label I = 0; I < 3; I++)
        {
            if(I==0){Info<<"min ";}
            if(I==1){Info<<"max ";}
            if(I==2){Info<<"ave ";}
            Info<<UPout[I]<<" "
                <<nutPout[I]<<" "
                <<kPout[I]<<" "
                <<epsilonPout[I]<<" "
                <<yPout[I]<<" "
                <<yIBPout[I]<<" "
                <<tauWallPout[I]<<" "<<endl;
        }
    }//End of serial run

    // Set the fields to calculated wall function values
    Uw.wallMask() = true;
    kw.wallMask() = wf;
    nutw.wallMask() = true;
    this->wallMask() = true;

    // Insert epsilon values into the internal field
    immersedBoundaryWallFunctionFvPatchScalarField::updateCoeffs();
}


void immersedBoundaryEpsilonRoughWallFunctionFvPatchScalarField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    // Insert epsilon values into the internal field
    this->setIbCellValues(this->wallValue());

    fvPatchScalarField::evaluate(commsType);
}


void immersedBoundaryEpsilonRoughWallFunctionFvPatchScalarField::
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
    os.writeKeyword("Ks") << Ks_ << token::END_STATEMENT << nl;
}

scalar immersedBoundaryEpsilonRoughWallFunctionFvPatchScalarField::yPlusLam
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
    immersedBoundaryEpsilonRoughWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

 
} // End namespace Foam

// ************************************************************************* //
