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

#include "immersedBoundaryVelocityWallFunctionFvPatchVectorField.H"
#include "immersedBoundaryWallFunctionFvPatchFields.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
 
void immersedBoundaryVelocityWallFunctionFvPatchVectorField::setIbCellValues
(
    const vectorField& ibcValues
) const
{
    const labelList& ibc = ibPatch().ibCells();
    //Info<< "setIbCellValues  field type "<<   this->dimensionedInternalField().name()<<endl;


    if (ibcValues.size() != ibc.size())
    {
        FatalErrorIn
        (
            "void immersedBoundaryVelocityWallFunctionFvPatchVectorField::"
            "setIbCellValues\n"
            "(\n"
            "    const vectorField& ibcValues\n"
            ") const"
        )   << "Size of ibcValues field not equal to the number of IB cells."
            << nl << "ibcValues: " << ibcValues.size()
            << " ibc: " << ibc.size()
            << abort(FatalError);
    }

    // Get non-const access to internal field
    vectorField& psiI = const_cast<vectorField&>(this->internalField());
 
    immersedBoundaryFvPatchVectorField::setIbCellValues(ibcValues);

    //wallValue = ibcValues;

    if (wallTangentialValue_.empty() || wallMask_.empty())
    {

        immersedBoundaryFvPatchVectorField::setIbCellValues(ibcValues);
    }
    else
    {

        const vectorField& n = ibPatch().ibNormals();

        // Calculate tangential component taking into account wall velocity
        //scalarField UtanOld = mag((I - sqr(n)) & this->ibCellValue());

        vectorField Uwall = this->ibValue();

        forAll (ibcValues, cellI)
        {
            // If mask is set, correct the velocity for the
            //  tangential wall value, otherwise use the fitted value
            if (wallMask_[cellI])
            {
                // Decompose fitted velocity into the normal and
                // tangential components
                const vector& curN = n[cellI];
                const vector curU = psiI[ibc[cellI]];

                scalar ibcNormal = curN & ibcValues[cellI];

                // Get tangential velocity and direction
                vector ibcTangential = (I - sqr(curN)) & curU;
                ibcTangential /= mag(ibcTangential) + SMALL;

                // Reconstruct the velocity, imposing the magnitude of
                // tangential value and add wall velocity
                psiI[ibc[cellI]] = curN*ibcNormal
                    + ibcTangential*wallTangentialValue_[cellI]
                    + Uwall[cellI];

             }
            else
            {
                psiI[ibc[cellI]] = ibcValues[cellI];

            }
            //wallValue[cellI] = psiI[ibc[cellI]];
        }
    }


}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

immersedBoundaryVelocityWallFunctionFvPatchVectorField::
immersedBoundaryVelocityWallFunctionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    immersedBoundaryWallFunctionFvPatchVectorField(p, iF),
    wallTangentialValue_(),
    tauWall_(),
    wallMask_(),
    depositionRate_(),
    entrainmentRate_(),
    sedFlux_()
{}


immersedBoundaryVelocityWallFunctionFvPatchVectorField::
immersedBoundaryVelocityWallFunctionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    immersedBoundaryWallFunctionFvPatchVectorField(p, iF, dict),
    wallTangentialValue_(),
    tauWall_(),
    wallMask_(),
    depositionRate_(),
    entrainmentRate_(),
    sedFlux_()
{}


immersedBoundaryVelocityWallFunctionFvPatchVectorField::
immersedBoundaryVelocityWallFunctionFvPatchVectorField
(
    const immersedBoundaryVelocityWallFunctionFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    immersedBoundaryWallFunctionFvPatchVectorField(ptf, p, iF, mapper),
    wallTangentialValue_(),
    tauWall_(),
    wallMask_(),
    depositionRate_(),
    entrainmentRate_(),
    sedFlux_()
{}


immersedBoundaryVelocityWallFunctionFvPatchVectorField::
immersedBoundaryVelocityWallFunctionFvPatchVectorField
(
    const immersedBoundaryVelocityWallFunctionFvPatchVectorField& ewfpsf
)
:
    immersedBoundaryWallFunctionFvPatchVectorField(ewfpsf),
    wallTangentialValue_(),
    tauWall_(),
    wallMask_(),
    depositionRate_(),
    entrainmentRate_(),
    sedFlux_()
{}


immersedBoundaryVelocityWallFunctionFvPatchVectorField::
immersedBoundaryVelocityWallFunctionFvPatchVectorField
(
    const immersedBoundaryVelocityWallFunctionFvPatchVectorField& ewfpsf,
    const DimensionedField<vector, volMesh>& iF
)
:
    immersedBoundaryWallFunctionFvPatchVectorField(ewfpsf, iF),
    wallTangentialValue_(),
    tauWall_(),
    wallMask_(),
    depositionRate_(),
    entrainmentRate_(),
    sedFlux_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void immersedBoundaryVelocityWallFunctionFvPatchVectorField::updateCoeffs()
{
        immersedBoundaryFvPatchVectorField::updateCoeffs();
}

const vectorField&
immersedBoundaryVelocityWallFunctionFvPatchVectorField::wallShearStress() const
{

    if (tauWall_.empty())
    {
        FatalErrorIn
        (
            "const vectorField& "
            "immersedBoundaryVelocityWallFunctionFvPatchVectorField::"
            "wallShearStress() const"
        )   << "tauWall not set for IB patch " << patch().name()
            << " for field " << dimensionedInternalField().name()
            << abort(FatalError);
    }

    return tauWall_;
}

Foam::scalarField&
immersedBoundaryVelocityWallFunctionFvPatchVectorField::wallTangentialValue() const
{
    // Bugfix 30/OCT/2015 - check if the mesh is moving

    const immersedBoundaryFvPatch& ibFvP =
        immersedBoundaryFvPatchVectorField::ibPatch();

    if
    (
        wallTangentialValue_.empty()
     || (ibFvP.movingIb() || ibFvP.boundaryMesh().mesh().moving())
    )
    {
		wallTangentialValue_.clear();
        wallTangentialValue_.setSize
        (
            this->ibPatch().ibCells().size(),
            0
        );
	}

    return wallTangentialValue_;
}


Foam::vectorField&
immersedBoundaryVelocityWallFunctionFvPatchVectorField::tauWall() const
{
    // Bugfix 30/OCT/2015 - check if the mesh is moving

//    const immersedBoundaryFvPatch& ibFvP =
//        immersedBoundaryFvPatchVectorField::ibPatch();
    if
    (
        tauWall_.empty()
//     || ibFvP.movingIb() || ibFvP.boundaryMesh().mesh().moving()
    )
    {
		tauWall_.clear();
        tauWall_.setSize
        (
            this->ibPatch().ibCells().size(),
            vector::zero
        );
    }

    return tauWall_;
}


Foam::boolList&
immersedBoundaryVelocityWallFunctionFvPatchVectorField::wallMask() const
{
    // Bugfix 30/OCT/2015 - check if the mesh is moving

    const immersedBoundaryFvPatch& ibFvP =
        immersedBoundaryFvPatchVectorField::ibPatch();
    if
    (
        wallMask_.empty()
     || (ibFvP.movingIb() || ibFvP.boundaryMesh().mesh().moving())
    )
    {
		wallMask_.clear();
        wallMask_.setSize
        (
            this->ibPatch().ibCells().size(),
            false
        );
    }

    return wallMask_;
}

Foam::scalarField&
immersedBoundaryVelocityWallFunctionFvPatchVectorField::depositionRate() const
{
    // Bugfix 30/OCT/2015 - check if the mesh is moving

    //const immersedBoundaryFvPatch& ibFvP =
        //immersedBoundaryFvPatchVectorField::ibPatch();
    if
    (
        depositionRate_.empty()
     //|| (ibFvP.movingIb() || ibFvP.boundaryMesh().mesh().moving())
    )
    {
		depositionRate_.clear();
        depositionRate_.setSize
        (
            this->ibPatch().ibCells().size(),
            0
        );
    }

    return depositionRate_;
}

Foam::scalarField&
immersedBoundaryVelocityWallFunctionFvPatchVectorField::entrainmentRate() const
{
    // Bugfix 30/OCT/2015 - check if the mesh is moving

    //const immersedBoundaryFvPatch& ibFvP =
        //immersedBoundaryFvPatchVectorField::ibPatch();
    if
    (
        entrainmentRate_.empty()
     //|| (ibFvP.movingIb() || ibFvP.boundaryMesh().mesh().moving())
    )
    {
		entrainmentRate_.clear();
        entrainmentRate_.setSize
        (
            this->ibPatch().ibCells().size(),
            0
        );
    }

    return entrainmentRate_;
}

Foam::scalarField&
immersedBoundaryVelocityWallFunctionFvPatchVectorField::sedFlux() const
{
    // Bugfix 30/OCT/2015 - check if the mesh is moving

    const immersedBoundaryFvPatch& ibFvP =
        immersedBoundaryFvPatchVectorField::ibPatch();
    if
    (
        sedFlux_.empty()
     || (ibFvP.movingIb() || ibFvP.boundaryMesh().mesh().moving())
    )
    {
		sedFlux_.clear();
        sedFlux_.setSize
        (
            this->ibPatch().ibCells().size(),
            false
        );
    }

    return sedFlux_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    immersedBoundaryVelocityWallFunctionFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
} // End namespace Foam

// ************************************************************************* //
