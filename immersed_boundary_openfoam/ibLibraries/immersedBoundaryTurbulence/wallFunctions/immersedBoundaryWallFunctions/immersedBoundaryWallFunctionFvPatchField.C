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

#include "immersedBoundaryWallFunctionFvPatchField.H"
#include "fvPatchFieldMapper.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
 

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //



template<class Type>
void immersedBoundaryWallFunctionFvPatchField<Type>::setIbCellValues
(
    const Field<Type>& ibcValues
) const
{
    const labelList& ibc = this->ibPatch().ibCells();
       //Info<< "setIbCellValues  field type "<<   this->dimensionedInternalField().name()<<endl;
      //Info<< "ibcValues  " << min(ibcValues) << " " << max(ibcValues) << " " << average(ibcValues) << endl;
      //Info<< "wallValue_  " << min(wallValue_) << " " << max(wallValue_) << " " << average(wallValue_) << endl;
    if (ibcValues.size() != ibc.size())
    {
        FatalErrorIn
        (
            "template<class Type>\n"
            "void immersedBoundaryWallFunctionFvPatchField<Type>::"
            "setIbCellValues\n"
            "(\n"
            "    const Field<Type>& ibcValues\n"
            ") const"
        )   << "Size of ibcValues field not equal to the number of IB cells."
            << nl << "ibcValues: " << ibcValues.size()
            << " ibc: " << ibc.size()
            << abort(FatalError);
    }

    // Get non-const access to internal field
    Field<Type>& psiI = const_cast<Field<Type>&>(this->internalField());

    //const Field<Type>& spv = this->ibSamplingPointValue();
    if (wallValue_.empty() || wallMask_.empty())
    {
        immersedBoundaryFvPatchField<Type>::setIbCellValues(ibcValues);
    }
    else
    {
        forAll (ibcValues, cellI)
        {
            // If mask is set use the wall value, otherwise use the
            // fitted value
            if (wallMask_[cellI])
            {
                psiI[ibc[cellI]] = wallValue_[cellI];
            }
            else
            {
                psiI[ibc[cellI]] = ibcValues[cellI];
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
immersedBoundaryWallFunctionFvPatchField<Type>::
immersedBoundaryWallFunctionFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    immersedBoundaryFvPatchField<Type>(p, iF),
    wallValue_(),
    wallMask_()
{}


template<class Type>
immersedBoundaryWallFunctionFvPatchField<Type>::
immersedBoundaryWallFunctionFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    immersedBoundaryFvPatchField<Type>(p, iF, dict),
    wallValue_(),
    wallMask_()
{}


template<class Type>
immersedBoundaryWallFunctionFvPatchField<Type>::
immersedBoundaryWallFunctionFvPatchField
(
    const immersedBoundaryWallFunctionFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    immersedBoundaryFvPatchField<Type>(ptf, p, iF, mapper),
    wallValue_(),
    wallMask_()
{}


template<class Type>
immersedBoundaryWallFunctionFvPatchField<Type>::
immersedBoundaryWallFunctionFvPatchField
(
    const immersedBoundaryWallFunctionFvPatchField& tkqrwfpf
)
:
    immersedBoundaryFvPatchField<Type>(tkqrwfpf),
    wallValue_(),
    wallMask_()
{}


template<class Type>
immersedBoundaryWallFunctionFvPatchField<Type>::
immersedBoundaryWallFunctionFvPatchField
(
    const immersedBoundaryWallFunctionFvPatchField& tkqrwfpf,
    const DimensionedField<Type, volMesh>& iF
)
:
    immersedBoundaryFvPatchField<Type>(tkqrwfpf, iF),
    wallValue_(),
    wallMask_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<Field<Type> >
immersedBoundaryWallFunctionFvPatchField<Type>::ibNewSamplingPointValue() const
{
    Info<<"ibNewSamplingPointValue"<<endl;
    const immersedBoundaryFvPatch& ibPatch_ = this->ibPatch();
    const labelList& ibc = ibPatch_.ibCells();
    tmp<Field<Type> > tIbPsi
    (
        new Field<Type>(ibc.size(), pTraits<Type>::zero)
    );
    Field<Type>& ibPsi = tIbPsi();

    word name = this->dimensionedInternalField().name();
    //if(name!="U" and name!="nut" and name!="p" and name!="omega" and name!="epsilon")
    //if(name=="U")
    {
    const Field<Type>& cellValues = this->internalField();

    const Field<Type>& ibValues = this->wallValue();
Info<<"ibValues[10] "<<ibValues[10]<<endl;
    // Get addressing
    const labelListList& ibcc = ibPatch_.ibCellCells();
    const labelListList& ibcp = ibPatch_.ibCellPts();
    const labelListList& procIbCells = ibPatch_.ibProcIbCells();
    const labelListList& procCells = ibPatch_.ibProcCells();
    const List<List<labelPair> >& ibcProcC = ibPatch_.ibCellProcCells();
    const List<List<labelPair> >& ibcProcP = ibPatch_.ibCellProcPts();

    // Get weightsw with corresponding IB points
    const scalarListList& cellWeights = ibPatch_.ibNewSamplingWeights();
    const scalarListList& cellProcWeights = ibPatch_.ibNewSamplingProcWeights();
    const scalarListList& ptsWeights = ibPatch_.ibNewPtsWeights();
    const scalarListList& ptsProcWeights = ibPatch_.ibNewPtsProcWeights();


    // Do interpolation, local cell data
    forAll (ibc, cellI)
    {
        const labelList& curAddr = ibcc[cellI];
        const scalarList& curWeights = cellWeights[cellI];

        // for Cells
        forAll (curAddr, ccI)
        {
            ibPsi[cellI] += curWeights[ccI]*cellValues[curAddr[ccI]];
        }

        // for IB points
        const labelList& curCellPts = ibcp[cellI];
        const scalarList& curPW = ptsWeights[cellI];

        forAll (curCellPts, ccI)
        {
            ibPsi[cellI] += curPW[ccI]*ibValues[curCellPts[ccI]];
    //if(cellI==10){Info<< this->dimensionedInternalField().name()<<" ibValues[curPtPts[ccI]] "<<ibValues[curPtPts[ccI]]<<endl;}
        }
    }

    // Create CellValues in different processors
    FieldField<Field, Type> procCellValues = ibPatch_.sendAndReceive(cellValues,procCells);
    FieldField<Field, Type> procIbValues = ibPatch_.sendAndReceive(cellValues,procIbCells);

    // Do interpolation, cell data from other processors
    forAll (ibc, cellI)
    {
        const List<labelPair>& curProcCells = ibcProcC[cellI];
        const scalarList& curProcWeights = cellProcWeights[cellI];

        // for Cells
        forAll (curProcCells, cpcI)
        {
            ibPsi[cellI] +=
                curProcWeights[cpcI]*
                procCellValues
                [
                    curProcCells[cpcI].first()
                ]
                [
                    curProcCells[cpcI].second()
                ];
        }


        const List<labelPair>& curProcPts = ibcProcP[cellI];
        const scalarList& curPtsProcWeights = ptsProcWeights[cellI];

        // for IB points
        forAll (curProcPts, cpcI)
        {
            ibPsi[cellI] +=
                curPtsProcWeights[cpcI]*
                procIbValues
                [
                    curProcPts[cpcI].first()
                ]
                [
                    curProcPts[cpcI].second()
                ];
        }
    }

    }
 
    return tIbPsi;
    //return ibPatch_.toSamplingPoints(this->internalField());
}


template<class Type>
void immersedBoundaryWallFunctionFvPatchField<Type>::updateCoeffs()
{
    //if (this->updated())
    //{
        //return;
    //}
    //this->setIbCellValues(this->wallValue());

    //fvPatchField<Type>::updateCoeffs();
    //Pout<<"immersedBoundaryWallFunctionFvPatchField "<< this->dimensionedInternalField().name()<<endl;
    //const immersedBoundaryFvPatch& ibFvP =
        //immersedBoundaryFvPatchField<Type>::ibPatch();

	//this->refGrad()=ibFvP.toTriFaces(this->wallGrad());
    immersedBoundaryFvPatchField<Type>::updateCoeffs();

}


template<class Type>
void immersedBoundaryWallFunctionFvPatchField<Type>::evaluate()
{
    //this->setIbCellValues(this->wallValue());
    immersedBoundaryFvPatchField<Type>::evaluate();

}

template<class Type>
Foam::Field<Type>& immersedBoundaryWallFunctionFvPatchField<Type>::wallValue() const
{
    // Note: on a moving mesh, the intersection has changed and
    // wallValue fields should be cleared and recalculated.
    // This should happen only once, but I cannot see the mechanism
    // HJ, 17/Oct/2012
    // Bugfix 30/OCT/2015 - check if the mesh is moving

    //const immersedBoundaryFvPatch& ibFvP =
        //immersedBoundaryFvPatchField<Type>::ibPatch();

    if
    (
        wallValue_.empty()
    // || (ibFvP.movingIb() || ibFvP.boundaryMesh().mesh().moving())
    )
    {
		wallValue_.clear();
        wallValue_.setSize
        (
            this->ibPatch().ibCells().size(),
            pTraits<Type>::zero
        );

        
    }

    return wallValue_;
}


template<class Type>
Foam::boolList& immersedBoundaryWallFunctionFvPatchField<Type>::wallMask() const
{
    // Note: on a moving mesh, the intersection has changed and
    // wallValue fields should be cleared and recalculated.
    // This should happen only once, but I cannot see the mechanism
    // HJ, 17/Oct/2012
    // Bugfix 30/OCT/2015 - check if the mesh is moving
    const immersedBoundaryFvPatch& ibFvP =
        immersedBoundaryFvPatchField<Type>::ibPatch();

    if
    (
        wallMask_.empty()
     || (ibFvP.movingIb() || ibFvP.boundaryMesh().mesh().moving())
    )
    {
        wallMask_.setSize
        (
            this->ibPatch().ibCells().size(),
            false
        );
    }

    return wallMask_;
}

template<class Type>
Foam::Field<Type>& immersedBoundaryWallFunctionFvPatchField<Type>::wallGrad() const
{
    // Note: on a moving mesh, the intersection has changed and
    // wallValue fields should be cleared and recalculated.
    // This should happen only once, but I cannot see the mechanism
    // HJ, 17/Oct/2012
    // Bugfix 30/OCT/2015 - check if the mesh is moving

    //const immersedBoundaryFvPatch& ibFvP =
        //immersedBoundaryFvPatchField<Type>::ibPatch();

    if
    (
        wallGrad_.empty()
//     || (ibFvP.movingIb() || ibFvP.boundaryMesh().mesh().moving())
    )
    {
        wallGrad_.setSize
        (
            this->ibPatch().ibCells().size(),
            pTraits<Type>::zero
        );

        
    }

    return wallGrad_;
}


template<class Type>
void immersedBoundaryWallFunctionFvPatchField<Type>::setIbCellValuesRef
(
    const Field<Type>& ibcValues
) const
{
    this->setIbCellValues(ibcValues);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
} // End namespace Foam

// ************************************************************************* //
