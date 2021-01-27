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

#include "immersedBoundaryFvPatch.H"
#include "fvMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::immersedBoundaryFvPatch::makeIbNewSamplingWeights() const
{
    if (debug)
    {
        Info<< "immersedBoundaryFvPatch::makeIbNewSamplingWeights() : "
            << "making sampling point weights"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (
           ibNewSamplingWeightsPtr_ 
        || ibNewSamplingProcWeightsPtr_
        || ibNewPtsWeightsPtr_
        || ibNewPtsProcWeightsPtr_
       )
    {
        FatalErrorIn("void immersedBoundaryFvPatch::makeIbNewSamplingWeights()")
            << "sampling point weights already exist"
            << abort(FatalError);
    }

    // Get addressing
    const labelList& ibc = ibCells();
    const labelListList& ibcc = ibCellCells();
    const labelListList& ibcp = ibCellPts();
    const List<List<labelPair> >& ibcProcC = ibCellProcCells();
    const List<List<labelPair> >& ibcProcP = ibCellProcPts();

    // Initialise the weights
    ibNewSamplingWeightsPtr_ = new scalarListList(ibc.size());
    scalarListList& cellWeights = *ibNewSamplingWeightsPtr_;

    forAll (cellWeights, cellI)
    {
        cellWeights[cellI].setSize(ibcc[cellI].size(), 0);
    }

    ibNewSamplingProcWeightsPtr_ = new scalarListList(ibc.size());
    scalarListList& cellProcWeights = *ibNewSamplingProcWeightsPtr_;

    forAll (cellProcWeights, cellI)
    {
        cellProcWeights[cellI].setSize(ibcProcC[cellI].size(), 0);
    }

    ibNewPtsWeightsPtr_ = new scalarListList(ibc.size());
    scalarListList& ptsWeights = *ibNewPtsWeightsPtr_;

    forAll (ptsWeights, cellI)
    {
        ptsWeights[cellI].setSize(ibcp[cellI].size(), 0);
    }

    ibNewPtsProcWeightsPtr_ = new scalarListList(ibc.size());
    scalarListList& ptsProcWeights = *ibNewPtsProcWeightsPtr_;

    forAll (ptsProcWeights, cellI)
    {
        ptsProcWeights[cellI].setSize(ibcProcP[cellI].size(), 0);
    }

    // Get sampling point locations
    const vectorField& samplingPoints = ibSamplingPoints();
    const scalarField& gammaIn = gamma().internalField();
    const vectorField& CIn = mesh_.C().internalField();
    const vectorField& pts = ibPoints();

    const scalarListList& gammaProc = ibProcGamma();
    const vectorListList& CProc = ibProcCentres();
    const vectorListList& CProcP = ibProcPts();

    // Go through all cellCells and calculate inverse distance for
    // all live points
    forAll (samplingPoints, cellI)
    {
        const vector& curP = samplingPoints[cellI];

        scalar sumW = 0;
    // for Cell part
        // Local weights
        scalarList& curCW = cellWeights[cellI];

        const labelList& curCells = ibcc[cellI];

        forAll (curCells, ccI)
        {
            // Only pick live cells
            if (gammaIn[curCells[ccI]] > SMALL)
            {
                curCW[ccI] = 1/mag(CIn[curCells[ccI]] - curP);

                sumW += curCW[ccI];
            }
            else
            {
                curCW[ccI] = 0;
            }
        }

        // Processor weights
        const List<labelPair>& interpProcCells = ibcProcC[cellI];

        scalarList& curProcCW = cellProcWeights[cellI];

        forAll (interpProcCells, cProcI)
        {
            if
            (
                gammaProc
                [
                    interpProcCells[cProcI].first()
                ]
                [
                    interpProcCells[cProcI].second()
                ] > SMALL
            )
            {
                curProcCW[cProcI] =
                    1/mag
                    (
                        CProc
                        [
                            interpProcCells[cProcI].first()
                        ]
                        [
                            interpProcCells[cProcI].second()
                        ] - curP
                    );

                sumW += curProcCW[cProcI];
            }
            else
            {
                curProcCW[cProcI] = 0;
            }
        }


    // For IB points part

        // Local weights
        scalarList& curPW = ptsWeights[cellI];

        const labelList& curPts = ibcp[cellI];

        forAll (curPts, ccI)
        {

            curPW[ccI] = 1.0/mag(pts[curPts[ccI]] - curP);

            sumW += curPW[ccI];
        }

        // Processor weights
        const List<labelPair>& interpProcPtss = ibcProcP[cellI];

        scalarList& curProcPW = ptsProcWeights[cellI];

        forAll (interpProcPtss, cProcI)
        {
 
            curProcPW[cProcI] =
                1/mag
                (
                   CProcP
                   [
                        interpProcPtss[cProcI].first()
                   ]
                   [
                       interpProcPtss[cProcI].second()
                   ] - curP
                );

             sumW += curProcPW[cProcI];
        }

        // Divide through by the sum
        if (sumW < SMALL)
        {
            InfoIn
            (
                "void immersedBoundaryFvPatch::makeIbSamplingWeights()"
            )   << "Insufficient live neighbourhood for IB cell "
                << ibc[cellI] << "." << nl
                << "Please adjust radiusFactor, angleFactor or "
                << "immersedBoundaryMaxCellCellRows "
                << "in immersedBoundaryFvPatch."
                << endl;

            // Reset sum and weights and use all points
             sumW = 0;
             curCW = 0;

             forAll (curCells, ccI)
             {
                 // Use all cells
                 curCW[ccI] = 1/mag(CIn[curCells[ccI]] - curP);
                 sumW += curCW[ccI];
             }

             forAll (curPts, ccI)
             {
                 curPW[ccI] = 1/mag(pts[curPts[ccI]] - curP);

                 sumW += curPW[ccI];
             }
        }

        forAll (curCW, ccI)
        {
            curCW[ccI] /= sumW;
        }

        forAll (curProcCW, cProcI)
        {
            curProcCW[cProcI] /= sumW;
        }

        forAll (curPW, ccI)
        {
            curPW[ccI] /= sumW;
        }

        forAll (curProcPW, cProcI)
        {
            curProcPW[cProcI] /= sumW;
        }
    }
}


const Foam::scalarListList&
Foam::immersedBoundaryFvPatch::ibNewSamplingWeights() const
{
    if (!ibNewSamplingWeightsPtr_)
    {
        makeIbNewSamplingWeights();
    }

    return *ibNewSamplingWeightsPtr_;
}


const Foam::scalarListList&
Foam::immersedBoundaryFvPatch::ibNewSamplingProcWeights() const
{
    if (!ibNewSamplingProcWeightsPtr_)
    {
        makeIbNewSamplingWeights();
    }

    return *ibNewSamplingProcWeightsPtr_;
}

const Foam::scalarListList&
Foam::immersedBoundaryFvPatch::ibNewPtsWeights() const
{
    if (!ibNewPtsWeightsPtr_)
    {
        makeIbNewSamplingWeights();
    }

    return *ibNewPtsWeightsPtr_;
}


const Foam::scalarListList&
Foam::immersedBoundaryFvPatch::ibNewPtsProcWeights() const
{
    if (!ibNewPtsProcWeightsPtr_)
    {
        makeIbNewSamplingWeights();
    }

    return *ibNewPtsProcWeightsPtr_;
}


// ************************************************************************* //
