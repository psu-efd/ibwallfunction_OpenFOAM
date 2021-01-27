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
#include "processorFvPatch.H"
#include "processorFvPatchFields.H"
//#include "foamTime.H"
#include "Time.H"
#include "fvBoundaryMesh.H"
#include "fvMesh.H"
#include "SortableList.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "treeBoundBox.H"
#include "treeDataCell.H"
#include "addToRunTimeSelectionTable.H"

#include "mathematicalConstants.H"
#include "meshSearch.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(immersedBoundaryFvPatch, 0);

    addToRunTimeSelectionTable(fvPatch, immersedBoundaryFvPatch, polyPatch);
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
/*
const Foam::scalar Foam::immersedBoundaryFvPatch::radiusFactor_
(
    debug::optimisationSwitch
    (
        "immersedBoundaryRadiusFactor",
        3.5
//         5
    )
);

const Foam::scalar Foam::immersedBoundaryFvPatch::angleFactor_
(
    debug::optimisationSwitch
    (
        "immersedBoundaryAngleFactor",
        80
//         170
    )
);
const Foam::scalar Foam::immersedBoundaryFvPatch::maxCellCellRows_
(
    debug::optimisationSwitch
    (
        "immersedBoundaryMaxCellCellRows",
        4
//         5
    )
);
const Foam::scalar Foam::immersedBoundaryFvPatch::distFactor_
(
    debug::optimisationSwitch
    (
        "immersedBoundaryDistFactor",
        1.5
    )
);
 
 */


 


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void Foam::immersedBoundaryFvPatch::makeConstant()
{
    IOdictionary ibmDict
    (
        IOobject
        (
            "ibmDict",
            mesh_.time().system(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    angleFactor_ = ibmDict.lookupOrDefault<scalar>("angleFactor",angleFactor_);
    radiusFactor_ = ibmDict.lookupOrDefault<scalar>("radiusFactor",radiusFactor_);
    maxCellCellRows_ = ibmDict.lookupOrDefault<scalar>("maxCellCellRows",maxCellCellRows_);
    distFactor_ = ibmDict.lookupOrDefault<scalar>("distFactor",distFactor_);
}

void Foam::immersedBoundaryFvPatch::makeGamma() const
{
    if (debug)
    {
        InfoIn("void immersedBoundaryFvPatch::makeGamma() const")
            << "creating fluid cells indicator "
            << "for immersed boundary" << name()
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (gammaPtr_)
    {
        FatalErrorIn("immersedBoundaryFvPatch::makeGamma() const")
            << "fluid cells indicator already exist"
            << "for immersed boundary" << name()
            << abort(FatalError);
    }

    gammaPtr_ =
        new volScalarField
        (
            IOobject
            (
                "ibGamma" + name(),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("1", dimless, 1)
        );

    scalarField& gammaI = gammaPtr_->internalField();

    // Start from gammaExt and mark IB cells and inactive
    gammaI = gammaExt().internalField();

    const labelList& ibc = ibCells();

    forAll (ibc, cellI)
    {
        gammaI[ibc[cellI]] = 0.0;
    }

    // Not allowed to call correctBoundaryConditions.  HJ, 16/Apr/2012
    // Evaluate coupled boundaries and copy out the uncoupled ones
    //gammaPtr_->boundaryField().evaluateCoupled();

    // Evaluate the coupled patchField, to replace evaluateCoupled()
    // Y.C. Xu 7/5/2017
    forAll (gammaPtr_->boundaryField(), patchI)
    {
        if (gammaPtr_->boundaryField()[patchI].coupled())
        {
            gammaPtr_->boundaryField()[patchI].initEvaluate
                (   
                    Pstream::blocking
                );
        }

        if (gammaPtr_->boundaryField()[patchI].coupled())
        {
            gammaPtr_->boundaryField()[patchI].evaluate
                (   
                    Pstream::blocking
                );
        }
    }


    forAll (gammaPtr_->boundaryField(), patchI)
    {
        if (!gammaPtr_->boundaryField()[patchI].coupled())
        {
            gammaPtr_->boundaryField()[patchI] =
                gammaPtr_->boundaryField()[patchI].patchInternalField();
        }
    }


}


void Foam::immersedBoundaryFvPatch::makeGammaExt() const
{
    if (debug)
    {
        InfoIn("void immersedBoundaryFvPatch::makeGammaExt() const")
            << "creating extended fluid cells indicator "
            << "for immersed boundary " << name()
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (gammaExtPtr_)
    {
        FatalErrorIn("immersedBoundaryFvPatch::makeGammaExt() const")
            << "extended fluid cells indicator already exist "
            << "for immersed boundary " << name()
            << abort(FatalError);
    }

    gammaExtPtr_ =
        new volScalarField
        (
            IOobject
            (
                "ibGammaExt"  + name(),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("1", dimless, 1)
        );

    IOdictionary ibmDict
        (
            IOobject
            (
                "ibmDict",
                mesh_.time().system(),
                mesh_,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        );

    pointField outsidePoints_ = ibmDict.lookup("outsidePoints");
    scalarField& gammaExtI = gammaExtPtr_->internalField();

    //const vectorField& C = mesh_.cellCentres();

    // Mark cells that are inside or outside of the triangular surface
    //boolList inside = ibPolyPatch_.triSurfSearch().calcInside(C);
    Info<<"start markInside"<<endl;

//weaken the enclosing surface condition so as to import bathymetry surface
//0-outside 1-inside 2-cut(with center inside) 3-cut(with center outside)

    labelList inside = ibPolyPatch_.markInside(mesh_,outsidePoints_);

    Info<<"finish markInside"<<endl;


    // Adjust selection of cells: inside or outside of immersed boundary
    if (internalFlow())
    {
        Info<< "Internal flow" << endl;
        forAll (gammaExtI, cellI)
        {
            if (inside[cellI]==1 or inside[cellI]==2)
            {
                gammaExtI[cellI] = 1;
            }
            if (inside[cellI]==0 or inside[cellI]==3)
            {
                gammaExtI[cellI] = 0;
            }            
        }
    }
    else
    {
        Info<< "External flow" << endl;
        forAll (gammaExtI, cellI)
        {
            if (inside[cellI]==1 or inside[cellI]==2)
            {
                gammaExtI[cellI] = 0;
            }
            if (inside[cellI]==0 or inside[cellI]==3)
            {
                gammaExtI[cellI] = 1;
            } 
        }
    }

    // Not allowed to call correctBoundaryConditions.  HJ, 16/Apr/2012
    // Evaluate coupled boundaries and copy out the uncoupled ones
    //gammaExtPtr_->boundaryField().evaluateCoupled();

    // Evaluate the coupled patchField, to replace evaluateCoupled()
    // Y.C. Xu 7/5/2017
    forAll (gammaExtPtr_->boundaryField(), patchI)
    {
        if (gammaExtPtr_->boundaryField()[patchI].coupled())
        {
            gammaExtPtr_->boundaryField()[patchI].initEvaluate
                (   
                    Pstream::blocking
                );
        }

        if (gammaExtPtr_->boundaryField()[patchI].coupled())
        {
            gammaExtPtr_->boundaryField()[patchI].evaluate
                (   
                    Pstream::blocking
                );
        }
    }


    forAll (gammaExtPtr_->boundaryField(), patchI)
    {
        if (!gammaExtPtr_->boundaryField()[patchI].coupled())
        {
            gammaExtPtr_->boundaryField()[patchI] =
                gammaExtPtr_->boundaryField()[patchI].patchInternalField();
        }
    } 
}


void Foam::immersedBoundaryFvPatch::makeSGamma() const
{
    if (debug)
    {
        InfoIn("void immersedBoundaryFvPatch::makeSGamma() const")
            << "creating fluid faces indicator "
            << "for immersed boundary " << name()
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (sGammaPtr_)
    {
        FatalErrorIn("immersedBoundaryFvPatch::makeSGamma() const")
            << "fluid faces indicator already exist "
            << "for immersed boundary " << name()
            << abort(FatalError);
    }

    // Note: change in algorithm
    // 1) First, mark all faces as dead
    // 2) Mark faces as live if they are between two live cells
    sGammaPtr_ =
        new surfaceScalarField
        (
            IOobject
            (
                "sGamma"  + name(),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0),
            calculatedFvsPatchScalarField::typeName
        );

    // Get access to components of sGamma
    scalarField& sGammaI = sGammaPtr_->internalField();

    surfaceScalarField::GeometricBoundaryField& sGammaPatches =
        sGammaPtr_->boundaryField();

    const unallocLabelList& owner = mesh_.owner();
    const unallocLabelList& neighbour = mesh_.neighbour();

    // Live cells indicator
    const volScalarField& g = gamma();

    // Extended live cells indicator
    const volScalarField& gExt = gammaExt();
    const scalarField& gExtIn = gExt.internalField();

    const volScalarField::GeometricBoundaryField& gExtPatches =
        gExt.boundaryField();

    volScalarField gIb = gExt - g;
    const scalarField& gIbIn = gIb.internalField();

    const volScalarField::GeometricBoundaryField& gIbPatches =
        gIb.boundaryField();

    // Internal faces: flux is live between all active and IB cells
    forAll (sGammaI, faceI)
    { 
        // If both cells are live, flux is live
        if
        (
            gExtIn[owner[faceI]] > SMALL
         && gExtIn[neighbour[faceI]] > SMALL
        )
        {
            sGammaI[faceI] = 1;
        }
    }

    // Kill fluxes between two IB cells
    forAll (sGammaI, faceI)
    {
        if
        (
            gIbIn[owner[faceI]] > SMALL
         && gIbIn[neighbour[faceI]] > SMALL
        )
        {
            sGammaI[faceI] = 0;
        }
    }

    forAll (gExtPatches, patchI)
    {
        if (gExtPatches[patchI].coupled())
        {
            scalarField& gP = sGammaPatches[patchI];

            // For coupled patches, check gammaExt
            scalarField gammaOwn = gExtPatches[patchI].patchInternalField();

            scalarField gammaNei = gExtPatches[patchI].patchNeighbourField();

            forAll (gammaOwn, faceI)
            {
                if
                (
                    gammaOwn[faceI] > SMALL
                 && gammaNei[faceI] > SMALL
                )
                {
                    gP[faceI] = 1;
                }
            }

            // For coupled patches, kill IB
            scalarField gammaIbOwn = gIbPatches[patchI].patchInternalField();

            scalarField gammaIbNei = gIbPatches[patchI].patchNeighbourField();

            forAll (gammaOwn, faceI)
            {
                if
                (
                    gammaIbOwn[faceI] > SMALL
                 && gammaIbNei[faceI] > SMALL
                )
                {
                    gP[faceI] = 0;
                }
            }
        }
        else
        {
            // For regular patches, check live cells only to achieve
            // correct global mass adjustment.
            // HJ, 21/May/2012
            scalarField gammaFc =
                g.boundaryField()[patchI].patchInternalField();

            scalarField& gP = sGammaPatches[patchI];

            forAll (gammaFc, faceI)
            {
                if (gammaFc[faceI] > SMALL)
                {
                   gP[faceI] = 1;
                }
            }
        }
    }
}


void Foam::immersedBoundaryFvPatch::makeIbCells() const
{
    if (debug)
    {
        InfoIn("void immersedBoundaryFvPatch::makeIbCells() const")
            << "create list of cells next to immersed boundary "
            << "for immersed boundary " << name()
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (ibCellsPtr_)
    {
        FatalErrorIn("immersedBoundaryFvPatch::makeIbCells() const")
            << "list of cells next to immersed boundary already exist "
            << "for immersed boundary " << name()
            << abort(FatalError);
    }

    labelHashSet ibCellSet;

    const unallocLabelList& owner = mesh_.owner();
    const unallocLabelList& neighbour = mesh_.neighbour();

    const volScalarField gE = gammaExt();
    const scalarField& gammaExtI = gE.internalField();

    forAll (neighbour, faceI)
    {
        if (mag(gammaExtI[neighbour[faceI]] - gammaExtI[owner[faceI]]) > SMALL)
        {
            if (gammaExtI[owner[faceI]] > SMALL)
            {
                if (!ibCellSet.found(owner[faceI]))
                {
                    ibCellSet.insert(owner[faceI]);
                }
            }
            else
            {
                if (!ibCellSet.found(neighbour[faceI]))
                {
                    ibCellSet.insert(neighbour[faceI]);
                }
            }
        }
    }

    // check cells next to the processor interface
    forAll (gE.boundaryField(), patchI)
    {
        if (gE.boundaryField()[patchI].coupled())
        {
            scalarField gammaExtOwn =
                gE.boundaryField()[patchI].patchInternalField();

            scalarField gammaExtNei =
                gE.boundaryField()[patchI].patchNeighbourField();

            const unallocLabelList& fCells =
                mesh_.boundary()[patchI].faceCells();

            forAll (gammaExtOwn, faceI)
            {
                if
                (
                    mag(gammaExtNei[faceI] - gammaExtOwn[faceI])
                  > SMALL
                )
                {
                    if (gammaExtOwn[faceI] > SMALL)
                    {
                        if (!ibCellSet.found(fCells[faceI]))
                        {
                            ibCellSet.insert(fCells[faceI]);
                        }
                    }
                    else if (2*gammaExtOwn.size() == fCells.size())
                    {
                        if
                        (
                           !ibCellSet.found
                            (
                                fCells[gammaExtOwn.size() + faceI]
                            )
                        )
                        {
                            ibCellSet.insert
                            (
                                fCells[gammaExtOwn.size() + faceI]
                            );
                        }
                    }
                }
            }
        }
    }


    ibCellsPtr_ = new labelList(ibCellSet.toc());

    //addIbCornerCells();

    //const labelList ibc = ibCellSet.toc();

    //addAdjacentCells(ibc);

    sort(*ibCellsPtr_);

    Info << "Number of IB cells: " << returnReduce(ibCellsPtr_->size(), sumOp<label>())<< endl;
    //Pout << "Number of IB cells: " << ibCellsPtr_->size() << endl;
}

#include "immersedBoundaryAddAdjacentCells.C"

void Foam::immersedBoundaryFvPatch::addIbCornerCells() const
{
    if (debug)
    {
        InfoIn("void immersedBoundaryFvPatch::addIbCornerCells() const")
            << "add cells next to sharp corner "
            << "for immersed boundary " << name()
            << endl;
    }

    const vectorField& C = mesh_.cellCentres();
    const vectorField& Cf = mesh_.faceCentres();
    const vectorField& Sf = mesh_.faceAreas();

    const cellList& meshCells = mesh_.cells();
    const unallocLabelList& own = mesh_.owner();
    const unallocLabelList& nei = mesh_.neighbour();

    label nCornerCells = 0;
    labelList cornerCells;

    const triSurfaceSearch& tss = ibPolyPatch_.triSurfSearch();



    do
    {
        // Access to derived data from do-loop: this will be re-calculated
        // HJ, 21/May/2012
        const labelList& ibc = ibCells();

        labelHashSet ibCellSet;

        forAll (ibc, cellI)
        {
            ibCellSet.insert(ibc[cellI]);
        }

        // Note: the algorithm is originally written with inward-facing normals
        // and subsequently changed: IB surface normals point outwards
        // HJ, 21/May/2012
        const vectorField& ibn = ibNormals();

        const scalarField& gammaI = gamma().internalField();

        labelHashSet cornerIbCellSet;

        const labelList& ibf = ibFaces();

        forAll (ibf, faceI)
        {
            if( ibf[faceI] > own.size() )
            {
                break;
            }

            const label& ownCell = own[ibf[faceI]];
            const label& neiCell = nei[ibf[faceI]];

//             label ibCell = -1;
            label liveCell = -1;

            if (gammaI[ownCell] > SMALL)
            {
//                 ibCell = neiCell;
                liveCell = ownCell;
            }
            else
            {
//                 ibCell = ownCell;
                liveCell = neiCell;
            }

            scalar delta = cellSize(liveCell);
            vector span(2*delta, 2*delta, 2*delta);

            pointIndexHit pih = tss.nearest(C[liveCell], span);

            if (pih.hit())
            {
                vector n =
                    triSurfaceTools::surfaceNormal
                    (
                        ibPolyPatch_.ibMesh(),
                        pih.index(),
                        pih.hitPoint()
                    );

                scalar totalArea = 0;
                {
                    const labelList& cellFaces = meshCells[liveCell];

                    forAll (cellFaces, faceI)
                    {
                        label curFace = cellFaces[faceI];

                        vector curSf = Sf[curFace];

                        if ((curSf & (Cf[curFace] - C[liveCell])) < 0)
                        {
                            curSf *= -1;
                        }

                        if ((curSf & n) > 0)
                        {
                            totalArea += (curSf & n);
                        }
                    }
                }

                scalar area = 0;
                {
                    const labelList& cellFaces = meshCells[liveCell];

                    forAll (cellFaces, faceI)
                    {
                        label curFace = cellFaces[faceI];

                        vector curSf = Sf[curFace];

                        label neiCell = -1;

                        //HJ bug fix?
                        if (mesh_.isInternalFace(curFace))
                        {
                            if (own[curFace] == liveCell)
                            {
                                neiCell = nei[curFace];
                            }
                            else
                            {
                                neiCell = own[curFace];
                            }
                        }

                        label ibCell = findIndex(ibc, neiCell);

                        if (ibCell != -1)
                        {
                            // Note that ibn points outwards
                            if ((-ibn[ibCell] & n) > 0)
                            {
                                area += mag(curSf & n);
                            }
                        }
                    }
                }

                if (area/totalArea < 0.5)
                {
                    if (!cornerIbCellSet.found(liveCell) and !ibCellSet.found(liveCell))
                    {
                        cornerIbCellSet.insert(liveCell);
                    }
                }
            }
        }

        cornerCells = cornerIbCellSet.toc();
        ibCellsPtr_->append(cornerCells);
        nCornerCells += cornerCells.size();

        //deleteDemandDrivenData(gammaPtr_);

        deleteDemandDrivenData(ibFacesPtr_);
        deleteDemandDrivenData(ibFaceCellsPtr_);
        deleteDemandDrivenData(ibFaceFlipsPtr_);

        deleteDemandDrivenData(ibPointsPtr_);
        deleteDemandDrivenData(ibNormalsPtr_);
        deleteDemandDrivenData(hitFacesPtr_);
        deleteDemandDrivenData(ibSamplingPointsPtr_);

        deleteDemandDrivenData(ibSamplingWeightsPtr_);
        deleteDemandDrivenData(ibSamplingProcWeightsPtr_);
    }
    while (cornerCells.size() > 0);
}


void Foam::immersedBoundaryFvPatch::makeIbFaces() const
{
    if (debug)
    {
        InfoIn("void immersedBoundaryFvPatch::makeIbFaces() const")
            << "create list of faces next to immersed boundary "
            << "for immersed boundary " << name()
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (ibFacesPtr_ || ibFaceCellsPtr_ || ibFaceFlipsPtr_)
    {
        FatalErrorIn("immersedBoundaryFvPatch::makeIbFaces() const")
            << "list of faces next to immersed boundary already exist "
            << "for immersed boundary " << name()
            << abort(FatalError);
    }

    // Mark IB cells with their index
    const labelList& ibc = ibCells();

    labelList ibCellIndicator(mesh_.nCells(), -1);

    forAll (ibc, ibcI)
    {
        ibCellIndicator[ibc[ibcI]] = ibcI;
    }

    dynamicLabelList ibF(6*ibc.size());
    dynamicLabelList ibFC(6*ibc.size());
    DynamicList<bool> ibFF(6*ibc.size());

    const unallocLabelList& owner = mesh_.owner();
    const unallocLabelList& neighbour = mesh_.neighbour();

    const scalarField& gammaI = gamma().internalField();
 
    forAll (neighbour, faceI)
    {
        // Old check done using mag gamma.  Changed for internal cells
        if
        (
//            ibCellIndicator[owner[faceI]] == -1
//         && ibCellIndicator[neighbour[faceI]] > -1
            gammaI[owner[faceI]] > SMALL
         && gammaI[neighbour[faceI]] < SMALL
        )
        {
            // Owner is live, neighbour IB.  Its IB index is in
            // ibCellIndicator
            ibF.append(faceI);
            ibFC.append(ibCellIndicator[neighbour[faceI]]);
            ibFF.append(false);
        }
        else if
        (
//            ibCellIndicator[owner[faceI]] > -1
//         && ibCellIndicator[neighbour[faceI]] == -1
             gammaI[owner[faceI]] < SMALL
          && gammaI[neighbour[faceI]] > SMALL
        )
        {
            // Neighbour is live, owner IB.  Its IB index is in
            // ibCellIndicator
            ibF.append(faceI);
            ibFC.append(ibCellIndicator[owner[faceI]]);
            ibFF.append(true);
        }
    }

    volScalarField::GeometricBoundaryField& gammaPatches =
        gammaPtr_->boundaryField();
/*
    tmp<volScalarField > tgamma
    (
        new volScalarField(gamma())
    );

    volScalarField& gamma = tgamma();
	
	volScalarField::GeometricBoundaryField& gammaPatches =
        gamma.boundaryField();*/

    forAll (gammaPatches, patchI)
    {
        // Note: take faceCells from fvPatch (because of empty)
        const labelList& fc = mesh_.boundary()[patchI].faceCells();
        const label start = mesh_.boundaryMesh()[patchI].start();

        if (gammaPatches[patchI].coupled())
        {

            gammaPatches[patchI].initEvaluate(Pstream::blocking);// not quite sure
            gammaPatches[patchI].evaluate(Pstream::blocking);// not quite sure

            scalarField gammaOwn =
                gammaPatches[patchI].patchInternalField();

            scalarField gammaNei =
                gammaPatches[patchI].patchNeighbourField();

            forAll (gammaOwn, patchFaceI)
            {
                if
                (
                    mag(gammaNei[patchFaceI] - gammaOwn[patchFaceI]) > SMALL
                )
                {
                    if (ibCellIndicator[fc[patchFaceI]] > -1)
                    {  
                        // Owner cell is IB
                        ibF.append(start + patchFaceI);
                        ibFC.append(ibCellIndicator[fc[patchFaceI]]);
                        ibFF.append(true);
                    }
                    else
                    {  
                        // Neighbour cell is IB
                        ibF.append(start + patchFaceI);
                        ibFC.append(-1);
                        ibFF.append(false);
                    }
                }
            }
        }
    }

    // Pack the data
    ibF.shrink();
    ibFC.shrink();
    ibFF.shrink();

    ibFacesPtr_ = new labelList(ibF.xfer());

    ibFaceCellsPtr_ = new labelList(ibFC.xfer());
    ibFaceFlipsPtr_ = new boolList(ibFF.xfer());
}


void Foam::immersedBoundaryFvPatch::makeIbInsideFaces() const
{
    if (debug)
    {
        InfoIn("void immersedBoundaryFvPatch::makeIbInsideFaces() const")
            << "create list of faces next to immersed boundary "
            << "for immersed boundary " << name()
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (ibInsideFacesPtr_)
    {
        FatalErrorIn("immersedBoundaryFvPatch::makeIbInsideFaces() const")
            << "list of faces next to immersed boundary already exist "
            << "for immersed boundary " << name()
            << abort(FatalError);
    }

    labelHashSet ibInsideFSet;

    const labelList& owner = mesh_.faceOwner();
    const labelList& neighbour = mesh_.faceNeighbour();

    const volScalarField gE = gammaExt();
    const scalarField& gammaExtI = gE.internalField();

    forAll (neighbour, faceI)
    {
        if (mag(gammaExtI[neighbour[faceI]] - gammaExtI[owner[faceI]]) > SMALL)
        {
            ibInsideFSet.insert(faceI);
        }
    }

    forAll (gE.boundaryField(), patchI)
    {
        if (gE.boundaryField()[patchI].coupled())
        {
            scalarField gammaOwn =
                gE.boundaryField()[patchI].patchInternalField();

            scalarField gammaNei =
                gE.boundaryField()[patchI].patchNeighbourField();

            label size = mesh_.boundaryMesh()[patchI].size();
            label start = mesh_.boundaryMesh()[patchI].start();

            forAll (gammaOwn, faceI)
            {
                if
                (
                    mag(gammaNei[faceI] - gammaOwn[faceI]) > SMALL
                )
                {
                    if (!ibInsideFSet.found(start + faceI))
                    {
                        ibInsideFSet.insert(start + faceI);
                    }

                    if (2*gammaOwn.size() == size)
                    {
                        if
                        (
                           !ibInsideFSet.found
                            (
                                start + size/2 + faceI
                            )
                        )
                        {
                            ibInsideFSet.insert
                            (
                                start + size/2 + faceI
                            );
                        }
                    }
                }
            }
        }
    }

    ibInsideFacesPtr_ = new labelList(ibInsideFSet.toc());
    sort(*ibInsideFacesPtr_);
}


void Foam::immersedBoundaryFvPatch::makeIbInternalFaces() const
{
    if (debug)
    {
        InfoIn("void immersedBoundaryFvPatch::makeIbInternalFaces() const")
            << "create list of faces next to immersed boundary"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (ibInternalFacesPtr_)
    {
        FatalErrorIn("immersedBoundaryFvPatch::makeIbInternalFaces() const")
            << "list of faces next to immersed boundary already exist"
            << abort(FatalError);
    }

    labelHashSet ibInternalFacesSet;

    const labelList& owner = mesh_.faceOwner();
    const labelList& neighbour = mesh_.faceNeighbour();

    volScalarField gammaTmp
    (
        "gammaTmp",
        gammaExt() - gamma()
    );
    const scalarField& gammaTmpI = gammaTmp.internalField();

    forAll (neighbour, faceI)
    {
        if
        (
            (gammaTmpI[neighbour[faceI]] > SMALL)
         && (gammaTmpI[owner[faceI]] > SMALL)
        )
        {
            ibInternalFacesSet.insert(faceI);
        }
    }

    forAll (gammaTmp.boundaryField(), patchI)
    {
        if (gammaTmp.boundaryField()[patchI].coupled())
        {
            scalarField gammaOwn =
                gammaTmp.boundaryField()[patchI].patchInternalField();

            scalarField gammaNei =
                gammaTmp.boundaryField()[patchI].patchNeighbourField();

            label size = mesh_.boundaryMesh()[patchI].size();
            label start = mesh_.boundaryMesh()[patchI].start();

            forAll (gammaOwn, faceI)
            {
                if
                (
                    (gammaNei[faceI] > SMALL)
                 && (gammaOwn[faceI] > SMALL)
                )
                {
                    if (!ibInternalFacesSet.found(start + faceI))
                    {
                        ibInternalFacesSet.insert(start + faceI);
                    }

                    if (2*gammaOwn.size() == size)
                    {
                        if
                        (
                           !ibInternalFacesSet.found
                            (
                                start + size/2 + faceI
                            )
                        )
                        {
                            ibInternalFacesSet.insert
                            (
                                start + size/2 + faceI
                            );
                        }
                    }
                }
            }
        }
    }

    ibInternalFacesPtr_ = new labelList(ibInternalFacesSet.toc());
    sort(*ibInternalFacesPtr_);
}


void Foam::immersedBoundaryFvPatch::makeIbPointsAndNormals() const
{
 


    if (debug)
    {
        InfoIn
        (
            "void immersedBoundaryFvPatch::makeIbPointsAndNormals() const"
        )   << "create immersed  boundary points and normals "
            << "for immersed boundary " << name()
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (ibPointsPtr_ || ibNormalsPtr_ || hitFacesPtr_ || ibSamplingPointsPtr_)
    {
        FatalErrorIn
        (
            "immersedBoundaryFvPatch::makeIbPointsAndNormals() const"
        )
            << "immersed boundary points and normals already exist"
            << "for immersed boundary " << name()
            << abort(FatalError);
    }

    // Find average cell dimension
    const labelList& ibc = ibCells();

    scalarField delta(ibc.size(), 0.0);

    forAll (delta, cellI)
    {
        delta[cellI] = cellSize(ibc[cellI]);
    }

    // Find nearest triSurface point for each interface cell centre
    ibPointsPtr_ = new vectorField(ibc.size(), vector::zero);
    ibNormalsPtr_ = new vectorField(ibc.size(), vector::zero);
    hitFacesPtr_ = new labelList(ibc.size(), -1);
    ibSamplingPointsPtr_ = new vectorField(ibc.size(), vector::zero);

    vectorField& ibPoints = *ibPointsPtr_;
    vectorField& ibNormals = *ibNormalsPtr_;
    labelList& ibHitFaces = *hitFacesPtr_;
    vectorField& ibSamplingPoints = *ibSamplingPointsPtr_;

    const vectorField& C = mesh_.C().internalField();

    // Get IB cell centres
    vectorField ibCellCentres(C, ibc);

    const triSurfaceSearch& tss = ibPolyPatch_.triSurfSearch();

    forAll (ibc, cellI)
    {
        // Adjust search span if needed.  HJ, 14/Dec/2012
        vector span
        (
            2*radiusFactor_*delta[cellI],
            2*radiusFactor_*delta[cellI],
            2*radiusFactor_*delta[cellI]
        );

        pointIndexHit pih = tss.nearest(ibCellCentres[cellI], span);

        if (pih.hit())
        {
            ibPoints[cellI] = pih.hitPoint();
            ibNormals[cellI] =
                triSurfaceTools::surfaceNormal
                (
                    ibPolyPatch_.ibMesh(),
                    pih.index(),
                    pih.hitPoint()
                );
            scalar indicator = ibNormals[cellI] & (ibCellCentres[cellI] - ibPoints[cellI]);
            indicator = indicator/mag(indicator);
            // Note: ibNormals point OUT of the domain
            if (!internalFlow() and indicator < 0)
            {
                ibNormals[cellI] *= -1;
            }

            ibHitFaces[cellI] = pih.index();
        }
        else
        {
            FatalErrorIn
            (
                "immersedBoundaryFvPatch::makeIbPointsAndNormals() const"
            )   << "Can't find nearest triSurface point for cell "
                << ibc[cellI] << ", "
                << mesh_.cellCentres()[ibc[cellI]]
                << ".  Hit data = " << pih << nl
                << abort(FatalError);
        }

        if
        (
            mesh_.nGeometricD() < 3
         && mag(ibCellCentres[cellI].z() - ibPoints[cellI].z()) > SMALL
        )
        {
            WarningIn
            (
                "immersedBoundaryFvPatch::makeIbPointsAndNormals() const"
            )   << "Intersection point is not on symmetry plane " << nl
                << "C = " << ibCellCentres[cellI]
                <<  " D = " <<  ibPoints[cellI] << nl
                << "for 2-D geometry.  Adjusting" << endl;

               ibPoints[cellI].z() = ibCellCentres[cellI].z();
        }
    }

    // Calculate sampling points locations


    //ibSamplingPoints = ibPoints + distFactor_*(ibCellCentres - ibPoints);
    scalarField rM(ibCellSizes());
    ibSamplingPoints = ibPoints + distFactor_*rM*ibNormals;
}


void Foam::immersedBoundaryFvPatch::makeIbCellCells() const
{



    if (debug)
    {
        InfoIn("void immersedBoundaryFvPatch::makeIbCellCells() const")
            << "create neighbour cells for ib cells"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if
    (
        ibCellCellsPtr_
     || ibProcCentresPtr_
     || ibProcGammaPtr_
     || ibCellProcCellsPtr_
    )
    {
        FatalErrorIn
        (
            "immersedBoundaryFvPatch::makeIbCellCells() const"
        )   << "cell-cell addressing already exists"
            << abort(FatalError);
    }

    const labelList& ibc = ibCells();

    ibCellCellsPtr_ = new labelListList(ibc.size());
    labelListList& cellCells = *ibCellCellsPtr_;

    ibProcCentresPtr_ = new vectorListList(Pstream::nProcs());
    vectorListList& procCentres = *ibProcCentresPtr_;

    ibProcGammaPtr_ = new scalarListList(Pstream::nProcs());
    scalarListList& procGamma = *ibProcGammaPtr_;

    ibCellProcCellsPtr_ = new List<List<labelPair> >(ibc.size());
    List<List<labelPair> >& cellProcCells = *ibCellProcCellsPtr_;

    const cellList& meshCells = mesh_.cells();

    // In geometry initialisation, fields are not available: use raw mesh data
    // HJ after ZT, 6/Dec/2012
    const vectorField& C = mesh_.cellCentres();

    scalarField rM(ibCellSizes());

    rM *= radiusFactor_;

    const vectorField& ibp = ibPoints();

    // Note: the algorithm is originally written with inward-facing normals
    // and subsequently changed: IB surface normals point outwards
    // HJ, 21/May/2012
    const vectorField& ibn = ibNormals();

    forAll (cellCells, cellI)
    {
        labelList curCells;

        findCellCells
        (
            C[ibc[cellI]],
            ibc[cellI],
            curCells
        ); 
 
        cellCells[cellI] = labelList(curCells.size(), -1);

        label cI = 0;

        for (label i = 0; i < curCells.size(); i++)
        {
            label curCell = curCells[i];
            scalar r = mag(C[curCell] - C[ibc[cellI]]);
 
            if (r <= rM[cellI])
            {
                scalar angleLimit =
                    -Foam::cos(angleFactor_*constant::mathematical::pi/180);

                vector dir = (C[curCell] - ibp[cellI]);
                dir /= mag(dir) + SMALL;

                if ((-ibn[cellI] & dir) >= angleLimit)
                {
                    cellCells[cellI][cI++] = curCell;
                }
            }
        }

        //cellCells[cellI].setSize(cI);
		cellCells[cellI].setSize(min(cI,45)); // to limit the number of cellCells
    }

    // Find immersed boundary cells in each processor 
    if (ibProcIbCellsPtr_)
    {
        FatalErrorIn("immersedBoundaryFvPatch::makeIbCellCells() const")
            << "procCells addressing already exists"
            << abort(FatalError);
    }

    ibProcIbCellsPtr_ = new labelListList(Pstream::nProcs());
    labelListList& procIbCells = *ibProcIbCellsPtr_;

    // Find cells needed by other processors
    if (ibProcCellsPtr_)
    {
        FatalErrorIn("immersedBoundaryFvPatch::makeIbCellCells() const")
            << "procCells addressing already exists"
            << abort(FatalError);
    }
    ibProcCellsPtr_ = new labelListList(Pstream::nProcs());
    labelListList& procCells = *ibProcCellsPtr_;

    if (Pstream::parRun())
    {
        // Find immersed boundary cells in each processor
        // Y.C. Xu 7/4/2017
        procIbCells[Pstream::myProcNo()] = ibc;

        Pstream::gatherList(procIbCells);
        Pstream::scatterList(procIbCells);

        // Find immersed boundary cells whose cellCells next to processor boundaries
        labelHashSet nextProcIbCellsSet;

        forAll (ibc, cellI)
        {
            const labelList& curCellCells = cellCells[cellI];

            if (curCellCells.size())
            {
                forAll (curCellCells, cI)
                {
                    const labelList& faces = meshCells[curCellCells[cI]];

                    bool foundProcessorFace = false;

                    forAll (faces, faceI)
                    {
                        label patchID =
                            mesh_.boundaryMesh().whichPatch(faces[faceI]);

                        if (patchID != -1)
                        {
                            if
                            (
                                isA<processorPolyPatch>
                                (
                                    mesh_.boundaryMesh()[patchID]
                                )
                            )
                            {
                                foundProcessorFace = true;
                            }
                        }
                    }

                    if (foundProcessorFace)
                    {
                        nextProcIbCellsSet.insert(cellI);
                        break;
                    }
                }
            }
            else
            {
                const labelList& faces = meshCells[ibc[cellI]];

                bool foundProcessorFace = false;

                forAll (faces, faceI)
                {
                    label patchID =
                        mesh_.boundaryMesh().whichPatch(faces[faceI]);

                    if (patchID != -1)
                    {
                        if
                        (
                            isA<processorPolyPatch>
                            (
                                mesh_.boundaryMesh()[patchID]
                            )
                        )
                        {
                            foundProcessorFace = true;
                        }
                    }
                }

                if (foundProcessorFace)
                {
                    nextProcIbCellsSet.insert(cellI);
                }
            }
        }
        labelList nextProcIbCells = nextProcIbCellsSet.toc();

        sort(nextProcIbCells);

        // Note: new gather-scatter operations
        // HJ, 11/Aug/2016

        // Send and receive ibc centres and radii
        vectorListList ctrs(Pstream::nProcs());

        ctrs[Pstream::myProcNo()].setSize(procIbCells[Pstream::myProcNo()].size());
        vectorList& centres = ctrs[Pstream::myProcNo()];

        centres = vectorList(C, procIbCells[Pstream::myProcNo()]);
 
        Pstream::gatherList(ctrs);
        Pstream::scatterList(ctrs);

        scalarListList rMax(Pstream::nProcs());

        rMax[Pstream::myProcNo()] = scalarField(rM);

        Pstream::gatherList(rMax);
        Pstream::scatterList(rMax); 

        // Find cells needed by other processors
        const scalarField& gammaExtI = gammaExt().internalField();

        //labelHashSet procCellSet;

        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {

            labelHashSet procCellSet;

            if (procI != Pstream::myProcNo())
            {

                forAll (ctrs[procI], cellI)
                {
                    label nearestCellID =
                        findNearestCell(ctrs[procI][cellI]);

                    if (nearestCellID == -1)
                    {
                        FatalErrorIn
                        (
                            "immersedBoundaryFvPatch::makeIbCellCells() const"
                        ) << "Can't find nearest cell."
                            << abort(FatalError);
                    }

                    scalar R = mag(C[nearestCellID] - ctrs[procI][cellI]);

                    if (R < rMax[procI][cellI])
                    {

                        if (!procCellSet.found(nearestCellID) and gammaExtI[cellI] > SMALL)
                        {
                            procCellSet.insert(nearestCellID);
                        }

                        labelList tmpCellList;

                        findCellCells
                        (
                            ctrs[procI][cellI],
                            nearestCellID,
                            tmpCellList
                        );

                        forAll (tmpCellList, cI)
                        {
                            scalar r =
                                mag
                                (
                                    C[tmpCellList[cI]]
                                  - ctrs[procI][cellI]
                                );

                            if (r <= rMax[procI][cellI])
                            {
                                if (!procCellSet.found(tmpCellList[cI]))
                                {
                                    procCellSet.insert(tmpCellList[cI]);
                                }
                            }
                        }
                    }
                }
            }
			procCells[procI] = procCellSet.toc();          
        }

        //procCells[Pstream::myProcNo()] = procCellSet.toc();

		forAll(procCells,i)
		{
			//Pout<<"procCells "<<procCells[i].size()<<" "<<i<<endl;
		}
/*
        Pstream::gatherList(procCells);
        Pstream::scatterList(procCells);
 
        // Send cell center
        procCentres[Pstream::myProcNo()] =
            vectorField
            (
                C,
                procCells[Pstream::myProcNo()]
            );

        Pstream::gatherList(procCentres);
        Pstream::scatterList(procCentres);

        // Send cell gamma
        procGamma[Pstream::myProcNo()] =
            scalarField
            (
                gamma().internalField(),
                procCells[Pstream::myProcNo()]
            );

        Pstream::gatherList(procGamma);
        Pstream::scatterList(procGamma);
*/
        // Send and receive sizes
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                // Parallel data exchange
                {
                    OPstream toProc
                    (
                        Pstream::blocking,
                        procI,
                        sizeof(label)
                    );

                    toProc << procCells[procI].size();
                }
            }
        }

        labelList procSizes(Pstream::nProcs(), 0);
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                // Parallel data exchange
                {
                    IPstream fromProc
                    (
                        Pstream::blocking,
                        procI,
                        sizeof(label)
                    );

                    fromProc >> procSizes[procI];
                }
            }
        }

        // Send cell centres
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                vectorField centres(C, procCells[procI]);

                // Parallel data exchange
                {
                    OPstream toProc
                    (
                        Pstream::blocking,
                        procI,
                        centres.size()*sizeof(vector)
                    );

                    toProc << centres;
                }
            }
        }

        // Receive cell centres
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                procCentres[procI].setSize(procSizes[procI]);

                // Parallel data exchange
                {
                    IPstream fromProc
                    (
                        Pstream::blocking,
                        procI,
                        procSizes[procI]*sizeof(vector)
                    );

                    fromProc >> procCentres[procI];
                }
            }
            // else: already set to zero-size field
        }
		forAll(procCentres,i)
		{
			//Pout<<"procCentres "<<procCentres[i].size()<<" "<<i<<endl;
		}
        // Send cell gamma
        const scalarField& gammaI = gamma().internalField();
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                scalarField gamma(gammaI, procCells[procI]);

                // Parallel data exchange
                {
                    OPstream toProc
                    (
                        Pstream::blocking,
                        procI,
                        gamma.size()*sizeof(scalar)
                    );

                    toProc << gamma;
                }
            }
        }

        // Receive cell gamma
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                procGamma[procI].setSize(procSizes[procI]);

                // Parallel data exchange
                {
                    IPstream fromProc
                    (
                        Pstream::blocking,
                        procI,
                        procSizes[procI]*sizeof(scalar)
                    );

                    fromProc >> procGamma[procI];
                }
            }
            // else: already set to zero-size field
        }
        // Reset own size to zero?  HJ, 11/Aug/2016

        // Cell-procCells addressing
        forAll (cellProcCells, cellI)
        {
            scalar rMax = rM[cellI];

            cellProcCells[cellI].setSize(1000);
            //cellProcCells[cellI].clear();

            label index = 0;
            forAll (procCentres, procI)
            {
                if(procI!=Pstream::myProcNo())
                {
                    forAll (procCentres[procI], pointI)
                    {
                        scalar r =
                            mag
                            (
                                procCentres[procI][pointI]
                              - C[ibc[cellI]]
                            );
                        if (r <= rMax)
                        {
                            vector dir = (procCentres[procI][pointI] - ibp[cellI]);
                            dir /= mag(dir) + SMALL;
 
                            scalar angleLimit =
                               -Foam::cos(angleFactor_*constant::mathematical::pi/180);

                            // Change of sign of normal.  HJ, 21/May/2012
                            if ((ibn[cellI] & dir) >= angleLimit)
                            { 
                                cellProcCells[cellI][index].first() = procI;
                                cellProcCells[cellI][index].second() = pointI;
                                //cellProcCells[cellI].append(tmp);
                                index++;
                            }
                        }
                    }
                }
            }

            cellProcCells[cellI].setSize(index);
        }

		// to sort all cellProcCells and cellCells, find 45th nearest proCentre distance maxDistance
		scalarField maxDistance(ibc.size(),0); 

        forAll (cellProcCells, cellI)
        {
			scalarField distances(cellProcCells[cellI].size()+cellCells[cellI].size(),0);

			label totalIndex = 0;

			forAll (cellCells[cellI], index)
	        {

				scalar r =
                    mag
                    (
                        C[cellCells[cellI][index]]
                      - C[ibc[cellI]]
                    );
				distances[totalIndex] = r;
				totalIndex++;	
			}

			forAll (cellProcCells[cellI], index)
	        {
                label procI = cellProcCells[cellI][index].first();
                label pointI = cellProcCells[cellI][index].second();
				scalar r =
                    mag
                    (
                        procCentres[procI][pointI]
                      - C[ibc[cellI]]
                    );
				distances[totalIndex] = r;
				totalIndex++;
			}

			SortableList<scalar> sortedDistances(distances);
			maxDistance[cellI] = sortedDistances[min(45,distances.size())-1];

			if(distances.size()<45)
			{
				Pout<<"maxCellCellRows and radiusFactor needs to be larger"<<endl;
			}
		}

        // make sure every cellCell and cellProcCells is smaller than maxDistance
        forAll (ibc, cellI)
        {
            label newIndex = 0;

			forAll (cellCells[cellI], index)
	        {

				scalar r =
                    mag
                    (
                        C[cellCells[cellI][index]]
                      - C[ibc[cellI]]
                    );
 
            	if (r <= maxDistance[cellI])
	            {
	                newIndex++;
	            }
			}
			cellCells[cellI].setSize(newIndex);
			newIndex = 0;
			forAll (cellProcCells[cellI], index)
	        {
                label procI = cellProcCells[cellI][index].first();
                label pointI = cellProcCells[cellI][index].second();
				scalar r =
                    mag
                    (
                        procCentres[procI][pointI]
                      - C[ibc[cellI]]
                    );
	            if (r <= maxDistance[cellI])
	            {
	                newIndex++;
	            }
			}
            cellProcCells[cellI].setSize(newIndex);
        }
    }
}


void Foam::immersedBoundaryFvPatch::makeDeadCells() const
{
    if (debug)
    {
        InfoIn("void immersedBoundaryFvPatch::makeDeadCells() const")
            << "create list of dead cells"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (deadCellsPtr_)
    {
        FatalErrorIn("immersedBoundaryFvPatch::makeDeadCells() const")
            << "list of dead cells already exist"
            << abort(FatalError);
    }

    const scalarField& gammaExtI = gammaExt().internalField();

    deadCellsPtr_ = new labelList(label(sum(scalar(1) - gammaExtI)), -1);
    labelList& deadCells = *deadCellsPtr_;

    label counter = 0;

    forAll (gammaExtI, cellI)
    {
        if (gammaExtI[cellI] < SMALL)
        {
            deadCells[counter++] = cellI;
        }
    }
}


void Foam::immersedBoundaryFvPatch::makeDeadCellsExt() const
{
    if (debug)
    {
        InfoIn("void immersedBoundaryFvPatch::makeDeadCellsExt() const")
            << "create extended list of dead cells"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (deadCellsExtPtr_)
    {
        FatalErrorIn("immersedBoundaryFvPatch::makeDeadCellsExt() const")
            << "extended list of dead cells already exist"
            << abort(FatalError);
    }

    const scalarField& gammaI = gamma().internalField();

    deadCellsExtPtr_ = new labelList(label(sum(scalar(1) - gammaI)), -1);
    labelList& deadCellsExt = *deadCellsExtPtr_;

    label counter = 0;
    forAll (gammaI, cellI)
    {
        if (gammaI[cellI] < SMALL)
        {
            deadCellsExt[counter++] = cellI;
        }
    }
}


void Foam::immersedBoundaryFvPatch::makeDeadFaces() const
{
    if (debug)
    {
        InfoIn("void immersedBoundaryFvPatch::makeDeadFaces() const")
            << "create list of dead cells"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (deadFacesPtr_)
    {
        FatalErrorIn("immersedBoundaryFvPatch::makeDeadFaces() const")
            << "list of dead cells already exist"
            << abort(FatalError);
    }

    deadFacesPtr_ = new labelList(mesh_.nFaces());
    labelList& df = *deadFacesPtr_;
    label nDf = 0;

    const volScalarField& gE = gammaExt();
    const scalarField& gammaExtI = gE.internalField();

    const volScalarField::GeometricBoundaryField& gEPatches =
        gE.boundaryField();

    const labelList& owner = mesh_.faceOwner();
    const labelList& neighbour = mesh_.faceNeighbour();

    forAll (neighbour, faceI)
    {
        if (gammaExtI[neighbour[faceI]] + gammaExtI[owner[faceI]] < SMALL)
        {
            df[nDf] = faceI;
            nDf++;
        }
    }

    forAll (gEPatches, patchI)
    {
        const label start = mesh_.boundaryMesh()[patchI].start();

        if (gEPatches[patchI].coupled())
        {
            scalarField gammaExtOwn =
                gEPatches[patchI].patchInternalField();

            scalarField gammaExtNei =
                gEPatches[patchI].patchNeighbourField();

            forAll (gammaExtOwn, patchFaceI)
            {
                if
                (
                    gammaExtNei[patchFaceI] + gammaExtOwn[patchFaceI] < SMALL
                )
                {
                    df[nDf] = start + patchFaceI;
                    nDf++;
                }
            }
        }
        else
        {
            scalarField gammaExtOwn =
                gEPatches[patchI].patchInternalField();

            forAll (gammaExtOwn, patchFaceI)
            {
                if (gammaExtOwn[patchFaceI] < SMALL)
                {
                    df[nDf] = start + patchFaceI;
                    nDf++;
                }
            }
        }
    }

    df.setSize(nDf);
}


void Foam::immersedBoundaryFvPatch::makeLiveCells() const
{
    if (debug)
    {
        InfoIn("void immersedBoundaryFvPatch::makeLiveCells() const")
            << "create list of live cells"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (liveCellsPtr_)
    {
        FatalErrorIn("immersedBoundaryFvPatch::makeLiveCells() const")
            << "list of live cells already exist"
            << abort(FatalError);
    }

    const scalarField& gammaI = gamma().internalField();

    liveCellsPtr_ = new labelList(label(sum(gammaI)), -1);
    labelList& liveCells = *liveCellsPtr_;

    label counter = 0;
    forAll (gammaI, cellI)
    {
        if (gammaI[cellI] > (1.0 - SMALL))
        {
            liveCells[counter++] = cellI;
        }
    }
}


void Foam::immersedBoundaryFvPatch::makeIbCellSizes() const
{
    if (debug)
    {
        InfoIn("void immersedBoundaryFvPatch::makeIbCellsSize() const")
            << "create average sizes of immersed boundary cells"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (ibCellSizesPtr_)
    {
        FatalErrorIn("immersedBoundaryFvPatch::makeIbCellsSize() const")
            << "average sizes of immersed boundary cells already exist"
            << abort(FatalError);
    }

    ibCellSizesPtr_ = new scalarField(ibPoints().size(), 0.0);
    scalarField& delta = *ibCellSizesPtr_;

    if (mesh_.nGeometricD() == 3)
    {
        // Create a list of volumes with mapping to contain only IB cells
        scalarField V(mesh_.V().field(), ibCells());

        delta = Foam::pow(V, 1.0/3.0);
    }
    else
    {
        // For 2-D simulations with the immersed boundary method the geometry
        // needs to be aligned with the z-direction.
        // Having the x- or y-direction as empty is not allowed because of
        // the way the polynomials are expanded
        const Vector<label>& directions = mesh_.geometricD();

        if (directions[0] == -1 || directions[1] == -1)
        {
            FatalErrorIn("immersedBoundaryFvPatch::makeIbCellsSize() const")
                << "For 2-D simulations with the immersed boundary method "
                << "the geometry needs to be aligned with the z-direction.  "
                << "Having the x- or y-direction as empty is not allowed "
                << "because of the way the polynomials are expanded."
                << abort(FatalError);
        }

        scalar thickness = 0.0;

        for (direction dir = 0; dir < directions.nComponents; dir++)
        {
            if (directions[dir] == -1)
            {
                thickness = mesh_.bounds().span()[dir];
                break;
            }
        }

        // Field created with mapping for IB cells only
        delta = sqrt(scalarField(mesh_.V().field(), ibCells())/thickness);
    }
}


void Foam::immersedBoundaryFvPatch::makeIbSf() const
{
    if (debug)
    {
        InfoIn("void immersedBoundaryFvPatch::makeIbSf() const")
            << "creating ibSf and ibMagSf field"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (ibSfPtr_ || ibMagSfPtr_)
    {
        FatalErrorIn("immersedBoundaryFvPatch::makeIbSf() const")
            << "ibSf and ibMagSf field already exist"
            << abort(FatalError);
    }

    const vectorField& areas = mesh_.faceAreas();

    // Field created with mapping for IB cells only
    ibSfPtr_ = new vectorField(areas, ibFaces());
    ibMagSfPtr_ = new scalarField(mag(*ibSfPtr_));
}


void Foam::immersedBoundaryFvPatch::makeIbDelta() const
{
    if (debug)
    {
        InfoIn("void immersedBoundaryFvPatch::makeIbDelta() const")
            << "creating delta field"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (ibDeltaPtr_)
    {
        FatalErrorIn("immersedBoundaryFvPatch::makeIbDelta() const")
            << "delta field already exist"
            << abort(FatalError);
    }

    const vectorField& C = mesh_.cellCentres();

    // Field created with mapping for IB cells only
    ibDeltaPtr_ =
        new scalarField(mag(ibPoints() - vectorField(C, ibCells())) + SMALL);
}


void Foam::immersedBoundaryFvPatch::makeIbSamplingPointDelta() const
{
    if (debug)
    {
        InfoIn
        (
            "void immersedBoundaryFvPatch::makeIbSamplingPointDelta() const"
        )   << "creating sampling point delta field"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (ibSamplingPointDeltaPtr_)
    {
        FatalErrorIn
        (
            "immersedBoundaryFvPatch::makeIbSamplingPointDelta() const"
        )   << "sampling point delta field already exist"
            << abort(FatalError);
    }

    // Field created with mapping for IB cells only
    ibSamplingPointDeltaPtr_ =
        new scalarField(mag(ibPoints() - ibSamplingPoints()) + SMALL);
}


void Foam::immersedBoundaryFvPatch::makeTriSf() const
{
    if (debug)
    {
        InfoIn("void immersedBoundaryFvPatch::makeTriSf() const")
            << "creating delta field"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (triSfPtr_)
    {
        FatalErrorIn("immersedBoundaryFvPatch::makeTriSf() const")
            << "triSf field already exist"
            << abort(FatalError);
    }


    const triSurface& triMesh = ibPolyPatch_.ibMesh();
    const pointField& triMeshPoints = triMesh.points();

    triSfPtr_ = new vectorField(triMesh.size());
    vectorField& Sf = *triSfPtr_;

    forAll (triMesh, faceI)
    {
        Sf[faceI] = triMesh[faceI].normal(triMeshPoints);
    }

    if (!ibPolyPatch_.internalFlow())
    {
        // Tri surface points the wrong way; flip all area vectors
        Sf *= -1;
    }
}

void Foam::immersedBoundaryFvPatch::makeTriNormals() const
{
    if (debug)
    {
        InfoIn
            (
                "void immersedBoundaryFvPatch::makeTriNormals() const"
            )   << "create Tri Surface normal vectors"
                << "for immersed boundary " << name()
                << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (triNormalsPtr_)
    {
        FatalErrorIn
            (
                "immersedBoundaryFvPatch::makeTriNormals() const"
            )
            << "Tri Surface normals already exist"
                << "for immersed boundary " << name()
                << abort(FatalError);
    }
 

    const vectorField& Sfb = triSf();
	const triSurface& triMesh = ibPolyPatch_.ibMesh();

    triNormalsPtr_ = new vectorField(Sfb.size());
    vectorField& Norms = *triNormalsPtr_;

    scalarField sA = mag(Sfb);

    //check zero face area
    forAll (sA, triI)
    {
		if(sA[triI]<SMALL)
		{
	        Info<<triI<<" "<<Sfb[triI]<<" "<<sA[triI]<<" "<<triMesh.faceCentres()[triI]<<endl;
		}
    }

    forAll (Sfb, triI)
    {
	    Norms[triI] = Sfb[triI]/sA[triI];
    }
 
}

//make maximum slope at each bed cell
void Foam::immersedBoundaryFvPatch::makeTriSlopes() const
{
    if (debug)
    {
        InfoIn
            (
                "void immersedBoundaryFvPatch::makeTriSlope() const"
            )   << "create Tri Surface maximum slope direction at each cell"
                << "for immersed boundary " << name()
                << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (triSlopesPtr_)
    {
        FatalErrorIn
            (
                "immersedBoundaryFvPatch::makeTriSlopes() const"
            )
            << "Tri Surface maximum slope direction already exist"
                << "for immersed boundary " << name()
                << abort(FatalError);
    }

    const vectorField& sNom = triNormals();
    triSlopesPtr_ = new vectorField(sNom.size());
    vectorField& Slopes = *triSlopesPtr_;

    IOdictionary ibmDict
        (
            IOobject
            (
                "ibmDict",
                mesh_.time().system(),
                mesh_,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        );

    // read g vector
    vector gDirection(ibmDict.lookup("gDirection"));
 
    forAll (Slopes, triI)
    {
        //calculate maximum slope       s=-(g*n)/(n*n)*n+g
        Slopes[triI] = -(gDirection&sNom[triI])/(sNom[triI]&sNom[triI])*sNom[triI]
                  +gDirection;

        //if(mag(Slopes[triI])!=0)
        {
             Slopes[triI]/=mag(Slopes[triI])+SMALL;
        }
        // tell if the angle between g and s is larger than 90 degrees
        scalar sig;
        sig=Slopes[triI]&gDirection;
        if(sig<0)
        {
             Slopes[triI]*=-1.0;
        }
    }
 
}

//
void Foam::immersedBoundaryFvPatch::makeTriFacesInMeshType() const
{ 	const double Oldtime0=mesh_.time().elapsedCpuTime();

    if (debug)
    {
        InfoIn
            (
                "void immersedBoundaryFvPatch::makeTriFacesInMeshType() const"
            )   << "create Tri Surface in mesh type: outside/inside/cut by patch/wall/coupled "
                << "for immersed boundary " << name()
                << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (triFacesInMeshTypePtr_ || triFacesInMeshTypeVecPtr_ || triFacesInMeshPtr_)
    {
        FatalErrorIn
            (
                "immersedBoundaryFvPatch::makeTriFacesInMeshType() const"
            )
            << "create Tri Surface in mesh type: outside/inside/cut by patch/wall/coupled"
                << "for immersed boundary " << name()
                << abort(FatalError);
    }

    const vectorField& triCf = this->triCf();
	const pointField& pts = ibMesh().localPoints();
	const List<labelledTri>& triFaces(ibMesh().localFaces());

    triFacesInMeshTypePtr_ = new scalarListList(triCf.size());
    scalarListList& triType = *triFacesInMeshTypePtr_;

    triFacesInMeshTypeVecPtr_ = new vectorListList(triCf.size());
    vectorListList& triTypeVec = *triFacesInMeshTypeVecPtr_;

    triFacesInMeshPtr_ = new labelList();
    labelList& triFacesInMesh = *triFacesInMeshPtr_;
	labelHashSet triFacesInMeshSet; 

	forAll(triType, triI)
	{
		label cellID=mesh_.findCell(triCf[triI]);

		scalarList& curTriType=triType[triI];
		vectorList& curTriTypeVec=triTypeVec[triI];
		curTriType.setSize(0);
		curTriTypeVec.setSize(0); 
		if(cellID<0)
		{

			const labelledTri& triFace=triFaces[triI];
			forAll(triFace,I)
			{
				curTriType.append(0);//outside
				curTriTypeVec.append(vector::zero);

				const label& ptI=triFace[I];
				label ptInCellID = mesh_.findCell(pts[ptI]);

				// check if points of each triangle is in mesh
				if(ptInCellID<0)
				{
					continue;
				}
				else
				{					
					const cell& ptInCell=mesh_.cells()[ptInCellID];
					forAll(ptInCell, I)
					{
						label patchID = mesh_.boundaryMesh().whichPatch(ptInCell[I]);
						

						// check the patch
						if (patchID<0)
						{
							continue;
							//curTriType[I]=1;//inside
						}
						const fvPatch& curPatch = mesh_.boundary()[patchID];
						tmp<vectorField> patchNf = curPatch.nf();
						if (isType<wallFvPatch>(curPatch))
						{
							label faceI=curPatch.patch().whichFace(ptInCell[I]);
							curTriTypeVec.append(patchNf()[faceI]);
							curTriType.append(3);//cut by wall

						}
						else if (curPatch.coupled())
						{
							label faceI=curPatch.patch().whichFace(ptInCell[I]);
							curTriTypeVec.append(patchNf()[faceI]);
							curTriType.append(4);//cut by coupled patches
							if(!triFacesInMeshSet.found(triI))
							{
								triFacesInMeshSet.insert(triI);
							}
						}
						else if (curPatch.type()!="empty")
						{
							label faceI=curPatch.patch().whichFace(ptInCell[I]);
							curTriTypeVec.append(patchNf()[faceI]);
							curTriType.append(2);//cut by patches
							if(!triFacesInMeshSet.found(triI))
							{
								triFacesInMeshSet.insert(triI);
							}
						}
					}
					if(curTriType.size()<2)
					{
						Pout<<"Triangle is too big at "<<pts[ptI]<<endl;
					}
				}
			}
			continue;// triFace center is not in mesh
		}

		if(!triFacesInMeshSet.found(triI))
		{
			triFacesInMeshSet.insert(triI);
		}
		else
		{
			FatalErrorIn
	        (
	            "immersedBoundaryFvPatch::makeTriFacesInMeshType() const"
	        )
	            << "something goes wrong when check if a triangle is in mesh "
	            << "for immersed boundary " << name()
	            << abort(FatalError);
		}
		const cell& cell=mesh_.cells()[cellID];
		forAll(cell, I)
		{
			label patchID = mesh_.boundaryMesh().whichPatch(cell[I]);

			//tmp<vectorField > tpatchNf(new vectorField(curPatch.nf()));
			//const vectorField& patchNf = tpatchNf();
			if (patchID<0) //triangle center is not inside boundary cells
			{
				continue;
				//curTriType[I]=1;//inside
			}
			/*
			// this is to find the adjacent triangle to boundary, 
			// however, we need triangle with center in adjacent cell to boundary
			else
			{
				const labelledTri& triFace=triFaces[triI];
				bool allInside=true;
				forAll(triFace,I)
				{
					const label& ptI=triFace[I];
					if (mesh_.findCell(pts[ptI])<0)
					{
						allInside=false;
					}
				}
				if(allInside)
				{
					continue;
				}
			}*/
			const fvPatch& curPatch = mesh_.boundary()[patchID];
			tmp<vectorField> patchNf = curPatch.nf();
			if (isType<wallFvPatch>(curPatch))
			{
				label faceI=curPatch.patch().whichFace(cell[I]);
				curTriTypeVec.append(patchNf()[faceI]);
				curTriType.append(3);//cut by wall
			}
			else if (curPatch.coupled())
			{
				label faceI=curPatch.patch().whichFace(cell[I]);
				curTriTypeVec.append(patchNf()[faceI]);
				curTriType.append(4);//cut by coupled patches
			}
			else if (curPatch.type()!="empty")
			{
				label faceI=curPatch.patch().whichFace(cell[I]);
				curTriTypeVec.append(patchNf()[faceI]);
				curTriType.append(2);//cut by patches
			}
		}
		if(curTriType.size()<SMALL)
		{
			curTriType.setSize(1);
			curTriTypeVec.setSize(1);
			curTriType[0]=1;//inside
			curTriTypeVec[0]=vector::zero;
			continue;
		}
	}

	triFacesInMesh = triFacesInMeshSet.toc();
	const double Oldtime2=mesh_.time().elapsedCpuTime();
	Info<<"makeTriFacesInMeshType Executation Time = "<<Oldtime2-Oldtime0<< " s"<<endl;
 
}

// Find the cell with the nearest cell centre
Foam::label Foam::immersedBoundaryFvPatch::findNearestCell
(
    const point& location
) const
{
    const vectorField& C = mesh_.cellCentres();

    const scalarField& gammaExtI = gammaExt().internalField();

    label nearestCellI = -1;
    scalar minProximity = GREAT;

    for (label cellI = 0; cellI < C.size(); cellI++)
    {
        if (gammaExtI[cellI] > SMALL)
        {
            scalar proximity = magSqr(C[cellI] - location);
            if (proximity < minProximity)
            {
                nearestCellI = cellI;
                minProximity = proximity;
            }
        }
    }

    return nearestCellI;
}

//- Return extended cell-cell addressing for one cell
void Foam::immersedBoundaryFvPatch::findCellCells
(
    const vector& pt,
    const label cellID,
    labelList& cellCells
) const
{

    const labelListList& cellPoints = mesh_.cellPoints();
    const labelListList& pointCells = mesh_.pointCells();

    const scalarField& gammaExtI = gammaExt().internalField();

    labelHashSet cellSet;
    cellSet.insert(cellID);    

    // First row
    const labelList& curCellPoints = cellPoints[cellID];

    forAll (curCellPoints, pointI)
    {
        label curPoint = curCellPoints[pointI];
        const labelList& curPointCells = pointCells[curPoint];

        forAll (curPointCells, cI)
        {
            if (gammaExtI[curPointCells[cI]] > SMALL)
            {
                if (!cellSet.found(curPointCells[cI]))
                {
                    cellSet.insert(curPointCells[cI]);
                }
            }
        }
    }

    labelList curCells = cellSet.toc();

    // Second and other rows
    for (label nRows = 1; nRows < maxCellCellRows_; nRows++)
    {
        curCells = cellSet.toc();

        forAll (curCells, cellI)
        {
            label curCell = curCells[cellI];
            const labelList& curCellPoints = cellPoints[curCell];

            forAll (curCellPoints, pointI)
            {
                label curPoint = curCellPoints[pointI];
                const labelList& curPointCells = pointCells[curPoint];

                forAll (curPointCells, cI)
                {
                    if (gammaExtI[curPointCells[cI]] > SMALL)
                    {
                        if (!cellSet.found(curPointCells[cI]))
                        {
                            cellSet.insert(curPointCells[cI]);
                        }
                    }
                }
            }
        }

		if(cellSet.size()>80)
		{
			break;
		}

    }

    // Erase current cell
    cellSet.erase(cellID);

    // Sorting cells
    const vectorField& C = mesh_.cellCentres();

    curCells = cellSet.toc();
    scalarField distances(curCells.size(), 0);

    forAll (distances, cI)
    {
        distances[cI] =
            mag(C[curCells[cI]] - pt);
    }

    SortableList<scalar> sortedDistances(distances);

    labelList sortedCells(curCells.size(), -1);
    //labelList sortedCells(min(45,curCells.size()), -1);

    for (label i = 0; i < sortedCells.size(); i++)
    {
        sortedCells[i] =
            curCells[sortedDistances.indices()[i]];
    }
 
    cellCells = sortedCells;
}


Foam::scalar Foam::immersedBoundaryFvPatch::cellSize(label cellID) const
{
    scalar delta;

    if (mesh_.nGeometricD() == 3)
    {
        delta = Foam::pow(mesh_.V().field()[cellID], 1.0/3.0);
    }
    else
    {
        scalar thickness = 0.0;
        const Vector<label>& directions = mesh_.geometricD();
        for (direction dir = 0; dir < directions.nComponents; dir++)
        {
            if (directions[dir] == -1)
            {
                thickness = mesh_.bounds().span()[dir];
                break;
            }
        }

        delta = sqrt(mesh_.V().field()[cellID]/thickness);
    }

    return delta;
}


Foam::scalar Foam::immersedBoundaryFvPatch::cellProjection
(
    label cellID,
    const vector& dir
) const
{
    const vectorField& C = mesh_.cellCentres();
    const vectorField& Cf = mesh_.faceCentres();
    const vectorField& Sf = mesh_.faceAreas();

    const labelList& cellFaces = mesh_.cells()[cellID];

    scalar area = 0;

    forAll (cellFaces, faceI)
    {
        label curFace = cellFaces[faceI];

        vector curSf = Sf[curFace];

        if ((curSf & (Cf[curFace] - C[cellID])) < 0)
        {
            curSf *= -1;
        }

        if ((curSf&dir) > 1)
        {
            area += (curSf & dir);
        }
    }

    return area;
}

void Foam::immersedBoundaryFvPatch::makeAdjacentIbCells() const
{
    labelHashSet adjacentIbCellSet;

    const unallocLabelList& owner = mesh_.owner();
    const unallocLabelList& neighbour = mesh_.neighbour();

    const volScalarField g = gamma();
    const scalarField& gammaI = g.internalField();

/* //find all fluid cells
    forAll (gammaI, cellI)
    {
        if (gammaI[cellI] > SMALL)
        { 
            if (!adjacentIbCellSet.found(cellI))
            {
                adjacentIbCellSet.insert(cellI);
            }
        }
    }
*/

    forAll (neighbour, faceI)
    {
        if (mag(gammaI[neighbour[faceI]] - gammaI[owner[faceI]]) > SMALL)
        {
            if (gammaI[owner[faceI]] > SMALL)
            {
                if (!adjacentIbCellSet.found(owner[faceI]))
                {
                    adjacentIbCellSet.insert(owner[faceI]);
                }
            }
            else
            {
                if (!adjacentIbCellSet.found(neighbour[faceI]))
                {
                    adjacentIbCellSet.insert(neighbour[faceI]);
                }
            }
        }
    }



    forAll (g.boundaryField(), patchI)
    {

        if (g.boundaryField()[patchI].coupled())
        {
            scalarField gammaOwn =
                g.boundaryField()[patchI].patchInternalField();

            scalarField gammaNei =
                g.boundaryField()[patchI].patchNeighbourField();

            const unallocLabelList& fCells =
                mesh_.boundary()[patchI].faceCells();

            forAll (gammaOwn, faceI)
            {
                if
                    (   
                        mag(gammaNei[faceI] - gammaOwn[faceI])
                        > SMALL
                    )
                {
                    if (gammaOwn[faceI] > SMALL)
                    {
                        if (!adjacentIbCellSet.found(fCells[faceI]))
                        {   
                            adjacentIbCellSet.insert(fCells[faceI]);
                        }
                    }
                    else if (2*gammaOwn.size() == fCells.size())
                    {

                        if
                            (
                                !adjacentIbCellSet.found
                                (
                                    fCells[gammaOwn.size() + faceI]
                                )
                            )
                        {
                            adjacentIbCellSet.insert
                                (
                                    fCells[gammaOwn.size() + faceI]
                                );
                        }
                    }
                }
            }
        }
        else
        {
            scalarField gammaOwn =
                g.boundaryField()[patchI].patchInternalField();

        }
    }

    adjacentIbCellsPtr_ = new labelList(adjacentIbCellSet.toc());
                                                                          
}


void Foam::immersedBoundaryFvPatch::makeAdjacentIbPointsAndNormals() const
{
 

    if (debug)
    {
        InfoIn
        (
            "void immersedBoundaryFvPatch::makeAdjacentIbPointsAndNormals() const"
        )   << "create immersed  boundary points and normals "
            << "for immersed boundary " << name()
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (adjacentIbPointsPtr_ || adjacentIbNormalsPtr_ || adjacentIbDeltaPtr_)
    {
        FatalErrorIn
        (
            "immersedBoundaryFvPatch::makeAdjacentIbPointsAndNormals() const"
        )
            << "immersed boundary points and normals already exist"
            << "for immersed boundary " << name()
            << abort(FatalError);
    }

    // Find average cell dimension
    const labelList& aibc = adjacentIbCells();

    scalarField delta(aibc.size(), 0.0);

    forAll (delta, cellI)
    {
        delta[cellI] = cellSize(aibc[cellI]);
    }

    // Find nearest triSurface point for each interface cell centre
    adjacentIbPointsPtr_ = new vectorField(aibc.size(), vector::zero);
    adjacentIbNormalsPtr_ = new vectorField(aibc.size(), vector::zero);
 

    vectorField& adjacentIbPoints = *adjacentIbPointsPtr_;
    vectorField& adjacentIbNormals = *adjacentIbNormalsPtr_;
 
    const vectorField& C = mesh_.C().internalField();

    // Get adjacent IB cell centres
    vectorField adjacentIbCellCentres(C, aibc);

    const triSurfaceSearch& tss = ibPolyPatch_.triSurfSearch();

    forAll (aibc, cellI)
    {
        // Adjust search span if needed.  HJ, 14/Dec/2012
        vector span
        (
            2*radiusFactor_*delta[cellI],
            2*radiusFactor_*delta[cellI],
            2*radiusFactor_*delta[cellI]
        );

        pointIndexHit pih = tss.nearest(adjacentIbCellCentres[cellI], span);

        if (pih.hit())
        {
            adjacentIbPoints[cellI] = pih.hitPoint();
            adjacentIbNormals[cellI] =
                triSurfaceTools::surfaceNormal
                (
                    ibPolyPatch_.ibMesh(),
                    pih.index(),
                    pih.hitPoint()
                );

            // Note: ibNormals point OUT of the domain
            if (!internalFlow())
            {
                adjacentIbNormals[cellI] *= -1;
            }
        }
        else
        {
            FatalErrorIn
            (
                "immersedBoundaryFvPatch::makeAdjacentIbPointsAndNormals() const"
            )   << "Can't find nearest triSurface point for cell "
                << aibc[cellI] << ", "
                << mesh_.cellCentres()[aibc[cellI]]
                << ".  Hit data = " << pih << nl
                << abort(FatalError);
        }

        if
        (
            mesh_.nGeometricD() < 3
         && mag(adjacentIbCellCentres[cellI].z() - adjacentIbPoints[cellI].z()) > SMALL
        )
        {
            WarningIn
            (
                "immersedBoundaryFvPatch::makeAdjacentIbPointsAndNormals() const"
            )   << "Intersection point is not on symmetry plane " << nl
                << "C = " << adjacentIbCellCentres[cellI]
                <<  " D = " <<  adjacentIbPoints[cellI] << nl
                << "for 2-D geometry.  Adjusting" << endl;

               adjacentIbPoints[cellI].z() = adjacentIbCellCentres[cellI].z();
        }
    }

    // Field created with mapping for IB cells only
    adjacentIbDeltaPtr_ =
        new scalarField(mag(adjacentIbPoints - adjacentIbCellCentres) + SMALL);
}


    // complicated due to MPI, by Y.C. Xu 2017/6
void Foam::immersedBoundaryFvPatch::makeIbCellPts() const  
{
    if (debug)
    {
        InfoIn("void immersedBoundaryFvPatch::makeIbCellPts() const")
            << "add IB pts to ib cells"
            << endl;
    }
    
    // It is an error to attempt to recalculate
    // if the pointer is already set
    if
    (
        ibCellPtsPtr_ 
     || ibCellProcPtsPtr_ 
     || ibProcPtsPtr_
    )
    {
        FatalErrorIn
        (
            "immersedBoundaryFvPatch::makeIbCellPts() const"
        )   << "cell-cell addressing already exists"
            << abort(FatalError);
    }
    const labelList& ibc = ibCells();
    const labelListList& ibcc = ibCellCells();
    const List<List<labelPair> >& cellProcCells = ibCellProcCells();
    const vectorField& pts = ibPoints();
 
    ibCellPtsPtr_ = new labelListList(ibc.size());
    labelListList& ibCellPts = *ibCellPtsPtr_;

    ibProcPtsPtr_ = new vectorListList(Pstream::nProcs());
    vectorListList& procPtsCentres = *ibProcPtsPtr_;

    ibCellProcPtsPtr_ = new List<List<labelPair> >(ibc.size());
    List<List<labelPair> >& cellProcPts = *ibCellProcPtsPtr_;

    labelList ibCellIndicator(mesh_.nCells(), -1);

    forAll (ibc, cellI)
    {
        ibCellIndicator[ibc[cellI]] = cellI;
    }

    forAll (ibc, cellI)
    {
        const labelList& curAddr = ibcc[cellI];

        labelHashSet ibPtsSet;
        labelList& curPtPts = ibCellPts[cellI];
 
        if(!ibPtsSet.found(cellI))
        {
            ibPtsSet.insert(cellI);
        }
 
        /*if(!ibPtsSet.found(cellI) and ibCellIndicator[ibc[cellI]]>SMALL)
        {
            ibPtsSet.insert(cellI);
        }*/

        forAll (curAddr, ccI)
        {
            const label& curAddrI=curAddr[ccI];
 
            const label& insertedPointI(ibCellIndicator[curAddrI]);

            if(!ibPtsSet.found(insertedPointI) and insertedPointI>SMALL)
            {

                ibPtsSet.insert(insertedPointI);

            }
        }
 

        curPtPts =  ibPtsSet.toc();
    }    
    
    if (Pstream::parRun())
    {   
        // Send and receive ibc labelHashSet
        /*labelHashSetList procIbCellSet(Pstream::nProcs());

        procIbCellSet[Pstream::myProcNo()] = ibCellSet;
   
        Pstream::gatherList(procIbCellSet);
        Pstream::scatterList(procIbCellSet);

        // Send and receive IB point centres labelHashSet
        vectorListList procIbPts(Pstream::nProcs());*/

        // Send and receive ibc labelHashSet
        labelListList procIbCellIndicator(Pstream::nProcs());

        procIbCellIndicator[Pstream::myProcNo()] = ibCellIndicator;
   
        Pstream::gatherList(procIbCellIndicator);
        Pstream::scatterList(procIbCellIndicator);

        // Send and receive IB point centres labelHashSet
        vectorListList procIbPts(Pstream::nProcs());
        


        procIbPts[Pstream::myProcNo()] = pts;
   
        Pstream::gatherList(procIbPts);
        Pstream::scatterList(procIbPts);

        procPtsCentres.setSize(cellProcCells.size());

        forAll (cellProcCells, cellI)
        {
            const List<labelPair>& interpProcCells = cellProcCells[cellI];
            List<labelPair>& interpProcPts = cellProcPts[cellI];

            vectorList& curProcPtsCentres = procPtsCentres[cellI];

            labelHashSet interpProcPtsFirstSet;
            labelHashSet interpProcPtsSecondSet;

            forAll (interpProcCells, cProcI)
            {
                label procI = interpProcCells[cProcI].first();
                label pointI = interpProcCells[cProcI].second();
                
                labelList& cellIndicator = procIbCellIndicator[procI];

                const label& insertedPointI(cellIndicator[pointI]);
                if(!interpProcPtsSecondSet.found(insertedPointI) and insertedPointI>SMALL)
                {
                    
                    interpProcPtsFirstSet.insert(procI);
                    interpProcPtsSecondSet.insert(insertedPointI);
                    curProcPtsCentres.append(procIbPts[procI][insertedPointI]);            
                }                              
            }

            labelList tmpInterpProcPtsFirst = interpProcPtsFirstSet.toc();
            labelList tmpInterpProcPtsSecond = interpProcPtsSecondSet.toc();
             
            forAll (interpProcPts, cProcI)
            {
                interpProcPts[cProcI].first() = tmpInterpProcPtsFirst[cProcI];
                interpProcPts[cProcI].second() = tmpInterpProcPtsSecond[cProcI];        
            }            
        }
    
    }
}

void Foam::immersedBoundaryFvPatch::clearOut()
{
    Info<<"clearOut"<<endl;
    deleteDemandDrivenData(gammaPtr_);
    deleteDemandDrivenData(gammaExtPtr_);
    deleteDemandDrivenData(sGammaPtr_);
    deleteDemandDrivenData(ibCellsPtr_);
    deleteDemandDrivenData(ibFacesPtr_);
    deleteDemandDrivenData(ibFaceCellsPtr_);
    deleteDemandDrivenData(ibFaceFlipsPtr_);
    deleteDemandDrivenData(ibInsideFacesPtr_);
    deleteDemandDrivenData(ibInternalFacesPtr_);
    deleteDemandDrivenData(ibPointsPtr_);
    deleteDemandDrivenData(ibNormalsPtr_);
    deleteDemandDrivenData(hitFacesPtr_);
    deleteDemandDrivenData(ibSamplingPointsPtr_);

    deleteDemandDrivenData(ibSamplingWeightsPtr_);
    deleteDemandDrivenData(ibSamplingProcWeightsPtr_);

    deleteDemandDrivenData(cellsToTriAddrPtr_);
    deleteDemandDrivenData(cellsToTriWeightsPtr_);

    deleteDemandDrivenData(ibCellCellsPtr_);
    deleteDemandDrivenData(ibProcCellsPtr_);
    deleteDemandDrivenData(ibProcIbCellsPtr_);
    deleteDemandDrivenData(ibProcCentresPtr_);
    deleteDemandDrivenData(ibProcGammaPtr_);
    deleteDemandDrivenData(ibCellProcCellsPtr_);
    deleteDemandDrivenData(deadCellsPtr_);
    deleteDemandDrivenData(deadCellsExtPtr_);
    deleteDemandDrivenData(deadFacesPtr_);
    deleteDemandDrivenData(liveCellsPtr_);
    deleteDemandDrivenData(ibCellSizesPtr_);

    deleteDemandDrivenData(invDirichletMatricesPtr_);
    deleteDemandDrivenData(invNeumannMatricesPtr_);

    deleteDemandDrivenData(ibSfPtr_);
    deleteDemandDrivenData(ibMagSfPtr_);
    deleteDemandDrivenData(ibDeltaPtr_);
    deleteDemandDrivenData(ibSamplingPointDeltaPtr_);

    deleteDemandDrivenData(adjacentIbCellsPtr_);
    deleteDemandDrivenData(adjacentIbPointsPtr_);
    deleteDemandDrivenData(adjacentIbNormalsPtr_);
    deleteDemandDrivenData(adjacentIbDeltaPtr_);
    deleteDemandDrivenData(ibCellPtsPtr_);
    deleteDemandDrivenData(ibCellProcPtsPtr_);
    deleteDemandDrivenData(ibProcPtsPtr_);
    deleteDemandDrivenData(ibNewSamplingWeightsPtr_);
    deleteDemandDrivenData(ibNewSamplingProcWeightsPtr_);
    deleteDemandDrivenData(ibNewPtsWeightsPtr_);
    deleteDemandDrivenData(ibNewPtsProcWeightsPtr_);

    deleteDemandDrivenData(triFacesToTriPointsWeightsPtr_);

    deleteDemandDrivenData(triSfPtr_);
    deleteDemandDrivenData(triNormalsPtr_);
    deleteDemandDrivenData(triSlopesPtr_);
	deleteDemandDrivenData(triFacesInMeshPtr_);
	deleteDemandDrivenData(triFacesInMeshTypePtr_);
	deleteDemandDrivenData(triFacesInMeshTypeVecPtr_);
 
	deleteDemandDrivenData(triFacesToTriEdgesWeightsPtr_);

}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::immersedBoundaryFvPatch::initMovePoints()
{}


void Foam::immersedBoundaryFvPatch::movePoints()
{
    clearOut();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immersedBoundaryFvPatch::immersedBoundaryFvPatch
(
    const polyPatch& patch,
    const fvBoundaryMesh& bm
)
:
    fvPatch(patch, bm),
    ibPolyPatch_(refCast<const immersedBoundaryPolyPatch>(patch)),
    mesh_(bm.mesh()),
    ibUpdateTimeIndex_(-1),
    angleFactor_(0.0),
    radiusFactor_(3.0),
    maxCellCellRows_(3.0),
    distFactor_(1.5),
    gammaPtr_(NULL),
    gammaExtPtr_(NULL),
    sGammaPtr_(NULL),
    ibCellsPtr_(NULL),
    ibFacesPtr_(NULL),
    ibFaceCellsPtr_(NULL),
    ibFaceFlipsPtr_(NULL),
    ibInsideFacesPtr_(NULL),
    ibInternalFacesPtr_(NULL),
    ibPointsPtr_(NULL),
    ibNormalsPtr_(NULL),
    hitFacesPtr_(NULL),
    ibSamplingPointsPtr_(NULL),
    ibSamplingWeightsPtr_(NULL),
    ibSamplingProcWeightsPtr_(NULL),
    cellsToTriAddrPtr_(NULL),
    cellsToTriWeightsPtr_(NULL),
    ibCellCellsPtr_(NULL),
    ibProcCellsPtr_(NULL),
    ibProcIbCellsPtr_(NULL),
    ibProcCentresPtr_(NULL),
    ibProcGammaPtr_(NULL),
    ibCellProcCellsPtr_(NULL),
    deadCellsPtr_(NULL),
    deadCellsExtPtr_(NULL),
    deadFacesPtr_(NULL),
    liveCellsPtr_(NULL),
    ibCellSizesPtr_(NULL),
    invDirichletMatricesPtr_(NULL),
    invNeumannMatricesPtr_(NULL),
    ibSfPtr_(NULL),
    ibMagSfPtr_(NULL),
    ibDeltaPtr_(NULL),
    ibSamplingPointDeltaPtr_(NULL),
    adjacentIbCellsPtr_(NULL),
    adjacentIbPointsPtr_(NULL),
    adjacentIbNormalsPtr_(NULL),
    adjacentIbDeltaPtr_(NULL),
    ibCellPtsPtr_(NULL),
    ibCellProcPtsPtr_(NULL),
    ibProcPtsPtr_(NULL),
    ibNewSamplingWeightsPtr_(NULL),
    ibNewSamplingProcWeightsPtr_(NULL),
    ibNewPtsWeightsPtr_(NULL),
    ibNewPtsProcWeightsPtr_(NULL),
    triFacesToTriPointsWeightsPtr_(NULL),
    triSfPtr_(NULL),
    triNormalsPtr_(NULL),
    triSlopesPtr_(NULL),
	triFacesInMeshPtr_(NULL),
	triFacesInMeshTypePtr_(NULL),
	triFacesInMeshTypeVecPtr_(NULL),
	triFacesToTriEdgesWeightsPtr_(NULL)

{
    this->makeConstant();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
//make triSurface related data
void Foam::immersedBoundaryFvPatch::makeTriSurfaceData() const
{
    makeTriSf();
    makeTriNormals();
    makeTriSlopes();
}

void Foam::immersedBoundaryFvPatch::makeImmersedBoundaryData() const
{


    Pout << "Start making immersed boundary data ... " << endl;
    makeTriSurfaceData();

    makeGammaExt();

    makeGamma();
    
    makeSGamma();

    makeIbCells();

    makeIbFaces();

    makeTriAddressing();

    makeIbInsideFaces();

    makeIbInternalFaces();

    makeIbPointsAndNormals();

    makeIbSamplingWeights();

    makeIbCellCells();

    makeDeadCells();

    makeLiveCells();

    Pout << "Finished making immersed boundary data ... " << endl;
}




const Foam::volScalarField& Foam::immersedBoundaryFvPatch::gamma() const
{
    if (!gammaPtr_)
    {
        makeGamma();
    }

    return *gammaPtr_;
}


const Foam::volScalarField& Foam::immersedBoundaryFvPatch::gammaExt() const
{
    if (!gammaExtPtr_)
    {
        makeGammaExt();
    }

    return *gammaExtPtr_;
}


const Foam::surfaceScalarField& Foam::immersedBoundaryFvPatch::sGamma() const
{
    if (!sGammaPtr_)
    {
        makeSGamma();
    }

    return *sGammaPtr_;
}


const Foam::labelList& Foam::immersedBoundaryFvPatch::ibCells() const
{
    if (!ibCellsPtr_)
    {
        makeIbCells();
    }

    return *ibCellsPtr_;
}


const Foam::labelList& Foam::immersedBoundaryFvPatch::ibFaces() const
{
    if (!ibFacesPtr_)
    {
        makeIbFaces();
    }

    return *ibFacesPtr_;
}


const Foam::labelList& Foam::immersedBoundaryFvPatch::ibFaceCells() const
{
    if (!ibFaceCellsPtr_)
    {
        makeIbFaces();
    }

    return *ibFaceCellsPtr_;
}


const Foam::boolList& Foam::immersedBoundaryFvPatch::ibFaceFlips() const
{
    if (!ibFaceFlipsPtr_)
    {
        makeIbFaces();
    }

    return *ibFaceFlipsPtr_;
}


const Foam::labelList& Foam::immersedBoundaryFvPatch::ibInsideFaces() const
{
    if (!ibInsideFacesPtr_)
    {
        makeIbInsideFaces();
    }

    return *ibInsideFacesPtr_;
}


const Foam::labelList& Foam::immersedBoundaryFvPatch::ibInternalFaces() const
{
    if (!ibInternalFacesPtr_)
    {
        makeIbInternalFaces();
    }

    return *ibInternalFacesPtr_;
}


const Foam::vectorField& Foam::immersedBoundaryFvPatch::ibPoints() const
{
    if (!ibPointsPtr_)
    {
        makeIbPointsAndNormals();
    }

    return *ibPointsPtr_;
}


const Foam::vectorField& Foam::immersedBoundaryFvPatch::ibNormals() const
{
    if (!ibNormalsPtr_)
    {
        makeIbPointsAndNormals();
    }

    return *ibNormalsPtr_;
}


const Foam::labelList& Foam::immersedBoundaryFvPatch::hitFaces() const
{
    if (!hitFacesPtr_)
    {
        makeIbPointsAndNormals();
    }

    return *hitFacesPtr_;
}


const Foam::vectorField&
Foam::immersedBoundaryFvPatch::ibSamplingPoints() const
{
    if (!ibSamplingPointsPtr_)
    {
        makeIbPointsAndNormals();
    }

    return *ibSamplingPointsPtr_;
}


const Foam::labelListList& Foam::immersedBoundaryFvPatch::ibCellCells() const
{
    if (!ibCellCellsPtr_)
    {
        makeIbCellCells();
    }

    return *ibCellCellsPtr_;
}


const Foam::vectorListList&
Foam::immersedBoundaryFvPatch::ibProcCentres() const
{
    if (!ibProcCentresPtr_)
    {
        makeIbCellCells();
    }

    return *ibProcCentresPtr_;
}


const Foam::scalarListList&
Foam::immersedBoundaryFvPatch::ibProcGamma() const
{
    if (!ibProcGammaPtr_)
    {
        makeIbCellCells();
    }

    return *ibProcGammaPtr_;
}


const Foam::List<Foam::List<Foam::labelPair> >&
Foam::immersedBoundaryFvPatch::ibCellProcCells() const
{
    if (!ibCellProcCellsPtr_)
    {
        makeIbCellCells();
    }

    return *ibCellProcCellsPtr_;
}


const Foam::labelListList& Foam::immersedBoundaryFvPatch::ibProcCells() const
{
    if (!ibProcCellsPtr_)
    {
        makeIbCellCells();
    }

    return *ibProcCellsPtr_;
}

const Foam::labelListList& Foam::immersedBoundaryFvPatch::ibProcIbCells() const
{
    if (!ibProcIbCellsPtr_)
    {
        makeIbCellCells();
    }

    return *ibProcIbCellsPtr_;
}


const Foam::labelList& Foam::immersedBoundaryFvPatch::deadCells() const
{
    if (!deadCellsPtr_)
    {
        makeDeadCells();
    }

    return *deadCellsPtr_;
}


const Foam::labelList& Foam::immersedBoundaryFvPatch::deadCellsExt() const
{
    if (!deadCellsExtPtr_)
    {
        makeDeadCellsExt();
    }

    return *deadCellsExtPtr_;
}


const Foam::labelList& Foam::immersedBoundaryFvPatch::deadFaces() const
{
    if (!deadFacesPtr_)
    {
        makeDeadFaces();
    }

    return *deadFacesPtr_;
}


const Foam::labelList& Foam::immersedBoundaryFvPatch::liveCells() const
{
    if (!liveCellsPtr_)
    {
        makeLiveCells();
    }

    return *liveCellsPtr_;
}


const Foam::scalarField& Foam::immersedBoundaryFvPatch::ibCellSizes() const
{
    if (!ibCellSizesPtr_)
    {
        makeIbCellSizes();
    }

    return *ibCellSizesPtr_;
}


const Foam::vectorField& Foam::immersedBoundaryFvPatch::ibSf() const
{
    if (!ibSfPtr_)
    {
        makeIbSf();
    }

    return *ibSfPtr_;
}


const Foam::scalarField& Foam::immersedBoundaryFvPatch::ibMagSf() const
{
    if (!ibMagSfPtr_)
    {
        makeIbSf();
    }

    return *ibMagSfPtr_;
}


const Foam::scalarField& Foam::immersedBoundaryFvPatch::ibDelta() const
{
    if (!ibDeltaPtr_)
    {
        makeIbDelta();
    }

    return *ibDeltaPtr_;
}


const Foam::scalarField&
Foam::immersedBoundaryFvPatch::ibSamplingPointDelta() const
{
    if (!ibSamplingPointDeltaPtr_)
    {
        makeIbSamplingPointDelta();
    }

    return *ibSamplingPointDeltaPtr_;
}


const Foam::vectorField& Foam::immersedBoundaryFvPatch::triSf() const
{
    if (!triSfPtr_)
    {
        makeTriSurfaceData();
    }

    return *triSfPtr_;
}

const Foam::vectorField& Foam::immersedBoundaryFvPatch::triNormals() const
{
    if (!triNormalsPtr_)
    {
        makeTriSurfaceData();
    }

    return *triNormalsPtr_;
}


const Foam::vectorField& Foam::immersedBoundaryFvPatch::triSlopes() const
{
    if (!triSlopesPtr_)
    {
        makeTriSurfaceData();
    }

    return *triSlopesPtr_;
}



const Foam::vectorField& Foam::immersedBoundaryFvPatch::triCf() const
{
    return ibMesh().faceCentres();
}


const Foam::scalarListList& Foam::immersedBoundaryFvPatch::triFacesInMeshType() const
{
	if(!triFacesInMeshTypePtr_)
	{
		makeTriFacesInMeshType();
	}
	
	return *triFacesInMeshTypePtr_;
}

const Foam::vectorListList& Foam::immersedBoundaryFvPatch::triFacesInMeshTypeVec() const
{
	if(!triFacesInMeshTypeVecPtr_)
	{
		makeTriFacesInMeshType();
	}
	
	return *triFacesInMeshTypeVecPtr_;
}

const Foam::labelList&
Foam::immersedBoundaryFvPatch::triFacesInMesh() const
{
	// update by Xu DEC2017, include point of triFace in mesh, to prevent extreme value
	if(!triFacesInMeshPtr_)
	{
		makeTriFacesInMeshType();
	}
	
	return *triFacesInMeshPtr_;
}

const Foam::labelList& Foam::immersedBoundaryFvPatch::adjacentIbCells() const
{
    if (!adjacentIbCellsPtr_)
    {
        makeAdjacentIbCells();
    }

    return *adjacentIbCellsPtr_;
}

const Foam::vectorField& Foam::immersedBoundaryFvPatch::adjacentIbPoints() const
{
    if (!adjacentIbPointsPtr_)
    {
        makeAdjacentIbPointsAndNormals();
    }

    return *adjacentIbPointsPtr_;
}


const Foam::vectorField& Foam::immersedBoundaryFvPatch::adjacentIbNormals() const
{
    if (!adjacentIbNormalsPtr_)
    {
        makeAdjacentIbPointsAndNormals();
    }

    return *adjacentIbNormalsPtr_;
}

const Foam::scalarField& Foam::immersedBoundaryFvPatch::adjacentIbDelta() const
{
    if (!adjacentIbDeltaPtr_)
    {
        makeAdjacentIbPointsAndNormals();
    }

    return *adjacentIbDeltaPtr_;
}

const Foam::labelListList& Foam::immersedBoundaryFvPatch::ibCellPts() const
{
    if (!ibCellPtsPtr_)
    {
        makeIbCellPts();
    }

    return *ibCellPtsPtr_;
}


const Foam::List<Foam::List<Foam::labelPair> >&
Foam::immersedBoundaryFvPatch::ibCellProcPts() const
{
    if (!ibCellProcPtsPtr_)
    {
        makeIbCellPts();
    }

    return *ibCellProcPtsPtr_;
}

const Foam::vectorListList&
Foam::immersedBoundaryFvPatch::ibProcPts() const
{
    if (!ibProcPtsPtr_)
    {
        makeIbCellPts();
    }

    return *ibProcPtsPtr_;
}



void Foam::immersedBoundaryFvPatch::updateTriMesh() 
{
    this->movePoints();
}

 
// ************************************************************************* //ppI
