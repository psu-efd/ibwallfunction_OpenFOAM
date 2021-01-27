/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "immersedBoundarySurfaceToCell.H"
#include "polyMesh.H"
#include "meshSearch.H"
#include "triSurface.H"
#include "triSurfaceSearch.H"
#include "cellClassification.H"
#include "cpuTime.H"
#include "demandDrivenData.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(immersedBoundarySurfaceToCell, 0);

addToRunTimeSelectionTable(topoSetSource, immersedBoundarySurfaceToCell, word);

addToRunTimeSelectionTable(topoSetSource, immersedBoundarySurfaceToCell, istream);

}


Foam::topoSetSource::addToUsageTable Foam::immersedBoundarySurfaceToCell::usage_
(
    immersedBoundarySurfaceToCell::typeName,
    "\n    Usage: immersedBoundarySurfaceToCell"
    "<surface> <outsidePoints> <cut> <inside> <outside> <near> <curvature>\n\n"
    "    <surface> name of triSurface\n"
    "    <outsidePoints> list of points that define outside\n"
    "    <cut> boolean whether to include cells cut by surface\n"
    "    <inside>   ,,                 ,,       inside surface\n"
    "    <outside>  ,,                 ,,       outside surface\n"
    "    <near> scalar; include cells with centre <= near to surface\n"
    "    <curvature> scalar; include cells close to strong curvature"
    " on surface\n"
    "    (curvature defined as difference in surface normal at nearest"
    " point on surface for each vertex of cell)\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::immersedBoundarySurfaceToCell::combine(topoSet& set, const bool add) const
{
    cpuTime timer;

    //
    // Cut cells with surface and classify cells
    //


    // Construct search engine on mesh

    const meshSearch queryMesh(mesh_);


    // Check all 'outside' points
    forAll(outsidePoints_, outsideI)
    {
            const point& outsidePoint = outsidePoints_[outsideI];

            // Find cell that the current point is in. Linear search.
            label cellI = queryMesh.findCell(outsidePoint, -1, false);
            if (returnReduce(cellI, maxOp<label>()) == -1)
            {
                FatalErrorIn("immersedBoundarySurfaceToCell::combine(topoSet&, const bool)")
                    << "outsidePoint " << outsidePoint
                    << " is not inside any cell"
                    << exit(FatalError);
            }
    }

    // Cut faces with surface and classify cells

    cellClassification cellType
    (
            mesh_,
            queryMesh,
            querySurf(),
            outsidePoints_
    );

    Info<< "    Marked inside/outside using surface intersection in = "
        << timer.cpuTimeIncrement() << " s" << endl << endl;

    // build inside, cut, and outside cell lists
    labelHashSet insideCells;
    labelHashSet outsideCells;
    labelHashSet cutCells;
 
    forAll(cellType, cellI)
    {
           label cType = cellType[cellI];
 
            if(cType == cellClassification::CUT)
            {
                if(!cutCells.found(cellI)) cutCells.insert(cellI);
            }
            else if(cType == cellClassification::INSIDE)
            {
                if(!insideCells.found(cellI)) insideCells.insert(cellI);
            }
            else
            {
                if(!outsideCells.found(cellI)) outsideCells.insert(cellI);
            }
    }
       
    labelHashSet oldCutCells(cutCells);

    ////////////////////////////////////////////////////////////////////
    //At this point, oldCutCells and cutCells have the same content.
    ////////////////////////////////////////////////////////////////////

    //add layers of cells inside
    for(label i=0; i<ibmRefineMeshInsideLayers_; i++)
    {
                forAllConstIter(labelHashSet, oldCutCells, cellI)
                {
                    label curCutCell = cellI.key();

                    forAll(mesh_.cellCells()[curCutCell], cellCellsI)
                    {
                        label curCellCells = mesh_.cellCells()[curCutCell][cellCellsI];
                        if(insideCells.found(curCellCells) &&  //it is inside and
                           !cutCells.found(curCellCells)     //not in cut cell list yet
                          )
                        {
                           cutCells.insert(curCellCells);
                           break;
                        }
                    }
                }

               /////////////////////////////////////////////////////////////////////
               //At this point, oldCutCells and cutCells have different content.
               //cutCells should be >= oldCutCells.
               /////////////////////////////////////////////////////////////////////

               //Now, assign cutCells to oldCutCells so they are equal and 
               //get ready for next round of cell addition.
               oldCutCells.clear();
               oldCutCells=cutCells;
    }


    //add layers of cells outside
    for(label i=0; i<ibmRefineMeshOutsideLayers_; i++)
    {
               forAllConstIter(labelHashSet, oldCutCells, cellI)
               {
                  label curCutCell = cellI.key();

                  forAll(mesh_.cellCells()[curCutCell], cellCellsI)
                  {
                      label curCellCells = mesh_.cellCells()[curCutCell][cellCellsI];
                      if(outsideCells.found(curCellCells) &&  //it is outside and
                         !cutCells.found(curCellCells)     //not in cut cell list yet
                        )
                      {
                          cutCells.insert(curCellCells);
                          break;
                      }
                 }
              }

              /////////////////////////////////////////////////////////////////////
              //At this point, oldCutCells and cutCells have different content.
              //cutCells should be >= oldCutCells.
              /////////////////////////////////////////////////////////////////////

             //Now, assign cutCells to oldCutCells so they are equal and
             //get ready for next round of cell addition.
             oldCutCells.clear();
             oldCutCells=cutCells;
    }


    //- Add cells using cutCells set which contains all
    //  the cells which need to be refined
    forAllConstIter(labelHashSet, cutCells, cellI)
    {
            label curCutCell = cellI.key();

            addOrDelete(set, curCutCell, add);
    }
}


void Foam::immersedBoundarySurfaceToCell::checkSettings() const
{
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::immersedBoundarySurfaceToCell::immersedBoundarySurfaceToCell
(
    const polyMesh& mesh,
    const fileName& surfName,
    const pointField& outsidePoints,
    const label ibmRefineMeshInsideLayers,
    const label ibmRefineMeshOutsideLayers
)
:
    topoSetSource(mesh),
    surfName_(surfName),
    outsidePoints_(outsidePoints),
    ibmRefineMeshInsideLayers_(ibmRefineMeshInsideLayers),
    ibmRefineMeshOutsideLayers_(ibmRefineMeshOutsideLayers),
    surfPtr_(new triSurface(surfName_)),
    querySurfPtr_(new triSurfaceSearch(*surfPtr_))
{
    checkSettings();
}


// Construct from components. Externally supplied surface.
Foam::immersedBoundarySurfaceToCell::immersedBoundarySurfaceToCell
(
    const polyMesh& mesh,
    const fileName& surfName,
    const pointField& outsidePoints,
    const label ibmRefineMeshInsideLayers,
    const label ibmRefineMeshOutsideLayers,
    const triSurface& surf,
    const triSurfaceSearch& querySurf
)
:
    topoSetSource(mesh),
    surfName_(surfName),
    outsidePoints_(outsidePoints),
    ibmRefineMeshInsideLayers_(ibmRefineMeshInsideLayers),
    ibmRefineMeshOutsideLayers_(ibmRefineMeshOutsideLayers),
    surfPtr_(&surf),
    querySurfPtr_(&querySurf)
{
    checkSettings();
}


// Construct from dictionary
Foam::immersedBoundarySurfaceToCell::immersedBoundarySurfaceToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    surfName_(dict.lookup("file")),
    outsidePoints_(dict.lookup("outsidePoints")),
    ibmRefineMeshInsideLayers_(readLabel(dict.lookup("ibmRefineMeshInsideLayers"))),
    ibmRefineMeshOutsideLayers_(readLabel(dict.lookup("ibmRefineMeshOutsideLayers"))),
    surfPtr_(new triSurface(surfName_)),
    querySurfPtr_(new triSurfaceSearch(*surfPtr_))
{
    checkSettings();
}


// Construct from Istream
Foam::immersedBoundarySurfaceToCell::immersedBoundarySurfaceToCell
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetSource(mesh),
    surfName_(checkIs(is)),
    outsidePoints_(checkIs(is)),
    ibmRefineMeshInsideLayers_(readLabel(checkIs(is))),
    ibmRefineMeshOutsideLayers_(readLabel(checkIs(is))),
    surfPtr_(new triSurface(surfName_)),
    querySurfPtr_(new triSurfaceSearch(*surfPtr_))
{
    checkSettings();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::immersedBoundarySurfaceToCell::~immersedBoundarySurfaceToCell()
{
    deleteDemandDrivenData(surfPtr_);
    deleteDemandDrivenData(querySurfPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::immersedBoundarySurfaceToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW))
    {
        Info<< "    Adding cells in relation to ibm surface " << surfName_
            << " ..." << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE || (action == topoSetSource::ADD))
    {
        FatalErrorIn
        (
            "immersedBoundarySurfaceToCell:applyToSet()"
        )   << "Illegal action specification. Olny new allowed."
            << exit(FatalError);
    }
}


// ************************************************************************* //
