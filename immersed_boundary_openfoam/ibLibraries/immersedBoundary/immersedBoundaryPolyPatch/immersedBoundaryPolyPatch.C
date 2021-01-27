/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "immersedBoundaryPolyPatch.H"
#include "polyBoundaryMesh.H"
#include "polyMesh.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

#include "cellClassification.H"
#include "meshSearch.H"
#include "triSurface.H"
#include "triSurfaceSearch.H"
#include "demandDrivenData.H"

#include "surfaceSets.H"
#include "treeDataFace.H"
#include "cellInfo.H"
#include "globalMeshData.H"
#include "MeshWave.H"

#include "argList.H"

#include "ibTriSurfaceTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(immersedBoundaryPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, immersedBoundaryPolyPatch, word);
    addToRunTimeSelectionTable
    (
        polyPatch,
        immersedBoundaryPolyPatch,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
Foam::scalar Foam::immersedBoundaryPolyPatch::cellSize(label cellID, const fvMesh& mesh_) const
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

// to see if the cut cell center is inside the solid part
bool Foam::immersedBoundaryPolyPatch::ifInside(label cellID, const fvMesh& mesh_) const
{ 
    const triSurfaceSearch& tss = triSurfSearch();

    vector C = mesh_.cellCentres()[cellID];

    scalar Delta = cellSize(cellID,mesh_);
 
    return ifInside(C, Delta, mesh_);
}

// to see if the cut cell center is inside the solid part
bool Foam::immersedBoundaryPolyPatch::ifInside(vector C, scalar Delta, const fvMesh& mesh_) const
{
    bool inside = true;
 
    const triSurfaceSearch& tss = triSurfSearch();

    vector span
        (
            30*Delta,
            30*Delta,
            30*Delta
        );
    pointIndexHit pih = tss.nearest(C, span);

    if (!pih.hit())
    {
        boundBox* meshBB = new boundBox(mesh_.points(), false);

        span = meshBB->span();

        pih = tss.nearest(C, span);
		
		delete meshBB;
    }

    if (pih.hit())
    {       
        point nearestPoint = pih.hitPoint();
        vector nearestPointNormal =
            triSurfaceTools::surfaceNormal
            (
                ibMesh(),
                pih.index(),
                pih.hitPoint()
            );
         // Note: iptsNormals point OUT of the domain
        if (!internalFlow())
        {
            nearestPointNormal *= -1;
        }

        scalar Indicator = nearestPointNormal & ( nearestPoint - C );

        if (Indicator>0)
        {
            inside=false;
        }
    }
    else
    {
        FatalErrorIn
            (
                "Foam::bool Foam::immersedBoundaryPolyPatch::ifInside(vector C) const"
            )   << "Can't find nearest triSurface point for point "
                << C<< ", "
                << "span = " << span
                << "\nYou could try to increase the search span. "
                << abort(FatalError);
    }


    return inside;
}

void Foam::immersedBoundaryPolyPatch::makeTriSurfSearch() const
{   
    if (debug)
    {
        Info<< "void immersedBoundaryPolyPatch::makeTriSurfSearch() const : "
            << "creating triSurface search algorithm"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already
    if (triSurfSearchPtr_)
    {
        FatalErrorIn("immersedBoundaryPolyPatch::makeTriSurfSearch() const")
            << "triSurface search algorithm already exist"
            << abort(FatalError);
    }

    surfPtr_= new triSurface(ibMesh_);

    triSurfSearchPtr_ = new triSurfaceSearch(*surfPtr_);  

	//delete surfPtr_;
//Pout<<" triSurfSearchPtr_ "<<endl;
//triSurfSearchPtr_->surface().writeStats(Pout);
   //triSurfSearchPtr_ = new triSurfaceSearch(*surfPtr_,     
   //                      indexedOctree<treeDataTriSurface>::perturbTol(),
   //                      6);
}


void Foam::immersedBoundaryPolyPatch::clearOut()
{
    deleteDemandDrivenData(triSurfSearchPtr_);
    deleteDemandDrivenData(surfPtr_);
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immersedBoundaryPolyPatch::immersedBoundaryPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    polyPatch(name, size, start, index, bm, patchType),
    ibMesh_
    (
        IOobject
        (
            name  + ".stl",
            bm.mesh().time().constant(), // instance
            "triSurface",                // local
            bm.mesh(),                   // registry
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    internalFlow_(false),
    movingIb_(false),
    triSurfSearchPtr_(NULL),
	surfPtr_(NULL)
{
}


Foam::immersedBoundaryPolyPatch::immersedBoundaryPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    polyPatch(name, dict, index, bm, patchType),
    ibMesh_
    (
        IOobject
        (
            name  + ".stl",
            bm.mesh().time().constant(), // instance
            "triSurface",                // local
            bm.mesh(),                   // registry
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    internalFlow_(dict.lookup("internalFlow")),
    movingIb_(false),
    triSurfSearchPtr_(NULL),
	surfPtr_(NULL)
{
    if (size() > 0)
    {
        FatalIOErrorIn
        (
            "immersedBoundaryPolyPatch::immersedBoundaryPolyPatch\n"
            "(\n"
            "    const word& name,\n"
            "    const dictionary& dict,\n"
            "    const label index,\n"
            "    const polyBoundaryMesh& bm\n"
            ")",
            dict
        )   << "Faces detected in the immersedBoundaryPolyPatch.  "
            << "This is not allowed: please make sure that the patch size "
            << "equals zero."
            << abort(FatalIOError);
    }
}


Foam::immersedBoundaryPolyPatch::immersedBoundaryPolyPatch
(
    const immersedBoundaryPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    polyPatch(pp, bm),
    ibMesh_
    (
        IOobject
        (
            pp.name() + ".stl",
            bm.mesh().time().constant(), // instance
            "triSurface",                // local
            bm.mesh(),                   // registry
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    internalFlow_(pp.internalFlow()),
    movingIb_(false),
    triSurfSearchPtr_(NULL),
	surfPtr_(NULL)
{}


Foam::immersedBoundaryPolyPatch::immersedBoundaryPolyPatch
(
    const immersedBoundaryPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    polyPatch(pp, bm, index, newSize, newStart),
    ibMesh_
    (
        IOobject
        (
            pp.name() + ".stl",
            bm.mesh().time().constant(), // instance
            "triSurface",                // local
            bm.mesh(),                   // registry
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    internalFlow_(pp.internalFlow()),
    movingIb_(false),
    triSurfSearchPtr_(NULL),
	surfPtr_(NULL)
{

}
 




// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

Foam::immersedBoundaryPolyPatch::~immersedBoundaryPolyPatch()
{
    clearOut();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::triSurfaceSearch&
Foam::immersedBoundaryPolyPatch::triSurfSearch() const
{


    if (!triSurfSearchPtr_)
    {

        makeTriSurfSearch();
    }

    return *triSurfSearchPtr_;
}


void Foam::immersedBoundaryPolyPatch::write(Ostream& os) const
{
    polyPatch::write(os);
    os.writeKeyword("internalFlow") << internalFlow_
        << token::END_STATEMENT << nl;
}




// build inside cell lists 
// 0-outside 1-inside 2-cut(with center inside) 3-cut(with center outside)
Foam::labelList Foam::immersedBoundaryPolyPatch::markInside
(
    const fvMesh& mesh,
    const pointField& outsidePts
)
const
{
 //   cpuTime timer;
 //   argList::noParallel();

    meshSearch queryMesh(mesh);
   
    const triSurfaceSearch& querySurf = triSurfSearch();

    labelListList ibProcOutsidePtsCellIndex(Pstream::nProcs());

    labelListList& procOutsidePtsCellIndex = ibProcOutsidePtsCellIndex;

    if (Pstream::parRun())//Start of mpi run
    {
        labelHashSet procOutsidePtsCellIndexSet;

        forAll(outsidePts, outsideI)
        {
            const point& outsidePoint = outsidePts[outsideI];

            label cellIndex = queryMesh.findCell(outsidePoint, -1, false);

            procOutsidePtsCellIndexSet.insert(cellIndex); //cellIndex for each outsidePt in each processor

        }
        procOutsidePtsCellIndex[Pstream::myProcNo()]=procOutsidePtsCellIndexSet.toc();

        Pstream::gatherList(procOutsidePtsCellIndex);
        Pstream::scatterList(procOutsidePtsCellIndex);

        forAll(outsidePts, outsideI)
        {
            const point& outsidePoint = outsidePts[outsideI];
               
            label cellIndex = -1;
            for (label procI = 0; procI < Pstream::nProcs(); procI++)
            {
                if(procOutsidePtsCellIndex[procI][outsideI]!=-1)
                {
                    cellIndex=procOutsidePtsCellIndex[procI][outsideI];  
                    // two processors may share one outside point, but it may not affect the result
                }
            }            
            if (cellIndex < SMALL)
            {
                FatalErrorIn("immersedBoundaryPolyPatch::markInside")
                    << "outsidePoint " << outsidePoint
                    << " is not inside any cell"
                    << exit(FatalError);
            }
        }

    }//End of mpi run
    else//Start of serial run
    {
        forAll(outsidePts, outsideI)
        {
            const point& outsidePoint = outsidePts[outsideI];

            // Find cell that the current point is in. Linear search.

            label cellIndex = queryMesh.findCell(outsidePoint, -1, false);

            if (cellIndex < SMALL)
            {
                FatalErrorIn("immersedBoundaryPolyPatch::markInside")
                    << "outsidePoint " << outsidePoint
                    << " is not inside any cell"
                    << exit(FatalError);
            }
        }
    }//End of serial run


/*
    List<label> cellType(mesh.cells().size());
    #include "markFaces.H"
    #include "markCells.H"
*/
 
    cellClassification cellType
    (
            mesh,
            queryMesh,
            querySurf,
            outsidePts
    );


    // build inside cell lists
    labelHashSet insideCells;
    labelHashSet outsideCells;
    labelHashSet cutCells;
 
    const triSurface& triS_ = ibMesh();

    scalar Delta = mag(outsidePts[0]-average(triS_.points()));

    labelList inside(cellType.size()); 

    scalar insiders=0;
    forAll(cellType, cellI)
    {
            label cType = cellType[cellI];
            inside[cellI]=0;
            if(cType == cellClassification::CUT)
            {
                if(!cutCells.found(cellI)) cutCells.insert(cellI);

                if(ifInside(cellI,mesh)!=ifInside(outsidePts[0],Delta,mesh)){inside[cellI]=2;insiders++;}

                else{inside[cellI]=3;}
            }
            else if(cType == cellClassification::INSIDE)
            {
                if(!insideCells.found(cellI)) insideCells.insert(cellI);inside[cellI]=1;insiders++;
            }
/*            else
            {
                if(!outsideCells.found(cellI)) outsideCells.insert(cellI);inside[cellI]=0;
            }
*/    }

    //cellType.writeStats(Pout);
    return inside;
}


void Foam::immersedBoundaryPolyPatch::moveTriSurfacePoints(const pointField& p) 
{
 // Record the motion of the patch
    movingIb_ = true;

    // Move points of the triSurface
    const pointField& oldPoints = ibMesh_.points();

    if (oldPoints.size() != p.size())
    {
        FatalErrorIn
        (
            "void immersedBoundaryPolyPatch::moveTriSurfacePoints\n"
            "(\n"
            "    const pointField& p\n"
            ")"
        )   << "Incorrect size of motion points for patch " << name()
            << ".  oldPoints = "
            << oldPoints.size() << " p = " << p.size()
            << abort(FatalError);
    }

    Info<< "Moving immersed boundary points for patch " << name()
        << endl;

/*    if(Pstream::parRun()) //If runnning parallel
    {
        typedef List<pointField> pointFieldList;
        pointFieldList procNewPoints(Pstream::nProcs());
        procNewPoints[Pstream::myProcNo()] = p;

        Pstream::gatherList(procNewPoints);
        Pstream::scatterList(procNewPoints);
        forAll(procNewPoints, procI)
        {
            ibMesh_.movePoints(procNewPoints[procI]);
        }
    }
    else
    {
        // Move points
        ibMesh_.movePoints(p);
    }
*/

    ibMesh_.movePoints(p);
    fileName path;
/*    if (Pstream::parRun())
    {
        path = boundaryMesh().mesh().time().path()/".."/"postProcessing"/"VTK";
    }
    else
    {
        path = boundaryMesh().mesh().time().path()/"postProcessing"/"VTK";
    }

        
    if(Pstream::myProcNo()==Pstream::masterNo())
    {
        mkDir(path);
        ibMesh_.triSurface::write
        (
            path/
            word
            (
                name() + "_"
              + Foam::name(boundaryMesh().mesh().time().timeName())
              + ".stl"
            )
        );
    }

    if(Pstream::myProcNo()==Pstream::masterNo())
    {
        fileName path(boundaryMesh().mesh().time().constant()/"triSurface"/name()+".stl");
        ibMesh_.triSurface::write(path);
        Info<<"write stl to "<<path<<endl;
    }*/

    clearOut();

}


void Foam::immersedBoundaryPolyPatch::sandSlide(const fvMesh& mesh_)
{
    IOdictionary ibmDict
    (
        IOobject
        (
            "ibmDict",
            mesh_.time().system(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // read constants from ibmDict
    vector gD(ibmDict.lookup("gDirection"));
    scalar reposeD(readScalar(ibmDict.lookup("staticFrictionCoefficient")));
 
    triSurface& triS_(ibMesh_);
    ibTriSurfaceTools ibT;
    ibT.sandSlide(triS_,mesh_,gD,reposeD);
}


void Foam::immersedBoundaryPolyPatch::operator=(const immersedBoundaryPolyPatch& p)
{
    clearAddressing();

    patchIdentifier::operator=(p);
    primitivePatch::operator=(p);
    //start_ = p.start_;
}


void Foam::immersedBoundaryPolyPatch::deleteAll() 
{
    clearOut();
    deleteDemandDrivenData(surfPtr_);
}


void Foam::immersedBoundaryPolyPatch::moveIb()
{
    movingIb_ = true;
}
// ************************************************************************* //
