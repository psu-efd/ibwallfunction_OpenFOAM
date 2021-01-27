/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

Application
    test inside or outside of each triangle of a STL triSurface 

Description
    by Xiaofeng Liu, PSU, 2016

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "triangle.H"
#include "triSurface.H"
#include "triSurfaceTools.H"
#include "triSurfaceSearch.H"
#include "argList.H"
#include "OFstream.H"
#include "IFstream.H"
#include "surfaceIntersection.H"
#include "SortableList.H"
#include "PatchTools.H"
#include "cellSet.H"
#include "PtrList.H"
#include "HashTable.H"
#include "triSurfaceFields.H"
//#include "itoa.H"
#include "triSurfaceMesh.H"

#include "meshSearch.H"

#include "mathematicalConstants.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //Read in the surface STL file
    fileName surfFileName("constant/triSurface/ibCylinder.stl");
    triSurface surf(surfFileName);

    Pout<< "Statistics:" << endl;
    surf.writeStats(Pout);
    Pout<< endl;
    //End of reading in STL file 

    const vectorField& triCentres = surf.faceCentres();

    meshSearch queryMesh(mesh, polyMesh::FACE_DIAG_TRIS);

    forAll(surf, triI)
    {

       if(queryMesh.findCell(triCentres[triI]) == -1) //outside
       {
            Pout << "triangle number = " << triI << " outside " << endl;
       }
       else
       {
            Pout << "triangle number = " << triI << " inside " <<  endl;
       }
    }

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
