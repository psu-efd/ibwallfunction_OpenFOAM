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
    test the motion of a STL triSurface 

Description
    by Xiaofeng Liu, PSU, 2014

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

#include "triSurfaceFieldToTecplot.H"

#include "mathematicalConstants.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //Read in the surface STL file
    fileName surfFileName("constant/triSurface/"+stlFile);
    triSurface surf(surfFileName);

    Pout<< "Statistics:" << endl;
    surf.writeStats(Pout);
    Pout<< endl;
    //End of reading in STL file 

    fileName outPathName(runTime.path()/"triSurface");
    mkDir(outPathName);

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
        
        //translation
        vector translation(vector::zero);
    
        //rotation - not done yet
        
        //set up some fake translation 
        translation.x() = 0.02*
           ::sin(constant::mathematical::pi*4.0*runTime.time().value());
        Info << "translation = " << translation << endl;

        pointField points(surf.points());
        
        points += translation;

        surf.movePoints(points);

        if(runTime.outputTime())
        {
/*
           fileName outPathName(runTime.timePath()/"triSurface");
           mkDir(outPathName);

           fileName outFileName(outPathName/stlFile);
           Info << outFileName << endl;
           surf.write(outFileName);
*/

           fileName outFileName(outPathName/"stl_"
                               +name(runTime.timeIndex())+".vtk");
           Info << outFileName << endl;
           surf.write(outFileName);

           runTime.write();
        }
    }    

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
