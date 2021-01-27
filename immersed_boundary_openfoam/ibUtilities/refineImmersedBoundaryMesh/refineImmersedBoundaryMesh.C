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

Application
    refineImmersedBoundaryMesh

Description
    Refine the background mesh around the immersed surface

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "refineImmersedBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "refine cells in multiple directions"
    );

    #include "addOverwriteOption.H"
    #include "addRegionOption.H"
    #include "addDictOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    runTime.functionObjects().off();
    #include "createMesh.H"
//    #include "createNamedPolyMesh.H"
    const word oldInstance = mesh.pointsInstance();


    const bool overwrite = args.optionFound("overwrite");

    refineImmersedBoundaryMesh rib(mesh);

    Info<< "Number of refinement cells = " << rib.refinementCells().size()
        << endl;

    string oldTimeName(runTime.timeName());

    if (!overwrite)
    {
        runTime++;
    }

    rib.refineMesh(rib.refinementCells());

    // Write resulting mesh
    if (overwrite)
    {
        mesh.setInstance(oldInstance);
    }
    mesh.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
