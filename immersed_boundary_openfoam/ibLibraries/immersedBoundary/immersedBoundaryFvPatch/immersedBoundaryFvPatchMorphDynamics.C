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
    
    immersedBoundaryFvPatchMorphDynamics.C
\*---------------------------------------------------------------------------*/

#include "immersedBoundaryFvPatch.H"
#include "fvMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * Puiblic Member Functions  * * * * * * * * * * * //
/*
void Foam::immersedBoundaryFvPatch::sandSlide() const
{
    makeConstants();
    triSurface& triS_(ibMeshRef());
    ibTriSurfaceTools ibT;
    ibT.sandSlide(triS_,mesh_,g_,repose_);
}

//movePoints  added by Y. Xu
void Foam::immersedBoundaryFvPatch::updatePoints( pointField& pts) const
{
    triSurface& triS_(ibMeshRef());    
    if (pts.size() != triS_.points().size())
    {
        FatalErrorIn
        (
            "void Foam::immersedBoundaryFvPatch::movePoints\n"
            "(\n"
            "    pointField& pts\n"
            ") const"
        )   << "Field size does not correspond to size of immersed boundary "
            << "triangulated surface for patch " << name() << nl
            << "Field size = " << pts.size()
            << " surface size = " << triS_.points().size()
            << abort(FatalError);
    }
    ibPolyPatchRef_.updateTriSurface(pts);
    //triS_.movePoints(pts);

}

void Foam::immersedBoundaryFvPatch::makeConstants() const
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

    g_=gD;
    repose_=reposeD;
}
*/
// ************************************************************************* //
