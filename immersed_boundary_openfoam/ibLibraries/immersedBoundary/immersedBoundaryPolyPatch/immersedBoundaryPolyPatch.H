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

Class
    Foam::immersedBoundaryPolyPatch

Description
    Immersed boundary patch

Author
    Zeljko Tukovic
    Reorganisation by Hrvoje Jasak

SourceFiles
    immersedBoundaryPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef immersedBoundaryPolyPatch_H
#define immersedBoundaryPolyPatch_H

#include "polyPatch.H"
#include "triSurfaceMesh.H"
#include "triSurfaceTools.H"
#include "triSurfaceSearch.H"
#include "Switch.H"
#include "boolList.H"
#include "volFields.H"
#include "topoSetSource.H"

#include "fvMesh.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class immersedBoundaryPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class immersedBoundaryPolyPatch
:
    public polyPatch
{
    // Private data

        //- Triangular surface representing immersed boundary.
        //  Name of tri surface will be identical to the name of the patch
       mutable triSurfaceMesh ibMesh_;

        //- Internal or external flow calculation
        Switch internalFlow_;


    // Demand-driven data

        //- Moving immersed boundary
        bool movingIb_;

        //- Triangular surface search algorithm
        mutable triSurfaceSearch* triSurfSearchPtr_;


    // Private Member Functions

        //- Find the size of a cell
        scalar cellSize(label cellID, const fvMesh& mesh) const;

        // to see if the cut cell center is inside the fluid part
        bool ifInside(label cellI, const fvMesh& mesh) const;
        bool ifInside(vector C, scalar Delta, const fvMesh& mesh) const;

        // Storage management

            //- Clear all demand-driven data
            void clearOut();


        // Make demand-driven data

            //- Make triSurface search algorithm
            void makeTriSurfSearch() const;

        //- triSurface to search on. On pointer since can be external.
        mutable triSurface* surfPtr_;

public:

    //- Runtime type information
    TypeName("immersedBoundary");


    // Constructors

        //- Construct from components
        immersedBoundaryPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct from dictionary
        immersedBoundaryPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct as copy, resetting the boundary mesh
        immersedBoundaryPolyPatch
        (
            const immersedBoundaryPolyPatch&,
            const polyBoundaryMesh&
        );

        //- Construct given the original patch and resetting the
        //  face list and boundary mesh information
        immersedBoundaryPolyPatch
        (
            const immersedBoundaryPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        );
        //- Construct as copy
        //immersedBoundaryPolyPatch(const polyPatch&);
        immersedBoundaryPolyPatch(polyPatch);
        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>
            (
                new immersedBoundaryPolyPatch(*this, bm)
            );
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        ) const
        {
            return autoPtr<polyPatch>
            (
                new immersedBoundaryPolyPatch
                (
                    *this,
                    bm,
                    index,
                    newSize,
                    newStart
                )
            );
        }

 



        
        //XL: added for whoever needs to call clone
        //    such as decomposePar
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const label newStart
        ) const
        {
            return autoPtr<polyPatch>
            (
                new immersedBoundaryPolyPatch
                (
                    *this, 
                    bm, 
                    index, 
                    0, 
                    newStart
                )
            );
        }



    //- Destructor
    virtual ~immersedBoundaryPolyPatch();


    // Member Functions

        // Access

            //- Return immersed boundary surface mesh
            const triSurfaceMesh& ibMesh() const
            {
                return ibMesh_;
            }
            triSurfaceMesh& ibMeshRef() 
            {
                return ibMesh_;
            }
            //- Return true if solving for flow inside the immersed boundary
            bool internalFlow() const
            {
                return internalFlow_;
            }

            //- Return triSurface search object
            const triSurfaceSearch& triSurfSearch() const;

            //- Return true if immersed boundary is moving
            bool movingIb() const
            {
                return movingIb_;
            }            

            //- Return true if immersed boundary is moving
            void moveIb();   

            //- mark inside cell list added by Y. Xu 0-outside 1-inside 2-cut(with center inside) 3-cut(with center outside)
            labelList markInside
            (
                 const fvMesh& mesh_,
                 const pointField& outsidePoints_
            ) const;

        //- Write
        virtual void write(Ostream&) const;

        // Edit

            //- Correct patches after moving points
            void moveTriSurfacePoints(const pointField& p);
  

        void operator=(const immersedBoundaryPolyPatch&);

        void deleteAll();

        void sandSlide(const fvMesh& mesh_);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
