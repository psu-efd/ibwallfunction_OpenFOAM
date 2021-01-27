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

Class
    Foam::immersedBoundarySurfaceToCell

Description
    A topoSetSource to select cells based on relation to immersed
    boundary surface.

    Selects:
    - add ibmRefineMeshInsideLayers of cells to inside
    - and add ibmRefineMeshOutsideLayers of cells to inside

    
SourceFiles
    immersedBoundarySurfaceToCell.C

\*---------------------------------------------------------------------------*/

#ifndef immersedBoundarySurfaceToCell_H
#define immersedBoundarySurfaceToCell_H

#include "topoSetSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
class triSurfaceSearch;
class triSurface;

/*---------------------------------------------------------------------------*\
                           Class immersedBoundarySurfaceToCell Declaration
\*---------------------------------------------------------------------------*/

class immersedBoundarySurfaceToCell
:
    public topoSetSource
{

    // Private data

        //- Add usage string
        static addToUsageTable usage_;

        //- Name of surface file
        const fileName surfName_;

        //- Points which are outside
        const pointField outsidePoints_;

        //- How many layers to add inside and outside
        const label ibmRefineMeshInsideLayers_;
        const label ibmRefineMeshOutsideLayers_;

        //- triSurface to search on. On pointer since can be external.
        const triSurface* surfPtr_;

        //- search engine on surface.
        const triSurfaceSearch* querySurfPtr_;
        
    // Private Member Functions

        //- Depending on surface add to or delete from cellSet.
        void combine(topoSet& set, const bool add) const;

        //- Check values at construction time.
        void checkSettings() const;

        const triSurfaceSearch& querySurf() const
        {
            return *querySurfPtr_;
        }
 

public:

    //- Runtime type information
    TypeName("immersedBoundarySurfaceToCell");

    // Constructors

        //- Construct from components
        immersedBoundarySurfaceToCell
        (
            const polyMesh& mesh,
            const fileName& surfName,
            const pointField& outsidePoints,
            const label ibmRefineMeshInsideLayers,
            const label ibmRefineMeshOutsideLayers
        );

        //- Construct from components (supplied surface, surfaceSearch)
        immersedBoundarySurfaceToCell
        (
            const polyMesh& mesh,
            const fileName& surfName,
            const pointField& outsidePoints,
            const label ibmRefineMeshInsideLayers,
            const label ibmRefineMeshOutsideLayers,
            const triSurface& surf,
            const triSurfaceSearch& querySurf
        );

        //- Construct from dictionary
        immersedBoundarySurfaceToCell
        (
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Construct from Istream
        immersedBoundarySurfaceToCell
        (
            const polyMesh& mesh,
            Istream&
        );


    //- Destructor
    virtual ~immersedBoundarySurfaceToCell();


    // Member Functions

        virtual sourceType setType() const
        {
            return CELLSETSOURCE;
        }

        virtual void applyToSet
        (
            const topoSetSource::setAction action,
            topoSet&
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
