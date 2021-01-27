/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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

Description

\*---------------------------------------------------------------------------*/
#ifndef ibTriSurfaceTools_H
#define ibTriSurfaceTools_H


#include "OFstream.H"
#include "DynamicList.H"
#include "word.H"
#include "triSurfaceFields.H"
#include "triSurface.H"
#include "fvMesh.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes


/*----------------------------------------------------------------------------*\
        Class ibTriSurfaceTools Declaration
\*----------------------------------------------------------------------------*/

class ibTriSurfaceTools
{
   private:
    // Private data
        mutable triSurface triSurface_;
        mutable vector g_;
        mutable scalar repose_;

        mutable bool sandSlidePtr_;
        mutable bool sandSlideoutsidePtr_;

        mutable scalar massRatio_;

        void makeSandSlide(const fvMesh& mesh, const bool&);
 

        void checkSandSlide(const fvMesh& mesh, const bool&);
 

        //- Make tri addressing added by Y. Xu
        void makeTriPointsToHitPoints() const;

        mutable labelListList triPointsToHitPointsAddrPtr_;
      
        mutable scalarListList triPointsToHitPointsWeightsPtr_;

		mutable boolList triInMesh_;

		// - calculate dz for different condition
        void calc_dz
        (
			const fvMesh& mesh,
			vectorField& cf,
            const labelList& pAddr,
            scalarField Ah,
            label triI,
			scalar& dz,
			scalar& totalMass,
			const bool&
        );

        bool inMesh
		(
			const fvMesh& mesh,
			const label& triI
		);

        //added by Y.Xu  
        //- Interpolation addressing from ghost cell hit points 
        //- to tri points
        const labelListList& triPointsToHitPointsAddr() const;

        //- Interpolation weights from ghost cell hit points to tri points
        const scalarListList& triPointsToHitPointsWeights() const;

        template<class Type>
        tmp<Field<Type> > toTriPoints
        (
            const Field<Type>& triFaceCentresValues
        ) const;

        template<class Type>
        tmp<Field<Type> > toTriPoints
        (
            const tmp<const Field<Type>& > ttriFaceCentresValues
        ) const;

        template<class Type>
        tmp<vectorField > gradTriFaces
        (
            const Field<Type>& triFaceCentresValues
        ) const;

        template<class Type>
        tmp<Field<Type> > gradTriFaces
        (
            const tmp<const Field<Type>& > ttriFaceCentresValues
        ) const;
    public:

    // Member Functions

        void sandSlide
        (
            triSurface& surf,
            const fvMesh& mesh,
            vector& gD,
            scalar& reposeD
        );

        void writeTriValue
        (
            const scalarField triangularValues,
            const fvMesh& mesh,
            const word& varName,
            const word& countName
        ) const;

        scalar Hausdorff_dist
        (
            const pointField& pts1,
            const pointField& pts2
        );


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

