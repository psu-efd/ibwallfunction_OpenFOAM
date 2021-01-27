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
    Foam::immersedBoundaryForces

Description
    Calculates the forces and moments by integrating the pressure and
    skin-friction forces over a given list of immersed boundary patches.

    Member function calcForcesMoment() calculates and returns the forces and
    moments on immersed boundary patches.

Author
    Hrvoje Jasak.  All rights reserved.

Modified  
    Only one element in the patchSet_ is allowed, i.e., in each patchSet, there
    is only one STL file (which can contain multiple objects).
 
    Xiaofeng Liu
    Penn State University

SourceFiles
    immersedBoundaryForces.C
    IOimmersedBoundaryForces.H

\*---------------------------------------------------------------------------*/

#ifndef immersedBoundaryForces_H
#define immersedBoundaryForces_H

#include "forces.H"
#include "volFieldsFwd.H"
#include "HashSet.H"
#include "Tuple2.H"
#include "OFstream.H"
#include "Switch.H"
#include "pointFieldFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class immersedBoundaryForces Declaration
\*---------------------------------------------------------------------------*/

class immersedBoundaryForces
:
    public forces
{
public: 

   // Tuple for pressure (.first()) and viscous (.second()) forces
   typedef Tuple2<vector, vector> pressureViscous;

   // Tuple for forces (.first()) and moment (.second())
   // pressure/viscous forces Tuples.
   typedef Tuple2<pressureViscous, pressureViscous> forcesMoments;

protected:

    // Private data
        //- Forces/moment file ptr
        mutable PtrList<OFstream>* filePtrs_;

        //- List of forcesMoment for all object in the STL file
        mutable PtrList<forcesMoments>* fmAllPtrs_;

    // Private Member Functions

        //- If the forces file has not been created create it
        virtual void makeFile();

        //- Allocate memory for fmPtrs_
        void makeFMList() const;

public:

    //- Runtime type information
    TypeName("immersedBoundaryForces");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        immersedBoundaryForces
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~immersedBoundaryForces();


    // Member Functions
        //- Write the forces
        virtual void write();

        //- Calculate and return list of forcesMoment for 
        //  all objects in the STL file
        virtual void calcForcesMomentAll() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
