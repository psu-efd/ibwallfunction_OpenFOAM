/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::immersedBoundaryKqRWallFunctionFvPatchField

Group
    grpWallFunctions

Description
    This boundary condition provides a suitable condition for turbulence
    \c k, \c q, and \c R fields for the case of high Reynolds number flow using
    wall functions.

    It is a simple wrapper around the zero-gradient condition.

    \heading Patch usage

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            immersedBoundaryKqRWallFunction;
    }
    \endverbatim

SeeAlso
    Foam::zeroGradientFvPatchField

SourceFiles
    immersedBoundaryKqRWallFunctionFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef immersedBoundaryKqRWallFunctionFvPatchField_H
#define immersedBoundaryKqRWallFunctionFvPatchField_H

#include "immersedBoundaryWallFunctionFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class immersedBoundaryKqRWallFunctionFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class immersedBoundaryKqRWallFunctionFvPatchField
:
    public immersedBoundaryWallFunctionFvPatchField<Type>
{

    // Private Member Functions

        //- Check the type of the patch
        void checkType();


public:

    //- Runtime type information
    TypeName("immersedBoundaryKqRWallFunction");


    // Constructors

        //- Construct from patch and internal field
        immersedBoundaryKqRWallFunctionFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        immersedBoundaryKqRWallFunctionFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        // immersedBoundaryKqRWallFunctionFvPatchField
        //  onto a new patch
        immersedBoundaryKqRWallFunctionFvPatchField
        (
            const immersedBoundaryKqRWallFunctionFvPatchField&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        immersedBoundaryKqRWallFunctionFvPatchField
        (
            const immersedBoundaryKqRWallFunctionFvPatchField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new immersedBoundaryKqRWallFunctionFvPatchField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        immersedBoundaryKqRWallFunctionFvPatchField
        (
            const immersedBoundaryKqRWallFunctionFvPatchField&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new immersedBoundaryKqRWallFunctionFvPatchField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Evaluate the patchField
            virtual void evaluate
            (
                const Pstream::commsTypes commsType=Pstream::Pstream::blocking
            );


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "immersedBoundaryKqRWallFunctionFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
