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
    Foam::immersedBoundaryNutkWallFunctionFvPatchScalarField

Group
    grpWallFunctions

Description
    This boundary condition provides a turbulent kinematic viscosity condition
    when using wall functions, based on turbulence kinetic energy.
    - replicates OpenFOAM v1.5 (and earlier) behaviour

    \heading Patch usage

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            immersedBoundaryNutkWallFunction;
    }
    \endverbatim

SeeAlso
    Foam::nutWallFunctionFvPatchScalarField

SourceFiles
    immersedBoundaryNutkWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef immersedBoundaryNutkWallFunctionFvPatchScalarField_H
#define immersedBoundaryNutkWallFunctionFvPatchScalarField_H

#include "immersedBoundaryWallFunctionFvPatchFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class immersedBoundaryNutkWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class immersedBoundaryNutkWallFunctionFvPatchScalarField
:
    public immersedBoundaryWallFunctionFvPatchField<scalar>
{
private:

        //- Calculate the turbulence viscosity
        virtual tmp<scalarField> calcNut() const;


        //- Cmu coefficient
        mutable  scalar Cmu_;

        //- Von Karman constant
        scalar kappa_;




        //- E coefficient
        scalar E_;


        //- Y+ at the edge of the laminar sublayer
        scalar yPlusLam_;
protected:

    // Protected Member Functions

public:

    //- Runtime type information
    TypeName("immersedBoundaryNutkWallFunction");


    // Constructors

        //- Construct from patch and internal field
        immersedBoundaryNutkWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        immersedBoundaryNutkWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  immersedBoundaryNutkWallFunctionFvPatchScalarField
        //  onto a new patch
        immersedBoundaryNutkWallFunctionFvPatchScalarField
        (
            const immersedBoundaryNutkWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        immersedBoundaryNutkWallFunctionFvPatchScalarField
        (
            const immersedBoundaryNutkWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new immersedBoundaryNutkWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        immersedBoundaryNutkWallFunctionFvPatchScalarField
        (
            const immersedBoundaryNutkWallFunctionFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new immersedBoundaryNutkWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

            //- Calculate the Y+ at the edge of the laminar sublayer
            static scalar yPlusLam(const scalar kappa, const scalar E);

       // Evaluation functions
            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Evaluate the patchField
            virtual void evaluate(const Pstream::commsTypes);

            //- Calculate and return the yPlus at the boundary
            //virtual tmp<scalarField> yPlus() const;


        // I-O

            //- Write
            virtual void write(Ostream& os) const;
        //- Calculate and return the yPlus at the boundary
        //virtual tmp<scalarField> yPlus() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
