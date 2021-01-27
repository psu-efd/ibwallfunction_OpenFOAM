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
    Foam::immersedBoundaryFixedFluxPressure

Description
    Foam::immersedBoundaryFixedFluxPressure

    No need for this B.C. for body force type of IBM

SourceFiles
    immersedBoundaryFixedFluxPressure.C

\*---------------------------------------------------------------------------*/

#ifndef immersedBoundaryFixedFluxPressure_H
#define immersedBoundaryFixedFluxPressure_H

#include "fvPatchFields.H"
#include "immersedBoundaryFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class immersedBoundaryFixedFluxPressure Declaration
\*---------------------------------------------------------------------------*/

class immersedBoundaryFixedFluxPressure
:
    public fvPatchScalarField
{
    // Private data

        //- Local reference cast into the immersedBoundaryFvPatch
        const immersedBoundaryFvPatch& ibPatch_;

        //- Local reference to fvMesh
        const fvMesh& mesh_;

        // Defining fields
        // Note: defining fields on faces of the IB patch
        // represented as a triangulated surface

            //- Defining value field
            mutable Field<scalar> refValue_;

            //- Defining normal gradient field
            mutable Field<scalar> refGrad_;


        //- Name of the predicted flux transporting the field
        word phiHbyAName_;

        //- Name of the flux transporting the field
        word phiName_;

        //- Name of the density field used to normalise the mass flux
        //  if neccessary
        word rhoName_;

        //- Name of the pressure diffusivity field
        word DpName_;

        //- Is the pressure adjoint, i.e. has the opposite sign
        Switch adjoint_;


        // Dead cell controls

            //- Set dead (fully inside) cell value
            mutable Switch setDeadCellValue_;

            //- Dead (fully inside) cell value
            mutable scalar deadCellValue_;


        // Fields on the ghost cells' nearest points to the triSurface
        // Field size equals the number of ghost cells

            //- Value field
            mutable Field<scalar> ibValue_;

            //- Normal gradient field
            mutable Field<scalar> ibGrad_;

        // Fields on the immersed boundary ghost cells
        // Field size equals the number of ghost cells

            //- Does the boundary condition has any ghost cells at all?
            mutable bool hasGhostCells_;

            //- Value field
            mutable Field<scalar> ibGhostCellValues_;




    // Private Member Functions

        //- Update the ghost cell values
        void updateGhostCellValues() const;

        //- Impose condition in ghost cells
        void imposeGhostCellsCondition() const;
 
        //- Impose condition in fully inside (dead) cells
        void imposeFullyInsideCellsCondition() const;
 
public:

    //- Runtime type information
    TypeName("immersedBoundaryFixedFluxPressure");


    // Constructors

        //- Construct from patch and internal field
        immersedBoundaryFixedFluxPressure
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        immersedBoundaryFixedFluxPressure
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given immersedBoundaryFixedFluxPressure
        //  onto a new patch
        immersedBoundaryFixedFluxPressure
        (
            const immersedBoundaryFixedFluxPressure&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        immersedBoundaryFixedFluxPressure
        (
            const immersedBoundaryFixedFluxPressure&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new immersedBoundaryFixedFluxPressure(*this)
            );
        }

        //- Construct as copy setting internal field reference
        immersedBoundaryFixedFluxPressure
        (
            const immersedBoundaryFixedFluxPressure&,
            const DimensionedField<Scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new immersedBoundaryFixedFluxPressure(*this, iF)
            );
        }


    // Member functions

        //- Return reference to immersed boundary patch
        const immersedBoundaryFvPatch& ibPatch() const
        {
            return ibPatch_;
        }


        // Return defining fields
        // Note: defining fields carry values on faces of the IB patch
        // represented as a triangulated surface

            //- Return access to reference value
            virtual Field<scalar>& refValue()
            {
                return refValue_;
            }

            //- Return reference value
            virtual const Field<scalar>& refValue() const
            {
                return refValue_;
            }

            //- Return access to reference gradient
            virtual Field<scalar>& refGrad()
            {
                return refGrad_;
            }

            //- Return reference gradient
            virtual const Field<scalar>& refGrad() const
            {
                return refGrad_;
            }

            //- Return true if this patch field fixes a value
            virtual bool fixesValue() const
            {
                return fixesValue_;
            }

            //- Return access to fixes value switch
            Switch& fixesValue()
            {
                return fixesValue_;
            }


        //- Return fields on ghost cells' nearest points to the IB

            //- Return IB field
            const Field<scalar>& ibValue() const;

            //- Return IB surface-normal gradient
            const Field<scalar>& ibGrad() const;

            //- Return ghost cell values
            const Field<scalar>& ibGhostCellValues() const;


        //- Return fields on triangular faces

            //- Return IB field
            tmp<Field<scalar> > triValue() const;

            //- Return IB surface-normal gradient
            tmp<Field<scalar> > triGrad() const;

            void writeTriValue(const word& varName) const;

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            )
            {}

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<scalar>&,
                const labelList&
            )
            {}


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            void updateCoeffs();

            //- Initialise the evaluation of the patch field
            virtual void initEvaluate
            (
                const Pstream::commsTypes commsType = Pstream::blocking
            );

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType = Pstream::blocking
            );

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<scalar> > valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const
            {
                return tmp<Field<scalar> >(new Field<scalar>(0));
            }

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<scalar> > valueBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const
            {
                return tmp<Field<scalar> >(new Field<scalar>(0));
            }

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            tmp<Field<scalar> > gradientInternalCoeffs() const
            {
                return tmp<Field<scalar> >(new Field<scalar>(0));
            }

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            tmp<Field<scalar> > gradientBoundaryCoeffs() const
            {
                return tmp<Field<scalar> >(new Field<scalar>(0));
            }

            //- impose the immersed boundary effect  
            void imposeImmersedBoundaryEffect() const;

            void copyBCParameters(const immersedBoundaryFixedFluxPressure& ptf) const; 


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "immersedBoundaryFixedFluxPressure.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
