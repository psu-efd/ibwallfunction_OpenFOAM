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

\*---------------------------------------------------------------------------*/

#include "immersedBoundaryMappedFixedValueFvPatchField.H"
#include "mappedPatchBase.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
immersedBoundaryMappedFixedValueFvPatchField<Type>::immersedBoundaryMappedFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    mappedPatchFieldBase<Type>(this->mapper(p, iF), *this)
{}


template<class Type>
immersedBoundaryMappedFixedValueFvPatchField<Type>::immersedBoundaryMappedFixedValueFvPatchField
(
    const immersedBoundaryMappedFixedValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    mappedPatchFieldBase<Type>(this->mapper(p, iF), *this, ptf)
{}


template<class Type>
immersedBoundaryMappedFixedValueFvPatchField<Type>::immersedBoundaryMappedFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict),
    mappedPatchFieldBase<Type>(this->mapper(p, iF), *this, dict)
{}


template<class Type>
immersedBoundaryMappedFixedValueFvPatchField<Type>::immersedBoundaryMappedFixedValueFvPatchField
(
    const immersedBoundaryMappedFixedValueFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    mappedPatchFieldBase<Type>(ptf)
{}


template<class Type>
immersedBoundaryMappedFixedValueFvPatchField<Type>::immersedBoundaryMappedFixedValueFvPatchField
(
    const immersedBoundaryMappedFixedValueFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    mappedPatchFieldBase<Type>(this->mapper(this->patch(), iF), *this, ptf)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const mappedPatchBase& immersedBoundaryMappedFixedValueFvPatchField<Type>::mapper
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
{
    if (!isA<mappedPatchBase>(p.patch()))
    {
        FatalErrorIn
        (
            "immersedBoundaryMappedFixedValueFvPatchField<Type>::mapper()"
        )   << "\n    patch type '" << p.patch().type()
            << "' not type '" << mappedPatchBase::typeName << "'"
            << "\n    for patch " << p.patch().name()
            << " of field " << iF.name()
            << " in file " << iF.objectPath()
            << exit(FatalError);
    }
    return refCast<const mappedPatchBase>(p.patch());
}


template<class Type>
void immersedBoundaryMappedFixedValueFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    //bool setAverage=this->setAverage_;
    const Type average_=this->average_;
//Info <<"average_ "<<average_<<endl;
   // this->operator==(this->mappedField());

    tmp<Field<Type> > tnewValues(new Field<Type>(0));
    Field<Type>& newValues = tnewValues();
    newValues = this->mappedField();

    const fvPatchField<Type> patchField_=this->patchField_;
    
    const fvPatch& patch_=fixedValueFvPatchField<Type>::patch();

   // const fvPatchScalarField& cellIbMaskExt = patch_.lookupPatchField<volScalarField, scalar>("cellIbMaskExt");

    const fvPatchScalarField& cellIbMask = patch_.lookupPatchField<volScalarField, scalar>("cellIbMaskExt");

    Type averagePsi =pTraits<Type>::zero;

    if(gSum(patchField_.patch().magSf())!=0)
    {

        averagePsi =
            gSum(patchField_.patch().magSf()*newValues*cellIbMask)
           /gSum(patchField_.patch().magSf()*cellIbMask);
    }
//Info <<"gSum(patchField_.patch().magSf()*newValues*cellIbMaskExt) "<<gSum(patchField_.patch().magSf()*newValues*cellIbMaskExt)<<endl;
    if (average_!=pTraits<Type>::zero && mag(averagePsi)/mag(average_) > 0.5)    
    {
            newValues *= mag(average_)/mag(averagePsi)*cellIbMask;
    }
    else
    {
        newValues += (average_ - averagePsi)*cellIbMask;
    }


    this->operator==(tnewValues);    

    if (debug)
    {
        Info<< "mapped on field:"
            << this->dimensionedInternalField().name()
            << " patch:" << this->patch().name()
            << "  avg:" << gAverage(*this)
            << "  min:" << gMin(*this)
            << "  max:" << gMax(*this)
            << endl;
    }

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void immersedBoundaryMappedFixedValueFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    mappedPatchFieldBase<Type>::write(os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
