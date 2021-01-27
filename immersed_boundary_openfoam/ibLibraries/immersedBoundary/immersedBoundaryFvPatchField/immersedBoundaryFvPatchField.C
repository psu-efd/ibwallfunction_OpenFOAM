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

\*---------------------------------------------------------------------------*/

#include "immersedBoundaryFvPatchField.H"
#include "fvPatchFieldMapper.H"
#include "fvMatrix.H"
#include "ibSurfaceWriter.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
/*
template<class Type>
const int
immersedBoundaryFvPatchField<Type>::nBcIter_
(
    debug::optimisationSwitch
    (
        "immersedBoundaryNBCIter",
        5
    )
);
 

template<class Type>
const Foam::scalar
immersedBoundaryFvPatchField<Type>::bcTolerance_
(
    debug::optimisationSwitch
    (
        "immersedBoundaryBCTolerance",
        0.01
    )
);

*/
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void immersedBoundaryFvPatchField<Type>::updateIbValues() const
{

    // Create IB cell values first
    Field<Type> ibcv
    (
        this->internalField(),
        ibPatch_.ibCells()
    );

    if (this->fixesValue())
    {
        ibValue_ = ibPatch_.toIbPoints(refValue_);

        ibGrad_ = (ibValue_ - ibcv)/ibPatch_.ibDelta();

        // Reverse normals for external flow
        if (!ibPatch_.internalFlow())
        {
            ibGrad_ *= -1;
        }
    }
    else
    {
        ibGrad_ = ibPatch_.toIbPoints(refGrad_);
        //ibGrad_ = refGrad_;
        ibValue_ = ibcv + ibGrad_*ibPatch_.ibDelta();
    }
 
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
immersedBoundaryFvPatchField<Type>::imposeDirichletCondition() const
{
    //Pout<<"this->imposeDirichletCondition();"<<endl;
	word name_ = this->dimensionedInternalField().name();
	bool adaptive_ =ibmDict_.lookupOrDefault("AdaptiveToStencil",true);
    // Get addressing
    const labelList& ibc = ibPatch_.ibCells();
    const labelListList& ibCellCells = ibPatch_.ibCellCells();
    //const labelListList& procCells = ibPatch_.ibProcCells();
    const List<List<labelPair> >& ibCellProcCells = ibPatch_.ibCellProcCells();
    const PtrList<scalarRectangularMatrix>& invMat =
        ibPatch_.invDirichletMatrices();

    const vectorField& ibp = ibPatch_.ibPoints();

    const scalarField& ibD = ibPatch_.ibDelta();
    const scalarField& ibSD = ibPatch_.ibSamplingPointDelta();

    // Note: the algorithm is originally written with inward-facing normals
    // and subsequently changed: IB surface normals point outwards
    // HJ, 21/May/2012
    const vectorField& ibn = ibPatch_.ibNormals();

    // Collect Dirichlet values from IB triagulation
    //if(refValue_.empty())
    //{
    //    refValue_ = ibPatch_.toTriFaces(ibValue_);
    //}
    //else
	if(name_!="C")
	{
        ibValue_ = ibPatch_.toIbPoints(refValue_);
    }

    // Reset the size and value of snGrad
    ibGrad_.setSize(ibc.size());
    ibGrad_ = pTraits<Type>::zero;
 
    // Get access to internal field
    const Field<Type>& psiI = this->internalField();

    // Collect polynomially interpolated values in IB cells
    tmp<Field<Type> > tpolyPsi(new Field<Type>(psiI, ibc));
    Field<Type>& polyPsi = tpolyPsi();
 
    const vectorField& C = mesh_.cellCentres();

    // Dimension the matrix
    label nCoeffs = 5;
    if (mesh_.nGeometricD() == 3)
    {
        nCoeffs += 4;
    }
 
    label counter = 0;
    scalarField error(ibc.size(), 0);

    // Note
    // In parallel, some processors will have IB cells and others will not.
    // Therefore, special reduce is needed instead of gMax.
    // HJ, 7/Dec/2012
    scalar maxError = 0;

    do
    {
        counter++;

        // Parallel communication for psi
        FieldField<Field, Type> procPsi = ibPatch_.sendAndReceive(psiI);

        // Prepare error normalisation
        scalar maxMagPolyPsi = 0;

        //Field<Type> ibSPV = ibSamplingPointValue();

        forAll (ibc, cellI)
        {
            
			label curCell = ibc[cellI];

            const labelList& curCells = ibCellCells[cellI];

            const List<labelPair>& curProcCells = ibCellProcCells[cellI];

            Field<Type> coeffs(nCoeffs, pTraits<Type>::zero);

            Field<Type> source
            (
                curCells.size() + curProcCells.size(),
                pTraits<Type>::zero
            );
 
            label pointID = 0;
            for (label i = 0; i < curCells.size(); i++)
            {
                source[pointID++] = psiI[curCells[i]] - ibValue_[cellI];

            }

            for (label i = 0; i < curProcCells.size(); i++)
            {
                source[pointID++] =
                    procPsi
                    [
                        curProcCells[i].first()
                    ]
                    [
                        curProcCells[i].second()
                    ]
                  - ibValue_[cellI];
            }
 
			// if the stencil has less points than nCoeffs, change to lower order interpolation
			if(adaptive_ && source.size()<nCoeffs+3+10)
			{
				Field<Type> ibSPV = ibSamplingPointValue();
				Type oldPolyPsi = polyPsi[cellI];
				ibGrad_[cellI] = (ibSPV[cellI]-ibValue_[cellI])/ibSD[cellI];
				polyPsi[cellI] = ibGrad_[cellI]*ibD[cellI]+ibValue_[cellI];
				error[cellI] = mag(polyPsi[cellI] - oldPolyPsi);
				//Info<<oldPolyPsi<<" "<<polyPsi[cellI]<<" "<<ibSPV[cellI]<<" "<<ibValue_[cellI]<<endl;
				continue;	
			}

            const scalarRectangularMatrix& curInvMatrix = invMat[cellI];

            for (label i = 0; i < nCoeffs; i++)
            {
                for (label j = 0; j < source.size(); j++)
                {
                    coeffs[i] += curInvMatrix[i][j]*source[j];
                }
            }
 
            Type oldPolyPsi = polyPsi[cellI];
 
            vector R =  C[curCell] - ibp[cellI];
 
            polyPsi[cellI] =
                ibValue_[cellI]
              + coeffs[0]*R.x()
              + coeffs[1]*R.y()
              + coeffs[2]*R.x()*R.y()
              + coeffs[3]*sqr(R.x())
              + coeffs[4]*sqr(R.y());

            if (mesh_.nGeometricD() == 3)
            {
                polyPsi[cellI] +=
                 //   coeffs[2]*R.z();
                    coeffs[5]*R.z()
                  + coeffs[6]*R.x()*R.z()
                  + coeffs[7]*R.y()*R.z()
                  + coeffs[8]*sqr(R.z());
            }

            // Change of sign of ibn
            ibGrad_[cellI] =
               -coeffs[0]*ibn[cellI].x()
              - coeffs[1]*ibn[cellI].y();

            if (mesh_.nGeometricD() == 3)
            {
                // Change of sign of ibn
                ibGrad_[cellI] +=
                    -coeffs[5]*ibn[cellI].z();
            }

            error[cellI] = mag(polyPsi[cellI] - oldPolyPsi);
        }


        // Insert polynomial values into the internal field
        setIbCellValues(polyPsi);

        // Parallelisation fixes.  HJ, 7/Dec/2012

        // Reduce max polyPsi
        if (!polyPsi.empty())
        {
            maxMagPolyPsi = max(mag(polyPsi));
        }
        else
        {
            // No IB cells
            maxMagPolyPsi = 0;
        }

        reduce(maxMagPolyPsi, maxOp<scalar>());

        error /= maxMagPolyPsi + SMALL;

        // Reduce max error
        if (!polyPsi.empty())
        {
            maxError = max(error);
        }
        else
        {
            // No IB cells
            maxError = 0;
        }
 //Info<<"imposeDirichletCondition error "<<gAverage(error)<<" "<<counter<<endl;
        reduce(maxError, maxOp<scalar>());
 
    }
    while (maxError > bcTolerance_ && counter < nBcIter_);

    if (counter == nBcIter_ && debug)
    {
        InfoIn
        (
            "template<class Type>\n"
            "tmp<Foam::Field<Type> >\n"
            "immersedBoundaryFvPatchField<Type>::"
            "imposeDirichletCondition() const"
        )   << this->dimensionedInternalField().name()
            << " for patch " << this->patch().name()
            << ", error, max: " << gMax(error)
            << ", min: " << gMin(error)
            << ", avg: "  << gAverage(error) << endl;
    }
//Info<<"gAverage(polyPsi) "<<gAverage(polyPsi)<<endl;
    return tpolyPsi;
}

template<class Type>
Foam::tmp<Foam::Field<Type> >
immersedBoundaryFvPatchField<Type>::imposeNeumannCondition() const
{
	word name_ = this->dimensionedInternalField().name();

	bool adaptive_ =ibmDict_.lookupOrDefault("AdaptiveToStencil",true);
    //Pout<<"this->imposeNeumannCondition(); "<<   this->dimensionedInternalField().name()<<endl;
    // Get addressing
    const labelList& ibc = ibPatch_.ibCells();

    const labelListList& ibCellCells = ibPatch_.ibCellCells();

    //const labelListList& procCells = ibPatch_.ibProcCells();

    const List<List<labelPair> >& ibCellProcCells = ibPatch_.ibCellProcCells();

    const PtrList<scalarRectangularMatrix>& invMat =
        ibPatch_.invNeumannMatrices();

    const vectorField& ibp = ibPatch_.ibPoints();
    const scalarField& ibD = ibPatch_.ibDelta();
    const scalarField& ibSD = ibPatch_.ibSamplingPointDelta();

    // Collect Neumann values from IB triagulation
	if(name_!="C")
    {
		ibGrad_ = ibPatch_.toIbPoints(refGrad_);
	}

    //ibGrad_ = refGrad_;
    // Reset the size and value of
    ibValue_.setSize(ibc.size());
    ibValue_ = pTraits<Type>::zero;

    // Get access to internal field
    const Field<Type>& psiI = this->internalField();

    // Collect polynomially interpolated values in IB cells
    tmp<Field<Type> > tpolyPsi(new Field<Type>(psiI, ibc));
    Field<Type>& polyPsi = tpolyPsi();

    const vectorField& C = mesh_.cellCentres();

    // Dimension the matrix
    label nCoeffs = 6;
    if (mesh_.nGeometricD() == 3)
    {
        nCoeffs += 4;
    }

    label counter = 0;
    scalarField error(ibc.size(), 0);

    // Initialise ibCell values using sampling point values to reduce
    // the number of iterations.  HJ, 26/Oct/2012
    setIbCellValues(ibSamplingPointValue());


    // Note
    // In parallel, some processors will have IB cells and others will not.
    // Therefore, special reduce is needed instead of gMax.
    // HJ, 7/Dec/2012
    scalar maxError = 0;

    do
    {
        counter++;

        // Parallel communication for psi
        FieldField<Field, Type> procPsi = ibPatch_.sendAndReceive(psiI);

        // Prepare error normalisation
        scalar maxMagPolyPsi = 0;

        const Field<Type> ibSPV = ibSamplingPointValue();

        forAll (ibc, cellI)
        {

            label curCell = ibc[cellI];

            const labelList& curCells = ibCellCells[cellI];
            const List<labelPair>& curProcCells = ibCellProcCells[cellI];

            Field<Type> coeffs(nCoeffs, pTraits<Type>::zero);

            Field<Type> source
            (
                curCells.size() + 1 + curProcCells.size(),
                pTraits<Type>::zero
            );
			// if the stencil has less points than nCoeffs, change to lower order interpolation
			if(adaptive_ && source.size()<nCoeffs+3+10)
			{
				//Field<Type> ibSPV = ibSamplingPointValue();
				Type oldPolyPsi = polyPsi[cellI];
				polyPsi[cellI] = ibSPV[cellI]-ibGrad_[cellI]*(ibSD[cellI]-ibD[cellI]);
				ibValue_[cellI] = ibSPV[cellI]-ibGrad_[cellI]*ibSD[cellI];
				error[cellI] = mag(polyPsi[cellI] - oldPolyPsi);
				//Info<<oldPolyPsi<<" "<<polyPsi[cellI]<<" "<<ibSPV[cellI]<<" "<<ibValue_[cellI]<<endl;
				continue;
			}
 
            const scalarRectangularMatrix& curInvMatrix = invMat[cellI];

            label pointID = 0;

            for (label i = 0; i < curCells.size(); i++)
            {
                source[pointID++] = psiI[curCells[i]];
            }

            source[pointID++] = ibGrad_[cellI];

            for (label i = 0; i < curProcCells.size(); i++)
            {
                source[pointID++] =
                    procPsi
                    [
                        curProcCells[i].first()
                    ]
                    [
                        curProcCells[i].second()
                    ];
            }

            for (label i = 0; i < nCoeffs; i++)
            {
                for (label j = 0; j < source.size(); j++)
                {
                    coeffs[i] += curInvMatrix[i][j]*source[j];
                }
            }

            Type oldPolyPsi = polyPsi[cellI];

            vector ibR =  C[curCell] - ibp[cellI];

            polyPsi[cellI] =
                coeffs[0]
              + coeffs[1]*ibR.x()
              + coeffs[2]*ibR.y()
              + coeffs[3]*ibR.x()*ibR.y()
              + coeffs[4]*sqr(ibR.x())
              + coeffs[5]*sqr(ibR.y());

            if (mesh_.nGeometricD() == 3)
            {
                polyPsi[cellI] +=
                    coeffs[6]*ibR.z()
                  + coeffs[7]*ibR.x()*ibR.z()
                  + coeffs[8]*ibR.y()*ibR.z()
                  + coeffs[9]*sqr(ibR.z());
            }

            ibValue_[cellI] = coeffs[0];

            error[cellI] = mag(polyPsi[cellI] - oldPolyPsi);
        }

        // Insert polynomial values into the internal field
        setIbCellValues(polyPsi);

        // Parallelisation fixes.  HJ, 7/Dec/2012

        // Reduce max polyPsi
        if (!polyPsi.empty())
        {
            maxMagPolyPsi = max(mag(polyPsi));
        }
        else
        {
            // No IB cells
            maxMagPolyPsi = 0;
        }

        reduce(maxMagPolyPsi, maxOp<scalar>());

        error /= maxMagPolyPsi + SMALL;

        // Reduce max error
        if (!polyPsi.empty())
        {
            maxError = max(error);
        }
        else
        {
            // No IB cells
            maxError = 0;
        }

        reduce(maxError, maxOp<scalar>());
 //Info<<"imposeNeumannCondition error "<<gAverage(error)<<" "<<counter<<endl;
    }
    while (maxError > bcTolerance_ && counter < nBcIter_);

    //Info <<maxError<<" "<<counter<<endl; 

    if (counter == nBcIter_ && debug)
    {
        InfoIn
        (
            "template<class Type>\n"
            "tmp<Foam::Field<Type> >\n"
            "immersedBoundaryFvPatchField<Type>::"
            "imposeNeumannCondition() const"
        )   << this->dimensionedInternalField().name()
            << " for patch " << this->patch().name()
            << ", error, max: " << gMax(error)
            << ", min: " << gMin(error)
            << ", avg: "  << gAverage(error) << endl;
    }

//     Info<< "Neumann condition for " << ibc.size() << " cells of field "
//         << this->dimensionedInternalField().name() << " = "
//         << polyPsi
//         << endl;
//Info<<"gAverage(polyPsi) "<<gAverage(polyPsi)<<endl;
    return tpolyPsi;
}

template<class Type>
Foam::tmp<Foam::Field<Type> >
immersedBoundaryFvPatchField<Type>::imposeNeumannCondition_low_order() const
{
	word name_ = this->dimensionedInternalField().name();
 
    // Get addressing
    const labelList& ibc = ibPatch_.ibCells();

    //const labelListList& procCells = ibPatch_.ibProcCells();

    const scalarField& ibD = ibPatch_.ibDelta();
    const scalarField& ibSD = ibPatch_.ibSamplingPointDelta();

    // Collect Neumann values from IB triagulation
	if(name_!="C")
    {
		ibGrad_ = ibPatch_.toIbPoints(refGrad_);
	}

    //ibGrad_ = refGrad_;
    // Reset the size and value of
    ibValue_.setSize(ibc.size());
    ibValue_ = pTraits<Type>::zero;

    // Get access to internal field
    const Field<Type>& psiI = this->internalField();

    // Collect polynomially interpolated values in IB cells
    tmp<Field<Type> > tpolyPsi(new Field<Type>(psiI, ibc));
    Field<Type>& polyPsi = tpolyPsi();

    // Dimension the matrix
    label nCoeffs = 6;
    if (mesh_.nGeometricD() == 3)
    {
        nCoeffs += 4;
    }

    label counter = 0;
    scalarField error(ibc.size(), 0);

    // Initialise ibCell values using sampling point values to reduce
    // the number of iterations.  HJ, 26/Oct/2012
    setIbCellValues(ibSamplingPointValue());


    // Note
    // In parallel, some processors will have IB cells and others will not.
    // Therefore, special reduce is needed instead of gMax.
    // HJ, 7/Dec/2012
    scalar maxError = 0;

    do
    {
        counter++;

        // Parallel communication for psi
        //FieldField<Field, Type> procPsi = ibPatch_.sendAndReceive(psiI,procCells);

        // Prepare error normalisation
        scalar maxMagPolyPsi = 0;

        const Field<Type> ibSPV = ibSamplingPointValue();

        forAll (ibc, cellI)
        {
			//Field<Type> ibSPV = ibSamplingPointValue();
			Type oldPolyPsi = polyPsi[cellI];
			polyPsi[cellI] = ibSPV[cellI]-ibGrad_[cellI]*(ibSD[cellI]-ibD[cellI]);
			ibValue_[cellI] = ibSPV[cellI]-ibGrad_[cellI]*ibSD[cellI];
			error[cellI] = mag(polyPsi[cellI] - oldPolyPsi);
			//Info<<oldPolyPsi<<" "<<polyPsi[cellI]<<" "<<ibSPV[cellI]<<" "<<ibValue_[cellI]<<endl;
        }

        // Insert polynomial values into the internal field
        setIbCellValues(polyPsi);

        // Parallelisation fixes.  HJ, 7/Dec/2012

        // Reduce max polyPsi
        if (!polyPsi.empty())
        {
            maxMagPolyPsi = max(mag(polyPsi));
        }
        else
        {
            // No IB cells
            maxMagPolyPsi = 0;
        }

        reduce(maxMagPolyPsi, maxOp<scalar>());

        error /= maxMagPolyPsi + SMALL;

        // Reduce max error
        if (!polyPsi.empty())
        {
            maxError = max(error);
        }
        else
        {
            // No IB cells
            maxError = 0;
        }

        reduce(maxError, maxOp<scalar>());
 
    }
    while (maxError > bcTolerance_ && counter < nBcIter_);

    if (counter == nBcIter_ && debug)
    {
        InfoIn
        (
            "template<class Type>\n"
            "tmp<Foam::Field<Type> >\n"
            "immersedBoundaryFvPatchField<Type>::"
            "imposeNeumannCondition() const"
        )   << this->dimensionedInternalField().name()
            << " for patch " << this->patch().name()
            << ", error, max: " << gMax(error)
            << ", min: " << gMin(error)
            << ", avg: "  << gAverage(error) << endl;
    }
    return tpolyPsi;
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
immersedBoundaryFvPatchField<Type>::imposeNeumannCondition2nd() const
{

	word name_ = this->dimensionedInternalField().name();
    const labelList& ibc = ibPatch_.ibCells();

    // Collect Neumann values from IB triagulation
	if(name_!="C")
    {
		ibGrad_ = ibPatch_.toIbPoints(refGrad_);
	}
    // Reset the size and value of
    ibValue_.setSize(ibc.size());
    ibValue_ = pTraits<Type>::zero;

    // Get access to internal field
    const Field<Type>& psiI = this->internalField();

    // Collect polynomially interpolated values in IB cells
    tmp<Field<Type> > tpolyPsi(new Field<Type>(psiI, ibc));

    Field<Type>& polyPsi = tpolyPsi();

    const Field<Type> ibSPV = ibSamplingPointValue();

    const scalarField& ibD = ibPatch_.ibDelta();
    const scalarField& ibSD = ibPatch_.ibSamplingPointDelta();
/*    const vectorField& ibN = ibPatch_.ibNormals();

    word name_ = this->dimensionedInternalField().name();

    typedef GeometricField<Type, fvPatchField, volMesh> volTypeField;

    volTypeField& this_old = const_cast<volTypeField&>
        (mesh_.lookupObject<volTypeField>(name_));
*/
    //if(name_ == "grad(p)" or name_ == "grad(k)" or name_ == "grad(epsilon)")
    //{
	forAll(ibSPV,ibCellI)
	{
    	ibValue_[ibCellI] = ibSPV[ibCellI]-ibGrad_[ibCellI]*ibSD[ibCellI];
    	polyPsi[ibCellI] = ibSPV[ibCellI]-ibGrad_[ibCellI]*(ibSD[ibCellI]-ibD[ibCellI]);
	}
    //}
/*
    else if(name_ == "p" or name_ == "k" or name_ == "epsilon" )
    {   
        const volScalarField& psiI1 = refCast<const volScalarField>(this_old);

        vectorField tempGrad = fvc::grad(psiI1);
 
        vectorField vectorGradField_old = ibPatch_.toSamplingPoints(tempGrad);

        forAll(vectorGradField_old, ibCellI)
        {
            const vector& vectorGrad_old = vectorGradField_old[ibCellI];
            const vector& N = ibN[ibCellI];

            scalar ibGrad = N & vectorGrad_old;

            polyPsi[ibCellI] = ibSPV[ibCellI] - (ibSD[ibCellI]-ibD[ibCellI])*(ibGrad*pTraits<Type>::one);
            ibValue_[ibCellI] = ibSPV[ibCellI] - ibSD[ibCellI]*(ibGrad*pTraits<Type>::one);
        }
    }
    else
    {
	    forAll(ibSPV,ibCellI)
	    {
	        ibValue_[ibCellI] = ibSPV[ibCellI];
	        polyPsi[ibCellI] = ibSPV[ibCellI];
	    }
    }
*/
    setIbCellValues(polyPsi);

    return tpolyPsi;
}



template<class Type>
void immersedBoundaryFvPatchField<Type>::imposeDeadCondition()
{
//    Pout<<"this->imposeDeadCondition();"<<endl;
    const labelList& dc = ibPatch_.deadCells();

//     Info<< "Dead condition for " << dc.size()  << " cells of field "
//         << this->dimensionedInternalField().name()
//         << " set to value " << deadCellValue_
//         << endl;

    // Get non-const access to internal field
    Field<Type>& psiI = const_cast<Field<Type>&>(this->internalField());

    forAll (dc, dcI)
    {
        psiI[dc[dcI]] = deadCellValue_;
    }
}


template<class Type>
void immersedBoundaryFvPatchField<Type>::correctDiag
(
    fvMatrix<Type>& eqn
) const
{
    scalarField& Diag = eqn.diag();

    const labelList& dce = ibPatch_.deadCellsExt();

    // Estimate diagonal in live cells
    scalar liveDiag = 1;

    // In case dce.size() == Diag.size(),  added by Y.C. Xu 2017/6
    scalar gSumMagDiag = gSumMag(Diag);
    scalar gMaxDiag = gMax(Diag);

    if (dce.size() < Diag.size())
    {
    // Problem: if dce.size() == Diag.size() in one processor, gSumMag(Diag) here may not be solved.
    // added by Y.C. Xu 2017/6
        //liveDiag = gSumMag(Diag)/(Diag.size() - dce.size()); 
        liveDiag = gSumMagDiag/(Diag.size() - dce.size()); 

        // Correct for sign
        //liveDiag *= sign(gMax(Diag));
        liveDiag *= sign(gMaxDiag);
    }


    forAll (dce, cellI)
    {
        if (mag(Diag[dce[cellI]]) < SMALL)
        {
            Diag[dce[cellI]] = liveDiag;
        }
    }

}


template<class Type>
void immersedBoundaryFvPatchField<Type>::correctOffDiag // impose fixed gradient condition, modified by Xu 7/2017
(
    fvMatrix<Type>& eqn
) const
{
    // Calculate gradient contribution
    const labelList& ibFaces = ibPatch_.ibFaces();
    const labelList& ibFaceCells = ibPatch_.ibFaceCells();
    //const labelList& ibc = ibPatch_.ibCells();


    //Info<<"correctOffDiag "<< ibFaces.size()<<" "<<ibCells.size()<<endl;

    const scalarField& ibGamma = ibPatch_.gamma().internalField();

    const unallocLabelList& own = mesh_.owner();
    const unallocLabelList& nei = mesh_.neighbour();

    // Get delta coefficients
    const surfaceScalarField& dc = mesh_.deltaCoeffs();
    const scalarField& dcI = dc.internalField();

    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();


    // assign new values of ibFaceGrad
    tmp<vectorField > tibFaceGrad
    (
        new vectorField(ibFaces.size(), pTraits<vector>::zero)
    );

    vectorField& ibFaceGrad = tibFaceGrad();

    word name_ = this->dimensionedInternalField().name();

    volScalarField& psiI = const_cast<volScalarField&>
        (mesh_.lookupObject<volScalarField>(name_));

    typedef GeometricField<Type,fvPatchField,volMesh> volTypeField;

    volTypeField& psiI2 = const_cast<volTypeField&>
        (mesh_.lookupObject<volTypeField>(name_));


    word tNname=pTraits<Type>::typeName;

    //if(name_ == "p" or name_ == "k" or name_ == "epsilon")
    //if(tNname == "scalar")
    {   
        Info<< "CorrectOffDiag for field "
           << this->dimensionedInternalField().name() << endl;         

        // calculate the gradient on the internal faces
		const vectorField gradFieldSI = linearInterpolate(fvc::grad(psiI));

        // calculate the gradient on the proc faces

				forAll(psiI.boundaryField(),patchi)
				{
					if(psiI.boundaryField()[patchi].coupled())
					{
			            psiI.boundaryField()[patchi].initEvaluate(Pstream::blocking);// not quite sure
					}
                }
				forAll(psiI.boundaryField(),patchi)
				{
					if(psiI.boundaryField()[patchi].coupled())
					{
			            psiI.boundaryField()[patchi].evaluate(Pstream::blocking);// not quite sure
					}
                }

        forAll (ibFaces, faceI)
        {
            const label curFace = ibFaces[faceI];

            if (curFace < nei.size())
            {   
                ibFaceGrad[faceI] = gradFieldSI[curFace];
            }
        }
    }
    

    if (eqn.symmetric())
    {
        scalarField& diag = eqn.diag();
        scalarField& upper = eqn.upper();
        Field<Type>& source = eqn.source();

        Info<< "Symmetric correctOffDiag for field "
             << this->dimensionedInternalField().name() << endl;

        forAll (ibFaces, faceI)
        {

            const label curFace = ibFaces[faceI];

            if (curFace < nei.size())
            {

                // calculate surface vector on ibFaces    added by Xu 7/2017
                const vector SfN = Sf[curFace]/magSf[curFace];

                // calculate value gradient on ibFaces o->n     added by Xu 7/2017
                const Type faceGrad = (SfN&ibFaceGrad[faceI]) * pTraits<Type>::one;

                // Internal face.  One side is an ibCell and another is a
                // live cell. Add gradient to the source of the live cell
                // and kill the off-diagonal coefficient
                if (ibGamma[own[curFace]] > SMALL)
                {
                    diag[own[curFace]] += upper[curFace];

                    //source[own[curFace]] +=
                        //upper[curFace]*ibGrad_[ibFaceCells[faceI]]
                        ///dcI[curFace];

                    source[own[curFace]] -=
                        upper[curFace]*faceGrad/dcI[curFace];

                }
                else
                {
                    diag[nei[curFace]] += upper[curFace];

                    //source[nei[curFace]] -=
                        //upper[curFace]*ibGrad_[ibFaceCells[faceI]]
                        ///dcI[curFace];

                    source[nei[curFace]] +=
                        upper[curFace]*faceGrad/dcI[curFace];

                }

                upper[curFace] = 0.0;
            }
            else
            {
                // calculate value gradient on ibFaces o->n     added by Xu 7/2017
                //const Type faceGrad = ibFaceGrad[faceI] * pTraits<Type>::one;

                // else MPI
                label patchi = mesh_.boundaryMesh().whichPatch(curFace);

                if (!eqn.internalCoeffs()[patchi].empty())
                {
                    label patchFacei =
                        mesh_.boundaryMesh()[patchi].whichFace(curFace);

                    const fvPatchField<Type>& ptf = psiI2.boundaryField()[patchi];
       		        const tmp<Field<Type> > tpnf = ptf.patchNeighbourField();
                    const Type& pnff = tpnf()[patchFacei];

       		        tmp<Field<scalar> > tpic = eqn.internalCoeffs()[patchi].component(0);
                    const scalar picf = tpic()[patchFacei];

       		        tmp<Field<Type> > tpbc = eqn.boundaryCoeffs()[patchi];
                    const Type pbcf = tpbc()[patchFacei];
                    //const Type pbcf=eqn.boundaryCoeffs()[patchi][patchFacei];

                    eqn.internalCoeffs()[patchi][patchFacei] =
                        pTraits<Type>::zero;

                    eqn.boundaryCoeffs()[patchi][patchFacei] =
                        pTraits<Type>::zero;

                    // Check if the live cell is on local or neighbour side
                    // HJ, 7/Dec/2012
                    if (ibFaceCells[faceI] > -1)
                    {
                        if (ibGamma[ibFaceCells[faceI]] > SMALL)
                        {
                            //source[ibFaceCells[faceI]] +=
                                //ibGrad_[ibFaceCells[faceI]]
                                ///dc.boundaryField()[patchi][patchFacei];

                            //source[ibFaceCells[faceI]] +=
                               // faceGrad
                               // /dc.boundaryField()[patchi][patchFacei];

			                diag[ibFaceCells[faceI]] += picf;

							//scalar dt=dc.boundaryField()[patchi][patchFacei];
                            //Type piff=psiI2[ibc[ibFaceCells[faceI]]];

                            source[ibFaceCells[faceI]] +=
                                 cmptMultiply(pbcf, pnff);
                        }
                    }
                }
            }
        }
    }
    else if (eqn.asymmetric())
    {
        scalarField& diag = eqn.diag();
        scalarField& upper = eqn.upper();
        scalarField& lower = eqn.lower();
        Field<Type>& source = eqn.source();

        Info<< "Asymmetric correctOffDiag for field "
           << this->dimensionedInternalField().name() << endl;

        forAll (ibFaces, faceI)
        {
            const label curFace = ibFaces[faceI];

            if (curFace < nei.size())
            {

                // calculate surface vector on ibFaces    added by Xu 7/2017
                const vector SfN = Sf[curFace]/magSf[curFace];

                // calculate value gradient on ibFaces o->n     added by Xu 7/2017
                const Type faceGrad = (SfN&ibFaceGrad[faceI]) * pTraits<Type>::one;

                // Internal face.  One side is an ibCell and another is a
                // live cell. Add gradient to the source of the live cell
                // and kill the off-diagonal coefficient
                if (ibGamma[own[curFace]] > SMALL)
                {
                    diag[own[curFace]] += upper[curFace];

                    source[own[curFace]] -=
                        upper[curFace]*faceGrad/dcI[curFace];
                    //source[own[curFace]] +=
                        //upper[curFace]*ibGrad_[ibFaceCells[faceI]]/dcI[faceI];
                }
                else
                {
                    diag[nei[curFace]] += lower[curFace];

                    source[nei[curFace]] +=
                        lower[curFace]*faceGrad/dcI[curFace];
                    //source[nei[curFace]] -=
                        //lower[curFace]*ibGrad_[ibFaceCells[faceI]]/dcI[faceI];

                }

                upper[curFace] = 0.0;
                lower[curFace] = 0.0;
            }
            else
            {

                // calculate value gradient on ibFaces o->n     added by Xu 7/2017
                //const Type faceGrad = ibFaceGrad[faceI] * pTraits<Type>::one;

                // else MPH
                label patchi = mesh_.boundaryMesh().whichPatch(curFace);

                if (!eqn.internalCoeffs()[patchi].empty())
                {
                    label patchFacei =
                        mesh_.boundaryMesh()[patchi].whichFace(curFace);

                    const fvPatchField<Type>& ptf = psiI2.boundaryField()[patchi];
       		        const tmp<Field<Type> > tpnf = ptf.patchNeighbourField();
                    const Type& pnff = tpnf()[patchFacei];

       		        tmp<Field<scalar> > tpic = eqn.internalCoeffs()[patchi].component(0);
                    const scalar picf = tpic()[patchFacei];

       		        tmp<Field<Type> > tpbc = eqn.boundaryCoeffs()[patchi];
                    const Type pbcf = tpbc()[patchFacei];

                    eqn.internalCoeffs()[patchi][patchFacei] =
                        pTraits<Type>::zero;

                    eqn.boundaryCoeffs()[patchi][patchFacei] =
                        pTraits<Type>::zero;

                    // Check if the live cell is on local or neighbour side
                    // HJ, 7/Dec/2012
                    if (ibFaceCells[faceI] > -1)
                    {
                        if (ibGamma[ibFaceCells[faceI]] > SMALL)
                        { 
                            //source[ibFaceCells[faceI]] +=
                                //ibGrad_[ibFaceCells[faceI]]
                                ///dc.boundaryField()[patchi][patchFacei];

                            //source[ibFaceCells[faceI]] +=
                               // faceGrad
                               // /dc.boundaryField()[patchi][patchFacei];

			                diag[ibFaceCells[faceI]] += picf;

							//scalar dt=dc.boundaryField()[patchi][patchFacei];
                            //Type piff=psiI2[ibc[ibFaceCells[faceI]]];

                            source[ibFaceCells[faceI]] +=
                                 cmptMultiply(pbcf, pnff);
                        }
                    }
                }
            }
        }
 
    }

    // Note: potentially deal with face flux correction ptr.
    // HJ, 16/Apr/2012
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class Type>
bool immersedBoundaryFvPatchField<Type>::motionUpdateRequired() const
{
    if
    (
        ibPatch_.movingIb()
     || ibPatch_.boundaryMesh().mesh().moving()
    )
    {
        if
        (
            ibUpdateTimeIndex_
         != ibPatch_.boundaryMesh().mesh().time().timeIndex()
        )
        {
            // Mesh is moving and current time has not been updated
            return true;
        }
    }

    return false;
}


template<class Type>
void immersedBoundaryFvPatchField<Type>::motionUpdate() const
{
    // Motion update, clear data related to immersed boundary points
    ibValue_.clear();
    ibGrad_.clear();

    // Record motion update time
    ibUpdateTimeIndex_ = ibPatch_.boundaryMesh().mesh().time().timeIndex();
}


template<class Type>
void immersedBoundaryFvPatchField<Type>::setIbCellValues
(
    const Field<Type>& ibcValues
) const
{
    const labelList& ibc = ibPatch_.ibCells();

    if (ibcValues.size() != ibc.size())
    {
        FatalErrorIn
        (
            "template<class Type>\n"
            "void immersedBoundaryFvPatchField<Type>::setIbCellValues\n"
            "(\n"
            "    const Field<Type>& ibcValues\n"
            ") const"
        )   << "Size of ibcValues field not equal to the number of IB cells."
            << nl << "ibcValues: " << ibcValues.size()
            << " ibc: " << ibc.size()
            << abort(FatalError);
    }

    // Get non-const access to internal field
    Field<Type>& psiI = const_cast<Field<Type>&>(this->internalField());

    forAll (ibcValues, cellI)
    {
        psiI[ibc[cellI]] = ibcValues[cellI];
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
immersedBoundaryFvPatchField<Type>::immersedBoundaryFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(p, iF),
    ibPatch_(refCast<const immersedBoundaryFvPatch>(p)),
    mesh_(p.boundaryMesh().mesh()),
    refValue_(ibPatch_.ibMesh().size(), pTraits<Type>::zero),
    refGrad_(ibPatch_.ibMesh().size(), pTraits<Type>::zero),
    fixesValue_(false),
    setDeadCellValue_(false),
    deadCellValue_(pTraits<Type>::zero),
    ibUpdateTimeIndex_(-1),
    ibValue_(),
    ibGrad_(),
	ibmDict_        
		(
            IOobject
            (
                "ibmDict",
                mesh_.time().system(),
                mesh_,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        )
{}


template<class Type>
immersedBoundaryFvPatchField<Type>::immersedBoundaryFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<Type>(p, iF, dict),
    ibPatch_(refCast<const immersedBoundaryFvPatch>(p)),
    mesh_(p.boundaryMesh().mesh()),
    refValue_("refValue", dict, ibPatch_.ibMesh().size()),
    refGrad_("refGradient", dict, ibPatch_.ibMesh().size()),
    fixesValue_(dict.lookup("fixesValue")),
    setDeadCellValue_(dict.lookup("setDeadCellValue")),
    deadCellValue_(pTraits<Type>(dict.lookup("deadCellValue"))),
    ibUpdateTimeIndex_(-1),
    ibValue_(),
    ibGrad_(),
	ibmDict_        
		(
            IOobject
            (
                "ibmDict",
                mesh_.time().system(),
                mesh_,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        )
{
    if (!isType<immersedBoundaryFvPatch>(p))
    {
        FatalIOErrorIn
        (
            "immersedBoundaryFvPatchField<Type>::"
            "immersedBoundaryFvPatchField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const Field<Type>& field,\n"
            "    const dictionary& dict\n"
            ")\n",
            dict
        )   << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }
}


template<class Type>
immersedBoundaryFvPatchField<Type>::immersedBoundaryFvPatchField
(
    const immersedBoundaryFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchField<Type>(ptf, p, iF, mapper),
    ibPatch_(refCast<const immersedBoundaryFvPatch>(p)),
    mesh_(p.boundaryMesh().mesh()),
    refValue_(ptf.refValue()),
    refGrad_(ptf.refGrad()),
    fixesValue_(ptf.fixesValue()),
    setDeadCellValue_(ptf.setDeadCellValue_),
    deadCellValue_(ptf.deadCellValue_),
    ibUpdateTimeIndex_(-1),
    ibValue_(),
    ibGrad_(),
	ibmDict_(ptf.ibmDict_)
{
    // Note: NO MAPPING.  Fields are created on the immersed boundary
    // HJ, 12/Apr/2012
    if (!isType<immersedBoundaryFvPatch>(p))
    {
        FatalErrorIn
        (
            "immersedBoundaryFvPatchField<Type>::"
            "immersedBoundaryFvPatchField\n"
            "(\n"
            "    const immersedBoundaryFvPatchField<Type>&,\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<Type, volMesh>& iF,\n"
            "    const fvPatchFieldMapper& mapper\n"
            ")\n"
        )   << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }
}


template<class Type>
immersedBoundaryFvPatchField<Type>::immersedBoundaryFvPatchField
(
    const immersedBoundaryFvPatchField<Type>& ptf
)
:
    fvPatchField<Type>(ptf),
    ibPatch_(ptf.ibPatch()),
    mesh_(ptf.patch().boundaryMesh().mesh()),
    refValue_(ptf.refValue()),
    refGrad_(ptf.refGrad()),
    fixesValue_(ptf.fixesValue()),
    setDeadCellValue_(ptf.setDeadCellValue_),
    deadCellValue_(ptf.deadCellValue_),
    ibUpdateTimeIndex_(-1),
    ibValue_(),
    ibGrad_(),
	ibmDict_(ptf.ibmDict_)
{}


template<class Type>
immersedBoundaryFvPatchField<Type>::immersedBoundaryFvPatchField
(
    const immersedBoundaryFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(ptf, iF),
    ibPatch_(ptf.ibPatch()),
    mesh_(ptf.patch().boundaryMesh().mesh()),
    refValue_(ptf.refValue()),
    refGrad_(ptf.refGrad()),
    fixesValue_(ptf.fixesValue()),
    setDeadCellValue_(ptf.setDeadCellValue_),
    deadCellValue_(ptf.deadCellValue_),
    ibUpdateTimeIndex_(-1),
    ibValue_(),
    ibGrad_(),
	ibmDict_(ptf.ibmDict_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Field<Type>& immersedBoundaryFvPatchField<Type>::ibValue() const
{
//Info<<"Field<Type>& immersedBoundaryFvPatchField<Type>::ibValue() const"<<endl;
    // Note: on a moving mesh, the intersection has changed and
    // ibValue and ibGrad fields should be cleared and recalculated.
    // HJ, 17/Oct/2012
    if (this->motionUpdateRequired())
    {
        motionUpdate();
    }

    if (ibValue_.empty())
    {
 
        this->updateIbValues();
    }
 
    return ibValue_;
}


template<class Type>
const Field<Type>& immersedBoundaryFvPatchField<Type>::ibGrad() const
{
    // Note: on a moving mesh, the intersection has changed and
    // ibValue and ibGrad fields should be cleared and recalculated.
    // HJ, 17/Oct/2012

    if (this->motionUpdateRequired())
    {
        motionUpdate();
    }

    //if (ibGrad_.empty())
    {
        this->updateIbValues();
    }

    return ibGrad_;
}


template<class Type>
tmp<Field<Type> > immersedBoundaryFvPatchField<Type>::ibCellValue() const
{
    // Collect IB cell values
    tmp<Field<Type> > tibcv
    (
        new Field<Type>
        (
            this->internalField(),
            ibPatch_.ibCells()
        )
    );

    return tibcv;
}


template<class Type>
tmp<Field<Type> >
immersedBoundaryFvPatchField<Type>::ibSamplingPointValue() const
{
    return ibPatch_.toSamplingPoints(this->internalField());
}


template<class Type>
tmp<Field<Type> > immersedBoundaryFvPatchField<Type>::triValue() const
{
    return ibPatch_.toTriFaces(this->ibValue());
}


template<class Type>
tmp<Field<Type> > immersedBoundaryFvPatchField<Type>::triGrad() const
{
    return ibPatch_.toTriFaces(this->ibGrad());
}



template<class Type>
void immersedBoundaryFvPatchField<Type>::updateCoeffs()
{
    nBcIter_=5;
    bcTolerance_=0.001;

	bool high_order_ =ibmDict_.lookupOrDefault("NeumannConditionHighOrder",true);
    word name_ = this->dimensionedInternalField().name();
    //Info<< " immersedBoundaryFvPatchField::updateCoeffs() "<<  name_<<" "<<this->fixesValue()<<endl;
    if (this->fixesValue() )
    {
    	//Info<<"this->imposeDirichletCondition(); "<<   name_<<endl;
        //this->imposeDirichletCondition();
        //this->imposeDirichletCondition2nd();
		//if(high_order_)
 		{
			this->imposeDirichletCondition();
		}
		//else
		{
			//this->imposeDirichletCondition2nd();
		}
    }
    else  
    {

        //Info<<"this->imposeNeumannCondition(); "<<   name_<<endl;

    	//this->imposeNeumannCondition();
	    //if(name_ == "grad(p)" or name_ == "grad(k)" or name_ == "grad(epsilon)")
		if(high_order_)
 		{
			this->imposeNeumannCondition();
		}
		else
		{
			if(name_!="C")
			{	
	     		this->imposeNeumannCondition_low_order();
			}
			else
			{
	     		this->imposeNeumannCondition();
			}
//const Field<Type>& dd = const_cast<Field<Type>&>
//        (mesh_.lookupObject<Field<Type> >(name_));
//Info<<name_<<" "<<min(dd)<<endl;
		}

    }

    // Fix the value in dead cells
    if (setDeadCellValue_)
    {
        this->imposeDeadCondition();
    }
 
    fvPatchField<Type>::updateCoeffs();
}
 

//Impose the effect of immersed boundary
template<class Type>
void immersedBoundaryFvPatchField<Type>::updateIbCellValues()
{ 
    fvPatchField<Type>::updateCoeffs();
}

template<class Type>
void immersedBoundaryFvPatchField<Type>::initEvaluate
(
    const Pstream::commsTypes
)
{
    this->updateCoeffs();
}


template<class Type>
void immersedBoundaryFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes
)
{

    // Note
    // Since the boundary condition is performed by data fitting with the
    // internal field, fitting must be performed both on updateCoeffs
    // and on evaluate (internal field has changed in the meantime).
    // Bug fix, Zeljko Tukovic, 21/Jun/2012
     this->updateCoeffs();

    fvPatchField<Type>::evaluate();
}


template<class Type>
void immersedBoundaryFvPatchField<Type>::manipulateMatrix
(
    fvMatrix<Type>& eqn
)
{
    this->initEvaluate();

    // Build matrix diagonal for cells where it is missing
    this->correctDiag(eqn);

    // For Neumann boundary condition, manipulate matrix off-diagonal
    if (!this->fixesValue())
    {
        this->correctOffDiag(eqn);
    }

    //if (this->fixesValue())
    {
    	// Set values in IB cells
    	Field<Type> polyPsi(eqn.psi(), ibPatch_.ibCells());
    	eqn.setValues(ibPatch_.ibCells(), polyPsi);
    }

    // Correct equation for dead cells
    Field<Type> deadCellsPsi
    (
        ibPatch_.deadCells().size(),
        deadCellValue_
    );
    eqn.setValues(ibPatch_.deadCells(), deadCellsPsi);

    fvPatchField<Type>::manipulateMatrix(eqn);
}


template<class Type>
void immersedBoundaryFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    refValue_.writeEntry("refValue", os);
    refGrad_.writeEntry("refGradient", os);
    os.writeKeyword("fixesValue") << fixesValue_ << token::END_STATEMENT << nl;
    os.writeKeyword("setDeadCellValue")
        << setDeadCellValue_ << token::END_STATEMENT << nl;
    os.writeKeyword("deadCellValue")
        << deadCellValue_ << token::END_STATEMENT << nl;

    this->writeEntry("value", os);
 /*
    // Write immersed boundary data as a vtk file
    autoPtr<ibSurfaceWriter<Type> > writerPtr =
        ibSurfaceWriter<Type>::New("vtk");

    const triSurface& ts = ibPatch_.ibMesh();

    // Make a face list for writing
    faceList f(ts.size());
    forAll (ts, faceI)
    {
        f[faceI] = ts[faceI].triFaceFace();
    }

    writerPtr->write
    (
        this->dimensionedInternalField().path(),
        ibPatch_.name(),
        ts.points(),
        f,
        this->dimensionedInternalField().name(),
        this->triValue() 
    );
 */
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
