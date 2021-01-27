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

\*---------------------------------------------------------------------------*/

#include "immersedBoundaryFixedFluxPressure.H"
#include "fvPatchFieldMapper.H"
#include "fvMatrix.H"
#include "surfaceWriter.H"
#include "addToRunTimeSelectionTable.H"

#include "triSurfaceGeoMesh.H"
#include "triSurfaceFieldToTecplot.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
//    This member function only update the values, grandient, and other
//    information related to immersed boundary. It does NOT change any
//    internal and boundary value of the volume field. This part is 
//    designed so it does not interfer with the boundary effect treatment
//    of the original OF (i.e., correctBoundaryConditions(), updateCoeff(),
//    evaluate(), and any other callers such as linear solver and matrix
//    construction). 

//
//    First, the values at global image points are calculated,
//    then, the values at global ghost cells are calculated.
//    Finally, the values at local ghost cells (if any) are just a slices of
//    of the global ghost cell values.
void immersedBoundaryFixedFluxPressure::updateGhostCellValues() const
{
    if(ibPatch_.local_ghostCells().size() == 0)
    {
        //Even if there is no ghost cell locally on current processor,
        //it still needs to work since it should send information needed by others.
//        Pout << "There is no local ghost cell on this processor." << endl;
        hasGhostCells_ = false;
    }

    // Create local image point values
    Field<Type> local_imagePointValues
    (
        ibPatch_.global_ghostCells().size(),
        pTraits<Type>::zero 
    );

    // Create global image point values
    Field<Type> global_imagePointValues
    (
        ibPatch_.global_ghostCells().size(),
        pTraits<Type>::zero
    );

    //get the surrounding cells (in local cell ID) for all image points
    const labelListList& local_surroundingCells =
               ibPatch_.local_surroundingCells();

    //get the interpolation weights for all the image points
    const PtrList<scalarField>& interp_Weight =  ibPatch_.interp_Weight();

    //get the image points
    const vectorField& local_ipts = ibPatch_.local_ipts();

    // Get image point values through interpolation
    // If parallel, each processor needs to provide the values
    // of local interpolation cells.
    if(Pstream::parRun()) 
    {
       //get the list of global image points
       const vectorField& global_ipts = ibPatch_.global_ipts();

       //get the clean, non-duplicate list of surrounding cells on
       //on this processor (in local cell ID)
       const labelList& local_interpolateCellsList = 
                     ibPatch_.local_interpolateCells(); 

       //number of local interpolate cells
       label nLocalInterpCells = local_interpolateCellsList.size();

       //get the global interpolate cell list (in global cell ID)
       const labelList& global_interpolateCellsList =
                     ibPatch_.global_interpolateCellsList();

       //get the mapping from local interpolate cell to its
       //position in the global interpolate cell list
       const labelList& local_to_global_interpolateCellsMap =
                     ibPatch_.local_to_global_interpolateCellsMap();

       //get the mapping from each image point's surrounding cells to 
       //their positions in the golbal_interpolateCellsList
       const PtrList<labelList>& surroundCells_To_InterpCells_map = 
                     ibPatch_.surroundCells_To_InterpCells_map();

       //field value for all global_interpolateCellsList
       Field<Type> global_interpolateCells_Field
            (global_interpolateCellsList.size(), pTraits<Type>::zero);

       //every processor do its share locally
       //loop through all local interpolate cells
       for(int cellI=0;cellI<nLocalInterpCells;cellI++)
       {
            //current interpolate cell's local cell ID
            label localCellID = local_interpolateCellsList[cellI];

            //current interpolate cell's position in the 
            //global interpolate cell list
            label globalListPosition = 
                    local_to_global_interpolateCellsMap[cellI];

            //fill the value to the position in the global value list
            global_interpolateCells_Field[globalListPosition]=
                      this->internalField()[localCellID];
       }

       if(debug)
       {
//            Pout << "global_interpolateCells_Field before combine = "
//                 << global_interpolateCells_Field << endl;
       }

       //then combine them together
       reduce(global_interpolateCells_Field, sumOp<Field<Type> >());

       if(debug)
       {
//           Pout << "global_interpolateCells_Field after combine = "
//                << global_interpolateCells_Field << endl;
       }
       //end of gathering global interpolate cell data

       //Now since we have the value of all global interpolateCells,
       //each processor can start to work. Here, each processor
       //will interpolate for ALL (global) image points and then get
       //its local slice. It is a bit of waste. Ideally, each processor
       //only needs to interpolate to its local image points. But that
       //requires one more mapping to get its local image points.

       //get the interpolation weights for all (global) image points
       const PtrList<scalarField>& interp_Weight =  ibPatch_.interp_Weight();
       const labelListList& global_surroundingCells = 
                              ibPatch_.global_surroundingCells();

       if(debug)
       {
//          Pout << "interp_Weight = " << interp_Weight.first() << endl;
       }

       //do interpolation to ALL (global) image points
       forAll(global_ipts, iptsI)
        {
            const labelList& tmpCellList = global_surroundingCells[iptsI];
            const labelList& tmp_surroundCells_To_InterpCells_map = 
                              surroundCells_To_InterpCells_map[iptsI];
            const scalarField& tmpWeight = interp_Weight[iptsI];


            //loop over current image point's surrounding cells
            forAll(tmpCellList, cellI)
            {
                //find current surrounding cell's position in the
                //global_interpolateCells list
                label position = tmp_surroundCells_To_InterpCells_map[cellI];

                //find the value of the surrounding cells and do average
                global_imagePointValues[iptsI] += tmpWeight[cellI]*
                           global_interpolateCells_Field[position];

                if(debug)
                {
/*
                   Pout << "iptsI = " << iptsI << " "
                        << "image point = " << global_ipts[iptsI] << " " 
                        << "cellI = " << cellI << " "
                        << "position = " << position << " " 
                        << "tmpWeight[cellI] = " << tmpWeight[cellI] << " "
                        << "global_interpolateCells_Field[position] = " 
                        << global_interpolateCells_Field[position] << endl;
*/

                }
            }
        }

   }
    else //Serial run
    {
       //do interpolation to the image point
       forAll(local_ipts, iptsI)
        {
            const labelList& tmpCellList = local_surroundingCells[iptsI];
            const scalarField& tmpWeight = interp_Weight[iptsI];
            

            forAll(tmpCellList, cellI)
            {
                if(debug)
                {
/*
                   Info << "iptsI = " << iptsI << " "
                        << "image point = " << local_ipts[iptsI] << " " 
                        << "cellI = " << cellI << " " 
                        << "tmpWeight[cellI] = " << tmpWeight[cellI] << " " 
                        << "this->internalField()[tmpCellList[cellI]] = " 
                        << this->internalField()[tmpCellList[cellI]] << endl;
*/
                }

                //find the value of the surrounding cells and do average
                global_imagePointValues[iptsI] += tmpWeight[cellI]*
                           this->internalField()[tmpCellList[cellI]];
            }
        }
    }//End of serial run

    //get the location code of global image points
    const labelList& global_imageLocation = ibPatch_.global_ipLocation();

    //if the global ghost cell's image point is not inside a
    //fluid cell, then treat the ghost cell as dead cell
/*    forAll(global_imagePointValues, iptsI)
    {
        if(global_imageLocation[iptsI] != 1)
        {
           global_imagePointValues[iptsI] = deadCellValue_;
        }
    }
*/
    if (this->fixesValue())
    {
//        Pout << "fixes value" << endl;
 
        //value at the ghost cell hit points 
        ibValue_ = ibPatch_.toGhostCellHitPoints(refValue_);

        //values at the global ghost cells
        Field<Type> ibGlobal_GhostCellValues = 2.0*ibValue_ - global_imagePointValues;

        //now get the value at local ghost cells (if any)
        const labelList& local_to_global_iptsMap =
                           ibPatch_.local_to_global_iptsMap();
 
        //loop over all local image points (ghost cells, if any)
        ibGhostCellValues_.setSize(local_ipts.size());
        forAll(local_ipts, iptsI)
        {
           //get the position of current local impage point (ghost cell) in
           //the global image point (ghost cell) list
           label position = local_to_global_iptsMap[iptsI];

           ibGhostCellValues_[iptsI] = ibGlobal_GhostCellValues[position];
        }

        //gradient at the ghost celll hit points
        ibGrad_ = (global_imagePointValues - ibValue_)/
                     ibPatch_.global_ghostCellsDistance();

        if(debug)
        {
/*
           Info << "ib fixes value ..." << endl;
           Info << "ibGhostCellValues_ = " << ibGhostCellValues_ << endl;
           Info << "ibValue_ = " << ibValue_ << endl;
           Info << "imagePointValues = " << global_imagePointValues << endl;
*/
        }
    }
    else
    {
//        Pout << "fixes gradient " << endl;

        //gradient at the ghost cell hit points
        ibGrad_ = ibPatch_.toGhostCellHitPoints(refGrad_);
//Pout << "ibGrad_ = " << ibGrad_ << endl;

        //value at the ghost cell hit points
        ibValue_ = global_imagePointValues - 
                     ibGrad_*ibPatch_.global_ghostCellsDistance();
//Pout << "ibValue_ = " << ibValue_ << endl;

        //values at the global ghost cells
        Field<Type> ibGlobal_GhostCellValues = 2.0*ibValue_ - global_imagePointValues;
//Pout << "ibGlobal_GhostCellValues = " << ibGlobal_GhostCellValues << endl;

        //now get the value at local ghost cells (if any)
        const labelList& local_to_global_iptsMap =
                           ibPatch_.local_to_global_iptsMap();

        //loop over all local image points (ghost cells, if any)
        ibGhostCellValues_.setSize(local_ipts.size());
        forAll(local_ipts, iptsI)
        {
           //get the position of current local impage point (ghost cell) in
           //the global image point (ghost cell) list
           label position = local_to_global_iptsMap[iptsI];

           ibGhostCellValues_[iptsI] = ibGlobal_GhostCellValues[position];
        }
    }
//Pout << "ibGhostCellValues_ = " << ibGhostCellValues_ << endl;
}

// set the ghost cell values
void immersedBoundaryFixedFluxPressure::imposeGhostCellsCondition() const
{
    const labelList& local_ghostCells = ibPatch_.local_ghostCells();

    // Get non-const access to internal field
    Field<scalar>& psiI = const_cast<Field<scalar>&>(this->internalField());

    const Field<scalar>& gcValues = ibGhostCellValues();
   
    forAll (local_ghostCells, dcI)
    {
        psiI[local_ghostCells[dcI]] = gcValues[dcI];
    }
}


// set the dead (totally inside) cell values
void immersedBoundaryFixedFluxPressure::imposeFullyInsideCellsCondition() const
{
    const labelList& local_fullyInsideCells = ibPatch_.local_fullyInsideCells();

    // Get non-const access to internal field
    Field<scalar>& psiI = const_cast<Field<scalar>&>(this->internalField());

    forAll (local_fullyInsideCells, dcI)
    {
        psiI[local_fullyInsideCells[dcI]] = deadCellValue_;
    }
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

immersedBoundaryFixedFluxPressure::immersedBoundaryFixedFluxPressure
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvPatchField<scalar>(p, iF, Field<scalar>(0)),
    ibPatch_(refCast<const immersedBoundaryFvPatch>(p)),
    mesh_(p.boundaryMesh().mesh()),
    refValue_(ibPatch_.ibMesh().size(), 0.0),
    refGrad_(ibPatch_.ibMesh().size(), 0.0),
    fixesValue_(false),
    setDeadCellValue_(false),
    deadCellValue_(0.0),
    ibValue_(),
    ibGrad_(),
    hasGhostCells_(true),
    ibGhostCellValues_(),
    phiHbyAName_("phiHbyA"),
    phiName_("phi"),
    rhoName_("rho"),
    DpName_("Dp"),
    adjoint_(false)
{  

    if(debug) 
    {
//      Info << "Construct immersedBoundaryFixedFluxPressure from default." << endl;
    }
}


immersedBoundaryFixedFluxPressure::immersedBoundaryFixedFluxPressure
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<scalar>(p, iF, Field<scalar>(0)),
    ibPatch_(refCast<const immersedBoundaryFvPatch>(p)),
    mesh_(p.boundaryMesh().mesh()),
    refValue_("refValue", dict, ibPatch_.ibMesh().size()),
    refGrad_("refGradient", dict, ibPatch_.ibMesh().size()),
    fixesValue_(dict.lookup("fixesValue")),
    setDeadCellValue_(dict.lookup("setDeadCellValue")),
    deadCellValue_(dict.lookup("deadCellValue"))),
    ibValue_(),
    ibGrad_(),
    hasGhostCells_(true),
    ibGhostCellValues_(),
    phiHbyAName_(dict.lookupOrDefault<word>("phiHbyA", "phiHbyA")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    DpName_(dict.lookupOrDefault<word>("Dp", "Dp")),
    adjoint_(dict.lookupOrDefault<Switch>("adjoint", false))
{
    if(debug)  
    {
//        Info << "Construction of immersedBoundaryFixedFluxPressure from dictionary ... " <<endl;
    }

    if (!isType<immersedBoundaryFvPatch>(p))
    {
        FatalIOErrorIn
        (
            "immersedBoundaryFixedFluxPressure<Type>::"
            "immersedBoundaryFixedFluxPressure\n"
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


immersedBoundaryFixedFluxPressure::immersedBoundaryFixedFluxPressure
(
    const immersedBoundaryFixedFluxPressure& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& m
)
:
    fvPatchField<scalar>(p, iF, Field<scalar>(0)),
    ibPatch_(refCast<const immersedBoundaryFvPatch>(p)),
    mesh_(p.boundaryMesh().mesh()),
    refValue_(ptf.refValue()),
    refGrad_(ptf.refGrad()),
    fixesValue_(ptf.fixesValue()),
    setDeadCellValue_(ptf.setDeadCellValue_),
    deadCellValue_(ptf.deadCellValue_),
    ibValue_(),
    ibGrad_(),
    hasGhostCells_(true),
    ibGhostCellValues_(),
    phiHbyAName_(ptf.phiHbyAName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    DpName_(ptf.DpName_),
    adjoint_(ptf.adjoint_)
{
    if(debug)
    {
//        Info << "Construction of immersedBoundaryFixedFluxPressure by mapping " << endl;
    }

    // Note: NO MAPPING.  Fields are created on the immersed boundary
    // HJ, 12/Apr/2012
    if (!isType<immersedBoundaryFvPatch>(p))
    {
        FatalErrorIn
        (
            "immersedBoundaryFixedFluxPressure<Type>::"
            "immersedBoundaryFixedFluxPressure\n"
            "(\n"
            "    const immersedBoundaryFixedFluxPressure<Type>&,\n"
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


immersedBoundaryFixedFluxPressure::immersedBoundaryFixedFluxPressure
(
    const immersedBoundaryFixedFluxPressure& ptf
)
:
    fvPatchField<scalar>
    (
        ptf.patch(),
        ptf.dimensionedInternalField(),
        Field<scalar>(0)
    ),
    ibPatch_(ptf.ibPatch()),
    mesh_(ptf.patch().boundaryMesh().mesh()),
    refValue_(ptf.refValue()),
    refGrad_(ptf.refGrad()),
    fixesValue_(ptf.fixesValue()),
    setDeadCellValue_(ptf.setDeadCellValue_),
    deadCellValue_(ptf.deadCellValue_),
    ibValue_(),
    ibGrad_(),
    hasGhostCells_(true),
    ibGhostCellValues_(),
    phiHbyAName_(ptf.phiHbyAName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    DpName_(ptf.DpName_),
    adjoint_(ptf.adjoint_)
{
  if(debug)
  {
//     Info << "Construct immersedBoundaryFixedFluxPressure from a copy" << endl;
  }
}


immersedBoundaryFixedFluxPressure::immersedBoundaryFixedFluxPressure
(
    const immersedBoundaryFixedFluxPressure& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvPatchField<scalar>(ptf.patch(), iF, Field<scalar>(0)),
    ibPatch_(ptf.ibPatch()),
    mesh_(ptf.patch().boundaryMesh().mesh()),
    refValue_(ptf.refValue()),
    refGrad_(ptf.refGrad()),
    fixesValue_(ptf.fixesValue()),
    setDeadCellValue_(ptf.setDeadCellValue_),
    deadCellValue_(ptf.deadCellValue_),
    ibValue_(),
    ibGrad_(),
    ibGhostCellValues_(),
    phiHbyAName_(ptf.phiHbyAName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    DpName_(ptf.DpName_),
    adjoint_(ptf.adjoint_)
{
  if(debug)
  {
//     Info << "Construct immersedBoundaryFixedFluxPressure from a copy and internal field" << endl;
  }
}

//*******************************Member Operators****************************//
void Foam::immersedBoundaryFixedFluxPressure::copyBCParameters
                   (const immersedBoundaryFixedFluxPressure& ptf) const
{
    //Assign the parameter values for immersedBoundaryFixedFluxPressure.
    refValue_ = ptf.refValue();
    refGrad_ = ptf.refGrad();
    fixesValue_ = ptf.fixesValue();
    setDeadCellValue_ = ptf.setDeadCellValue_;
    deadCellValue_ = ptf.deadCellValue_;

    phiHbyAName_ = ptf.phiHbyAName_;
    phiName_ = ptf.phiName_;
    rhoName_ = ptf.rhoName_;
    DpName_ = ptf.DpName_;
    adjoint_ = ptf.adjoint_;

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Field<scalar>& immersedBoundaryFixedFluxPressure::ibValue() const
{
//    if (ibValue_.empty() && hasGhostCells_)
//    if (ibValue_.empty())
//    {
//        this->updateGhostCellValues();
//    }

    return ibValue_;
}


const Field<scalar>& immersedBoundaryFixedFluxPressure::ibGrad() const
{
//    if (ibGrad_.empty() && hasGhostCells_)
//    if (ibGrad_.empty())
//    {
//        this->updateGhostCellValues();
//    }

    return ibGrad_;
}


const Field<scalar>& immersedBoundaryFixedFluxPressure::ibGhostCellValues() const
{
//    if (ibGhostCellValues_.empty() && hasGhostCells_)
//    if (ibGhostCellValues_.empty())
//    {
//        this->updateGhostCellValues();
//    }

    return ibGhostCellValues_;
}


tmp<Field<scalar> > immersedBoundaryFixedFluxPressure::triValue() const
{
//    Pout << "requesting immersedBoundaryFixedFluxPressure::triValue()" << endl;
//    Pout << "asking for ibValue() " << ibValue() << endl;
    return ibPatch_.toTriFaces(this->ibValue());
}


tmp<Field<scalar> > immersedBoundaryFixedFluxPressure<scalar>::triGrad() const
{
    return ibPatch_.toTriFaces(this->ibGrad());
}

//- Write the triangular center value field to vtk for visualization
void immersedBoundaryFixedFluxPressure::writeTriValue(const word& varName) const
{
//    Pout << "requesting immersedBoundaryFixedFluxPressure::writeTriValue()" << endl;
//    Pout << "asking for ibValue() " << ibValue() << endl;
    Field<scalar> triangularValues = ibPatch_.toTriFaces(this->ibValue());

    //- triSurfaceMesh from triSurface
    const triSurfaceMesh& stlMesh = ibPatch_.ibMesh();

    DimensionedField<scalar, triSurfaceGeoMesh> stlField
    (
         IOobject
         (
          "stlField",
          mesh_.time().timeName(),
          mesh_.time(),
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
         ),
         stlMesh,
         dimensioned<scalar>("zero", dimless, 0.0)
     );

//Pout << "stlField = " << stlField << endl;
//Pout << "stlField.size() = " << stlField.size() << endl;
     // fill the values to the triangle's centres
     forAll(stlField, triI)
     {
        stlField[triI] = triangularValues[triI];
     }

     fileName tpPath(mesh_.time().path()/"triValues");   
 
//     triSurfaceScalarFieldToTecplot(stlMesh, stlField, tpPath, "testField");
     triSurfaceScalarFieldToVTK
           (stlMesh, stlField, tpPath, varName+mesh_.time().timeName());
}


void immersedBoundaryFixedFluxPressure::updateCoeffs()
{
//   Pout << "immersedBoundaryFixedFluxPressure::updateCoeff() " << endl;

//   imposeImmersedBoundaryEffect();
   updateGhostCellValues();

   fvPatchScalarField::updateCoeffs();
}


void immersedBoundaryFixedFluxPressure::initEvaluate
(
    const Pstream::commsTypes
)
{
//   Pout << "immersedBoundaryFixedFluxPressure::initEvaluate() " << endl;
}


void immersedBoundaryFixedFluxPressure::evaluate
(
    const Pstream::commsTypes
)
{
    // Note
    // Since the boundary condition is performed by data fitting with the
    // internal field, fitting must be performed both on updateCoeffs
    // and on evaluate (internal field has changed in the meantime).
//    Pout << "immersedBoundaryFixedFluxPressure::evaluate() " << endl;

//    this->updateCoeffs();
    updateGhostCellValues();

    fvPatchScalarField::evaluate();
}


//Impose the effect of immersed boundary
void immersedBoundaryFixedFluxPressure::imposeImmersedBoundaryEffect() const
{
//    Pout << "imposeImmersedBoundaryEffect ..." << endl;

    updateGhostCellValues();

    //impose the ghost cells
    imposeGhostCellsCondition();

    //impose the fully inside (dead) cells
    if (setDeadCellValue_)
    {
        imposeFullyInsideCellsCondition();
    }
}


void immersedBoundaryFixedFluxPressure::write(Ostream& os) const
{
    fvPatchScalarField::write(os);

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
    autoPtr<surfaceWriter> writerPtr =
        surfaceWriter::New("vtk");

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
        this->triValue()(),
        true
    );
*/
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        immersedBoundaryFixedFluxPressure
    );
}

// ************************************************************************* //
