Info<<"Start IB cells correction..."<< endl;    

int iter =0;
scalar number=0;

do
{

	forAll (mesh.boundary(), patchI)
	{
	    if (isA<immersedBoundaryFvPatch>(mesh.boundary()[patchI]))
	    {
	        const immersedBoundaryFvPatch& ibPatch =
	            refCast<const immersedBoundaryFvPatch>
	            (
	                mesh.boundary()[patchI]
	            );

            //#include "test.H"
	        const scalar Cmu25 = Foam::pow(0.09, 0.25);
	
            const vectorField& n = ibPatch.ibNormals();
		    // Velocity
		    const fvPatchVectorField& Uwg =
      			ibPatch.lookupPatchField<volVectorField, vector>("U");
    		const immersedBoundaryVelocityWallFunctionFvPatchVectorField& Uw =
        		refCast<const immersedBoundaryVelocityWallFunctionFvPatchVectorField>
        		(
            		Uwg
        		);

		    // Calculate tangential component, taking into account wall velocity
    		const vectorField UtanOld =
        		(I - sqr(n)) & (Uw.ibSamplingPointValue() - Uw.ibValue());
		    const scalarField magUtanOld = mag(UtanOld);

        	// Laminar viscosity
        	const fvPatchScalarField& nuwg =
        	    ibPatch.lookupPatchField<volScalarField, scalar>("nu");
        	const immersedBoundaryFvPatchScalarField& nuw =
        	    refCast<const immersedBoundaryFvPatchScalarField>(nuwg);
        	scalarField nu = nuw.ibCellValue();

        	// Turbulence kinetic energy
        	const fvPatchScalarField& kg =
        	    ibPatch.lookupPatchField<volScalarField, scalar>("k");
        	const immersedBoundaryWallFunctionFvPatchScalarField& kw =
        	    refCast<const immersedBoundaryWallFunctionFvPatchScalarField>(kg);

        	const labelList& ibc=ibPatch.ibCells();

        	scalarField k = kw.ibSamplingPointValue();

        	// Distance to sampling point
        	const scalarField& ySample = ibPatch.ibSamplingPointDelta();

        	// Calculate yPlus for sample points
        	scalarField ypd = Cmu25*ySample*sqrt(k)/nu;

        	labelHashSet ibCellsToChangeListSet;

        	forAll (ibc, cellI)
	        {
    	        label ibCellI=ibc[cellI];

    	        if(ypd[cellI]<22)
    	        {
    	            ibCellsToChangeListSet.insert(ibCellI);

    	        }
    	    }

    	    labelList ibCellsToChangeList=ibCellsToChangeListSet.toc();

			number=returnReduce(ibCellsToChangeList.size(), sumOp<label>());			
			
    	    Info<<"Iter. "<<iter++<<" IB cells correction number: "<< number<< endl;
	        
			//if(ibCellsToChangeList.size()>0)
			{
	            ibPatch.addAdjacentCells(ibCellsToChangeList);
	        }
	        cellIbMask = ibPatch.gamma();
	        cellIbMaskExt = ibPatch.gammaExt();
	        faceIbMask = ibPatch.sGamma();
	        const labelList& ibCells =  ibPatch.ibCells();
	        const labelListList& ibCellCells =  ibPatch.ibCellCells();
	        ibCellMask = 0;
	
	        forAll(ibCells, I)
	        {
	            ibCellMask.internalField()[ibCells[I]]=1.0;

	            forAll(ibCellCells[I],II)
	            {
	                samplingCellMask[ibCellCells[I][II]]=1.0;
	            }
	        }

	    }
	}
	
    forAll (cellIbMask.boundaryField(), patchI)
    {
        if (cellIbMask.boundaryField()[patchI].coupled())
        {
            cellIbMask.boundaryField()[patchI].initEvaluate
                (
                    Pstream::blocking
                );
            cellIbMaskExt.boundaryField()[patchI].initEvaluate
                (
                    Pstream::blocking
                );
        }

        if (cellIbMask.boundaryField()[patchI].coupled())
        {
            cellIbMask.boundaryField()[patchI].evaluate
                (
                    Pstream::blocking
                );
            cellIbMaskExt.boundaryField()[patchI].evaluate
                (
                    Pstream::blocking
                );
        }
    }

}while (number>0);

