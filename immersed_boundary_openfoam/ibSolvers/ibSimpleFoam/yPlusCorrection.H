
if(correctYPlus)
{

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

	        const scalar Cmu25 = Foam::pow(0.09, 0.25);
	/*
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
		    const scalarField& magUtanOld = mag(UtanOld);*/

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

        	const labelList& ibc = ibPatch.ibCells();

        	const scalarField k = kw.ibSamplingPointValue();
        	//const scalarField k = kw.ibCellValue();

        	// Distance to sampling point
        	//const scalarField& ySample = ibPatch.ibSamplingPointDelta();
        	const scalarField& ySample = ibPatch.ibDelta();
 
        	// Calculate yPlus for sample points
    	    tmp<scalarField > typd(new scalarField(Cmu25*ySample*sqrt(k)/nu));
 
	        scalarField& ypd = typd();

        	labelHashSet ibCellsToChangeListSet;

        	forAll (ibc, cellI)
	        {
    	        label ibCellI=ibc[cellI];
                //Info<<ypd[cellI]<<endl;
    	        if(ypd[cellI]<20)
    	        {
//    	            Info<<cellI<<tab<<ySample[cellI]<<endl;
    	            ibCellsToChangeListSet.insert(ibCellI);

    	        }
    	    }

    	    labelList ibCellsToChangeList = ibCellsToChangeListSet.toc();

			number = returnReduce(ibCellsToChangeList.size(), sumOp<label>());			
			
    	    Info<<"Iter. "<<iter++<<" IB cells correction number: "<< number<< endl;
	        
			if(number>0)
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

}

