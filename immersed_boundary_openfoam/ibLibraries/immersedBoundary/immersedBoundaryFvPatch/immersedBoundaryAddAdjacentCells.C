
void Foam::immersedBoundaryFvPatch::addAdjacentCells(const labelList& ibCellsToChangeList) const
{
    if (debug)
    {
        InfoIn("void immersedBoundaryFvPatch::addAdjacentCells() const")
            << "add adjacent cells next to ib cells whose yPlus is too small"
            << "for immersed boundary " << name()
            << endl;
    }

    //const vectorField& C = mesh_.cellCentres();
    //const vectorField& Cf = mesh_.faceCentres();
    //const vectorField& Sf = mesh_.faceAreas();

    //const cellList& meshCells = mesh_.cells();
    //const unallocLabelList& own = mesh_.owner();
    //const unallocLabelList& nei = mesh_.neighbour();

    const labelListList& meshCellCells = mesh_.cellCells();
    const cellList& meshCellFaces = mesh_.cells();

    const labelList& ibc = ibCells();

	volScalarField tmpgIn = gamma(); 
	tmpgIn.clear();// to make sure gammaPtr_ is initialized
	volScalarField tmpgExt = gammaExt(); 
	tmpgExt.clear();// to make sure gammaPtr_ is initialized
    scalarField& gIn = gammaPtr_->internalField();
    scalarField& gExtIn = gammaExtPtr_->internalField();


/*
    labelListList procFacesToChange(Pstream::nProcs());

	if (Pstream::parRun())
	{
		procFacesToChange[Pstream::myProcNo()].setSize(0);
   
        Pstream::gatherList(procFacesToChange);
        Pstream::scatterList(procFacesToChange);
	}
*/

    labelHashSet ibCellSet;

    labelHashSet ibCellOldSet;

    labelHashSet addAdjacentCellSet;

    forAll (ibc, cellI)
    {
        ibCellSet.insert(ibc[cellI]);
        ibCellOldSet.insert(ibc[cellI]);
    }

    forAll (ibCellsToChangeList, cellI)
    {
        //the id label of the IB cell which needs to be replaced by adjacent cell(s)
        label ibCellID = ibCellsToChangeList[cellI];

        //gIn[ibCellID] = 0;

        ibCellOldSet.unset(ibCellID);

        //neigbour cells of the IB cell
        labelList neiCells = meshCellCells[ibCellID];   

        forAll (neiCells, cellII)
        {
            label neiCellI = neiCells[cellII];            
         
            // if the neigbour cell is not IB cells & it is live (fluid cell)
            if(!ibCellSet.found(neiCellI)&& !addAdjacentCellSet.found(neiCellI) && gIn[neiCellI]>SMALL)
            {
                addAdjacentCellSet.insert(neiCellI);
                gIn[neiCellI] = 0;
            }
        }

        //gIn[ibCellID] = 0;

//        gExtIn[ibCellID] = 0;

	    if (Pstream::parRun())// find ibCells next to coupled patch
    	{   
			volScalarField::GeometricBoundaryField& gammaPatches = gammaPtr_->boundaryField();

        	cell neiFaces = meshCellFaces[ibCellID];		 

			forAll (neiFaces, faceI)
            {
                label curFace=neiFaces[faceI];
				label patchi = mesh_.boundaryMesh().whichPatch(curFace);
                if(patchi>-1)
				{	
					if(gammaPatches[patchi].coupled())
	                {
            			scalarField gammaOwn =
                			gammaPatches[patchi].patchInternalField();

            			scalarField gammaNei =
                			gammaPatches[patchi].patchNeighbourField();

						const label start = mesh_.boundaryMesh()[patchi].start();

						if(gammaNei[curFace-start]>SMALL)	
						{
							gIn[ibCellID] = -2;  
						}						
					}
				}
			}
		} 
    }    

	if (Pstream::parRun())
	{
		forAll (mesh_.boundary(), patchi)
		{
			if(mesh_.boundary()[patchi].coupled())
			{
				volScalarField::GeometricBoundaryField& gammaPatches = gammaPtr_->boundaryField();

            	gammaPatches[patchi].initEvaluate(Pstream::blocking);// not quite sure
            	gammaPatches[patchi].evaluate(Pstream::blocking);// not quite sure

	            scalarField gammaOwn =
	                gammaPatches[patchi].patchInternalField();

	            scalarField gammaNei =
	                gammaPatches[patchi].patchNeighbourField();

				forAll(gammaOwn, facei)
				{					
					if(gammaNei[facei]<-1 and gammaOwn[facei]>0)
					{
						label cellI=mesh_.boundary()[patchi].faceCells()[facei];

						addAdjacentCellSet.insert(cellI);

						gIn[cellI] = 0;
					}
				}


			}
		}
	}

    forAll (ibCellsToChangeList, cellI)
    {
        //the id label of the IB cell which needs to be replaced by adjacent cell(s)
        label ibCellID = ibCellsToChangeList[cellI];

        gIn[ibCellID] = 0;

        gExtIn[ibCellID] = 0;
  
    }  

    //labelList IbCellOld = ibCellOldSet.toc();
    labelList addAdjacentCells = addAdjacentCellSet.toc();

    ibCellsPtr_ = new labelList(ibCellOldSet.toc());
    //Info << "After yPlus correction, number of IB cells: " << returnReduce(ibCellsPtr_->size(), sumOp<label>())<< endl;

    ibCellsPtr_->append(addAdjacentCells);

    sort(*ibCellsPtr_);

    Info << "After yPlus correction, number of IB cells: " << returnReduce(ibCellsPtr_->size(), sumOp<label>())<< endl;

    // Get non-const reference to patch
    immersedBoundaryPolyPatch& ref_ibPolyPatch =
        const_cast<immersedBoundaryPolyPatch&>(ibPolyPatch_);
    ref_ibPolyPatch.moveIb();

    //deleteDemandDrivenData(gammaPtr_);

    deleteDemandDrivenData(ibFacesPtr_);
    deleteDemandDrivenData(ibFaceCellsPtr_);
    deleteDemandDrivenData(ibFaceFlipsPtr_);
    deleteDemandDrivenData(ibInsideFacesPtr_);
    deleteDemandDrivenData(ibInternalFacesPtr_);
    deleteDemandDrivenData(ibPointsPtr_);
    deleteDemandDrivenData(ibNormalsPtr_);
    deleteDemandDrivenData(hitFacesPtr_);
    deleteDemandDrivenData(ibSamplingPointsPtr_);

    deleteDemandDrivenData(ibSamplingWeightsPtr_);
    deleteDemandDrivenData(ibSamplingProcWeightsPtr_);

    deleteDemandDrivenData(cellsToTriAddrPtr_);
    deleteDemandDrivenData(cellsToTriWeightsPtr_);

    deleteDemandDrivenData(ibCellCellsPtr_);
    deleteDemandDrivenData(ibProcCellsPtr_);
    deleteDemandDrivenData(ibProcIbCellsPtr_);
    deleteDemandDrivenData(ibProcCentresPtr_);
    deleteDemandDrivenData(ibProcGammaPtr_);
    deleteDemandDrivenData(ibCellProcCellsPtr_);
    deleteDemandDrivenData(deadCellsPtr_);
    deleteDemandDrivenData(deadCellsExtPtr_);
    deleteDemandDrivenData(deadFacesPtr_);
    deleteDemandDrivenData(liveCellsPtr_);
    deleteDemandDrivenData(ibCellSizesPtr_);

    deleteDemandDrivenData(invDirichletMatricesPtr_);
    deleteDemandDrivenData(invNeumannMatricesPtr_);

    deleteDemandDrivenData(ibSfPtr_);
    deleteDemandDrivenData(ibMagSfPtr_);
    deleteDemandDrivenData(ibDeltaPtr_);
    deleteDemandDrivenData(ibSamplingPointDeltaPtr_);

    deleteDemandDrivenData(adjacentIbCellsPtr_);
    deleteDemandDrivenData(adjacentIbPointsPtr_);
    deleteDemandDrivenData(adjacentIbNormalsPtr_);
    deleteDemandDrivenData(adjacentIbDeltaPtr_);
    deleteDemandDrivenData(ibCellPtsPtr_);
    deleteDemandDrivenData(ibCellProcPtsPtr_);
    deleteDemandDrivenData(ibProcPtsPtr_);
    deleteDemandDrivenData(ibNewSamplingWeightsPtr_);
    deleteDemandDrivenData(ibNewSamplingProcWeightsPtr_);
    deleteDemandDrivenData(ibNewPtsWeightsPtr_);
    deleteDemandDrivenData(ibNewPtsProcWeightsPtr_);

    deleteDemandDrivenData(triFacesToTriPointsWeightsPtr_);

    deleteDemandDrivenData(triSfPtr_);
    deleteDemandDrivenData(triNormalsPtr_);
    deleteDemandDrivenData(triSlopesPtr_);
	deleteDemandDrivenData(triFacesInMeshPtr_);
	deleteDemandDrivenData(triFacesInMeshTypePtr_);
	deleteDemandDrivenData(triFacesInMeshTypeVecPtr_);
 
	deleteDemandDrivenData(triFacesToTriEdgesWeightsPtr_);
}
