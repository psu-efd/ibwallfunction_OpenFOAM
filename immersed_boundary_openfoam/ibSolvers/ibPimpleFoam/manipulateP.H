forAll (mesh.boundary(), patchI)
    {
        if (isA<immersedBoundaryFvPatch>(mesh.boundary()[patchI]))
        {

            const immersedBoundaryFvPatch& ibPatch_ =
                refCast<const immersedBoundaryFvPatch>
                (
                    mesh.boundary()[patchI]
                );    
 
			//Info<<"manipulate pEqn"<<endl;
 	
			fvMatrix<scalar>& eqn = pEqn;
			fvMesh& mesh_=mesh;
	
		    // Calculate gradient contribution
    		const labelList& ibFaces = ibPatch_.ibFaces();
    		const labelList& ibFaceCells = ibPatch_.ibFaceCells();
   			//const labelList& ibCells = ibPatch_.ibCells();

	
    		//Info<<"correctOffDiag "<< ibFaces.size()<<" "<<ibCells.size()<<endl;
	
    		const scalarField& ibGamma = ibPatch_.gamma().internalField();

    		const unallocLabelList& own = mesh_.owner();
    		const unallocLabelList& nei = mesh_.neighbour();

    		// Get delta coefficients
    		const surfaceScalarField& dc = mesh_.deltaCoeffs();
    		const scalarField& dcI = dc.internalField();
    		//const scalarField& ibD = ibPatch_.ibDelta();

    		const surfaceVectorField& Sf = mesh_.Sf();
		    const surfaceScalarField& magSf = mesh_.magSf();
 

		    // assign new values of ibFaceGrad
		    tmp<vectorField > tibFaceGrad
    		(
        		//new vectorField(ibCells.size(), pTraits<vector>::zero)
		        new vectorField(ibFaces.size(), pTraits<vector>::zero)
    		);

		    vectorField& ibFaceGrad = tibFaceGrad();

 

  
        	//Info<<name_<<endl;
         
	        volScalarField& psiI1 = const_cast<volScalarField&>
    	        (mesh_.lookupObject<volScalarField>("p"));
        
    	    // calcualte the gradient field
    	    //const GeometricField<vector,fvPatchField,volMesh>& gradField = fvc::grad(psiI1);

    	    // calculate the gradient on the internal faces
//    	    const vectorField& gradFieldSI = linearInterpolate(fvc::grad(psiI1));

    	    // insert the gradient on ibFaces to ibFaceGrad
    	    //ibFaceGrad.map(linearInterpolate(fvc::grad(psiI1)),ibFaces);  
    	    forAll (ibFaces, faceI)
    	    {
    	        const label curFace = ibFaces[faceI];

    	        if (curFace < nei.size())
    	        {
    	            //ibFaceGrad[faceI] = gradFieldSI[curFace];
    	        }
    	        //else
    	        //{
    	        //    label patchi = mesh_.boundaryMesh().whichPatch(curFace);

    	            //label patchFacei =
                  //  mesh_.boundaryMesh()[patchi].whichFace(curFace);
 
	//                ibFaceGrad[faceI] = gradField.boundaryField()[patchi][patchFacei];
 		//           }
    		    }
 
    	    scalarField& diag = eqn.diag();

    	    scalarField& upper = eqn.upper();
    	    Field<scalar>& source = eqn.source();

    	    forAll (ibFaces, faceI)
    	    {
    	        const label curFace = ibFaces[faceI];
    
	            if (curFace < nei.size())
	            {

                // calculate surface vector on ibFaces    added by Xu 7/2017
                const vector SfN = Sf[curFace]/magSf[curFace];

                // calculate value gradient on ibFaces o->n     added by Xu 7/2017
                const scalar faceGrad = (SfN&ibFaceGrad[faceI]) * pTraits<scalar>::one;
 
                // Internal face.  One side is an ibCell and another is a
                // live cell. Add gradient to the source of the live cell
                // and kill the off-diagonal coefficient
                if (ibGamma[own[curFace]] > SMALL)
                {
                    diag[own[curFace]] += upper[curFace];


                    source[own[curFace]] -=
                        upper[curFace]*faceGrad/dcI[curFace];

                }
                else
                {
                    diag[nei[curFace]] += upper[curFace];

                    source[nei[curFace]] +=
                        upper[curFace]*faceGrad/dcI[curFace];

                }


                upper[curFace] = 0.0;
            }
            else
            {
 
                // else MPH
                label patchi = mesh_.boundaryMesh().whichPatch(curFace);

                if (!eqn.internalCoeffs()[patchi].empty())
                {
                    label patchFacei =
                        mesh_.boundaryMesh()[patchi].whichFace(curFace);

                    // calculate value gradient on ibFaces o->n     added by Xu 8/2017
                    const scalarField& deltaCoeffs = mesh_.boundary()[patchi].deltaCoeffs();

                    tmp<scalarField> snGrad = psiI1.boundaryField()[patchi].snGrad(deltaCoeffs);

                    const scalar faceGrad = snGrad()[patchFacei];

                    eqn.internalCoeffs()[patchi][patchFacei] =
                        pTraits<scalar>::zero;

                    eqn.boundaryCoeffs()[patchi][patchFacei] =
                        pTraits<scalar>::zero;

                    // Check if the live cell is on local or neighbour side
                    // HJ, 7/Dec/2012
                    if (ibFaceCells[faceI] > -1)
                    {
                        if (ibGamma[ibFaceCells[faceI]] > SMALL)
                        {
                            source[ibFaceCells[faceI]] +=
                                faceGrad
                                /dc.boundaryField()[patchi][patchFacei];
                        }
                    }
                }
            }
        }

 

        Field<scalar> polyPsi(pEqn.psi(), ibPatch_.ibCells());

        //pEqn.setValues(ibPatch_.ibCells(), polyPsi);


        Field<scalar> deadCellsPsi
        (
            pEqn.psi(),
            ibPatch_.deadCells()
        ); 
        //Info<<polyPsi[11]<<endl;
        pEqn.setValues(ibPatch_.deadCells(), deadCellsPsi);
        Info<<"Finished manipulate pEqn"<<endl;
        }
    }

