  volVectorField bodyForce
    (
        IOobject
        (
            "bodyForce",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        //dimensionedVector("zero", dimVelocity/dimTime, vector::zero)
        dimensionedVector("zero", dimVelocity/dimTime, vector(1E-5,0,0))
    );


forAll (mesh.boundary(), patchI)
{
    if (isA<immersedBoundaryFvPatch>(mesh.boundary()[patchI]))
    {
        const immersedBoundaryFvPatch& ibPatch =
            refCast<const immersedBoundaryFvPatch>
            (
                mesh.boundary()[patchI]
            );

        const triSurface& triS_ = ibPatch.ibMesh();

        fvPatchVectorField& Up = U.boundaryField()[patchI];
        //update wall shear stress and Shields parameters
        tmp<vectorField > tshearStress(new vectorField);
        vectorField& wallShearStress = tshearStress();

        if(isA<immersedBoundaryVelocityWallFunctionFvPatchVectorField>(U.boundaryField()[patchI]))
        {
            // shearStress based on ibCells
            immersedBoundaryVelocityWallFunctionFvPatchVectorField& ibPatchField_U = 
                    refCast<immersedBoundaryVelocityWallFunctionFvPatchVectorField>(Up);

            //const vectorField& ibShearStress = ibPatchField_U.wallShearStress();
            wallShearStress = ibPatch.toTriFaces(ibPatchField_U.wallShearStress());

            vectorField wss=ibPatchField_U.wallShearStress();
            labelList ibc=ibPatch.ibCells();
            vectorField& bf=bodyForce.internalField();

            forAll(ibc, cellI)
            {
                bf[ibc[cellI]]=wss[cellI];
            }

        }
        else if((isA<immersedBoundaryFvPatchVectorField>(U.boundaryField()[patchI])))
        {
           FatalErrorIn("update_wall_shear.H")
              << "U should use immersedBoundaryVelocityWallFunctionFvPatchVectorField" << nl
              << abort(FatalError);

        }

        scalarField wallShearStressMag(mag(wallShearStress));

        fileName path;
        if (Pstream::parRun())
        {
            path = mesh.time().path()/".."/"postProcessing"/"VTK";
        }
        else
        {
            path = mesh.time().path()/"postProcessing"/"VTK";
        }

        if (runTime.outputTime())
        {
            reduce(wallShearStress, sumOp<vectorList>());
                       
            if(Pstream::myProcNo()==Pstream::masterNo())
            {
                const triSurface& ts = ibPatch.ibMesh();

                // Make a face list for writing
                faceList f(ts.size());
                forAll (ts, faceI)
                {
                    f[faceI] = ts[faceI].triFaceFace();
                }

                // Write immersed boundary data as a vtk file
                autoPtr<ibSurfaceWriter<vector> > writerPtr1 =
                ibSurfaceWriter<vector>::New("vtk");

                fileName path2 = path/"wallShearStress";
                mkDir(path);
                writerPtr1->write
                (
                    path2,
                    ibPatch.name()+"_"+mesh.time().timeName(),
                    ts.points(),
                    f,
                    "wallShearStress",
                    wallShearStress
                );

            }
        }
    }
}
