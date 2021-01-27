//
// createIbMask.H
// ~~~~~~~~~~~~

    Info<< "Create immersed boundary cell mask" << endl;

    volScalarField cellIbMask
    (
        IOobject
        (
            "cellIbMask",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("one", dimless, 1)
    );

    volScalarField cellIbMaskExt
    (
        IOobject
        (
            "cellIbMaskExt",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("one", dimless, 1)
    );

    volScalarField ibCellMask
    (
        IOobject
        (
            "ibCellMask",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0)
    );

    volScalarField samplingCellMask
    (
        IOobject
        (
            "samplingCellMask",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0)
    );

    Info<< "Create immersed boundary face mask" << endl;
    surfaceScalarField faceIbMask
    (
        IOobject
        (
            "faceIbMask",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("one", dimless, 1)
    );

    IOdictionary ibmDict
        (
            IOobject
            (
                "ibmDict",
                mesh.time().system(),
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        );


    // read constants from ibmDict
    vector gD(ibmDict.lookup("gDirection"));

	//if correct yplus in yPlusCorrection.H
	bool correctYPlus(ibmDict.lookupOrDefault("correctYPlus",false));	
	
	bool writeStlEachTimeStep(ibmDict.lookupOrDefault("writeStlEachTimeStep",false));


    forAll (mesh.boundary(), patchI)
    {
        if (isA<immersedBoundaryFvPatch>(mesh.boundary()[patchI]))
        {
            Info<< "Found immersed boundary patch " << patchI
                << " named " << mesh.boundary()[patchI].name()
                << endl;

            const immersedBoundaryFvPatch& ibPatch =
                refCast<const immersedBoundaryFvPatch>
                (
                    mesh.boundary()[patchI]
                );

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

    // Evaluate boundary conditions for IB masks
     //cellIbMask.boundaryField()[patchI].evaluate();
     //cellIbMaskExt.boundaryField()[patchI].evaluate();
    // Evaluate the coupled patchField, to replace evaluateCoupled()
    // Y.C. Xu 7/5/2017
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




 

