    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

#   include "createPhi.H"
    //create a body force field for the immersed boundary which is only
    //relevant for fluid cells
/*    volVectorField bodyForce
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
*/

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, simple.dict(), pRefCell, pRefValue);

    mesh.setFluxRequired(p.name());

    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );

    bool pEqnBoundaryManipulate(ibmDict.lookupOrDefault("pEqnBoundaryManipulate",false));
	
	
