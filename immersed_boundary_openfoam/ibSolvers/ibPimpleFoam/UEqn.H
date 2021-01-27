/*    tmp<volTensorField> tgradU = fvc::grad(U);

    volTensorField gradU
    (
        IOobject
        (
            "gradU",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedTensor("0",dimensionSet(0, 0, -1, 0, 0, 0, 0),tensor::zero)
    );

    gradU=tgradU;

    volScalarField procCell
    (
        IOobject
        (
            "procCell",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0",dimensionSet(0, 0, 0, 0, 0, 0, 0),0)
    );

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells =
            mesh.boundary()[patchi].faceCells();
 
    	if(isA<processorPolyPatch>(mesh.boundaryMesh()[patchi]))
        {
        forAll(mesh.boundary()[patchi], facei)
        {
            procCell[pFaceCells[facei]]=1;
        }
        }
    }

    volVectorField UU
    (
        IOobject
        (
            "UU",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("0",dimensionSet(0, 1, -1, 0, 0, 0, 0),vector::zero)
    );
    U.correctBoundaryConditions();
    UU=U;
 
    volScalarField nutt   
    (
        "nutt",
        turbulence->nut()
    );
    tmp<volTensorField> tgradU = fvc::grad(U);

    volScalarField G
    (
        "G",
        nutt*(tgradU() && dev(twoSymm(tgradU())))
    );

    volScalarField GG
    (
        IOobject
        (
            "GG",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero",G.dimensions(),0)
    );

    GG=G;
*/

// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

//calculate the body force
//volVectorField U_bak(U);

 

//bodyForce = fvc::ddt(U);
/*
surfaceScalarField nuEff
(
    "nuEff",
    fvc::interpolate(turbulence->nuEff())
);
 
bodyForce +=
    (
        fvc::div(phi, U)
        - fvc::laplacian(nuEff, U)
        + fvc::grad(p)
    );*/
//finished calculating the body force

tmp<fvVectorMatrix> UEqn
(
    fvm::ddt(U)
    + fvm::div(phi, U)
    + turbulence->divDevReff(U)
 ==
    fvOptions(U)
);

UEqn().relax();

fvOptions.constrain(UEqn());

if (pimple.momentumPredictor())
{
    tmp<fvVectorMatrix> UEqn1
    (
      UEqn() == -fvc::grad(p)
    );
    UEqn1().boundaryManipulate(U.boundaryField());

    //#include "manipulateU.H"

    solve(UEqn1());
    fvOptions.correct(U);
}
