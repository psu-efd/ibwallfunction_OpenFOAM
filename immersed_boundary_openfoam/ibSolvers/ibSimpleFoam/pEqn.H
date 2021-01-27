{
    volScalarField rAU(1.0/UEqn().A());
    volVectorField HbyA("HbyA", U);
    HbyA = rAU*UEqn().H();

    // Immersed boundary update
    HbyA.correctBoundaryConditions();

    surfaceScalarField phiHbyA("phiHbyA", faceIbMask*fvc::interpolate(HbyA) & mesh.Sf());
    phiHbyA *= faceIbMask;

    // Adjust immersed boundary fluxes
    immersedBoundaryAdjustPhi(phiHbyA, U);
    MRF.makeRelative(phiHbyA);
    adjustPhi(phiHbyA, U, p);

    tmp<volScalarField> rAtU(rAU);

    if (simple.consistent())
    {
        rAtU = 1.0/(1.0/rAU - UEqn().H1());
        phiHbyA +=
            faceIbMask*fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf();
        HbyA -= (rAU - rAtU())*fvc::grad(p);
    }

    UEqn.clear();

    // Non-orthogonal pressure corrector loop
    while (simple.correctNonOrthogonal())
    {   
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
        );
        //pRefCell=7597;
        pEqn.setReference(pRefCell, pRefValue);
		if(pEqnBoundaryManipulate)
		{
			pEqn.boundaryManipulate(p.boundaryField());
		}
        #include "manipulateP.H"
        pEqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - pEqn.flux();
        }
    }

#   include "immersedBoundaryContinuityErrs.H"

    // Explicitly relax pressure for momentum corrector
    p.relax();

    // Momentum corrector
    U = HbyA - rAtU()*fvc::grad(p);
    U.correctBoundaryConditions();
    fvOptions.correct(U);
}
