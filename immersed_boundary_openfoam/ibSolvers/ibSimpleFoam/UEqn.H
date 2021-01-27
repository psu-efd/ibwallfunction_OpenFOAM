    // Solve the Momentum equation

    MRF.correctBoundaryVelocity(U);

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + MRF.DDt(U)
      + turbulence->divDevReff(U)
      ==
        fvOptions(U)
    );

    UEqn().relax();

    fvOptions.constrain(UEqn());

    tmp<fvVectorMatrix> UEqn1
    (
      UEqn() == -fvc::grad(p)
    );

    //UEqn().boundaryManipulate(U.boundaryField());
    UEqn1().boundaryManipulate(U.boundaryField()); //In official version, preconditioner needs to be turned off
    //solve(UEqn() == -fvc::grad(p));
    //#include "manipulateU.H"
    solve(UEqn1());

    fvOptions.correct(U);

