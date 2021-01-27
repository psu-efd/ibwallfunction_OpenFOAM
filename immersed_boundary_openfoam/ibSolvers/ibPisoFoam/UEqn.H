// Solve the Momentum equation


fvVectorMatrix UEqn
(
    fvm::ddt(U) 
    + fvm::div(phi, U)
    + turbulence->divDevReff(U)
);

UEqn.relax();

if (piso.momentumPredictor())
{
    tmp<fvVectorMatrix> UEqn1
    (
      UEqn == -fvc::grad(p)
    );
//#include "checkPsi.H"
    //UEqn().boundaryManipulate(U.boundaryField());
    UEqn1().boundaryManipulate(U.boundaryField()); //In official version, preconditioner needs to be turned off
//#include "checkPsi.H"
    solve(UEqn1());
    //solve(UEqn == -fvc::grad(p));
}
