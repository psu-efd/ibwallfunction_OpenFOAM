forAll (mesh.boundary(), patchI)
    {
        if (isA<immersedBoundaryFvPatch>(mesh.boundary()[patchI]))
        {

            const immersedBoundaryFvPatch& ibPatch_ =
                refCast<const immersedBoundaryFvPatch>
                (
                    mesh.boundary()[patchI]
                );
        //const immersedBoundaryFvPatch& ibPatch= ibPatch_;
        //#include "test.H"       

        fvMatrix<vector>& eqn = UEqn1();
 
        Field<vector> polyPsi(eqn.psi(), ibPatch_.ibCells());

        eqn.setValues(ibPatch_.ibCells(), polyPsi);


        Field<vector> deadCellsPsi
        (
            eqn.psi(),
            ibPatch_.deadCells()
        ); 

        eqn.setValues(ibPatch_.deadCells(), deadCellsPsi);
        Info<<"Finished manipulate UEqn1"<<endl;
        }     
    }

