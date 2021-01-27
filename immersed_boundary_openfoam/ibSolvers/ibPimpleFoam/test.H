 
     //const   volScalarField& psiI = const_cast<volScalarField&>
      //  (db().lookupObject<volScalarField>("p"));
const scalarField& psi = p.internalField();
    const labelListList& procCells = ibPatch.ibProcCells(); 
 

    tmp<FieldField<Field, scalar> > tprocPsi
    (
        new FieldField<Field, scalar>(Pstream::nProcs())
    );
    FieldField<Field, scalar>& procPsi = tprocPsi();

    //const labelListList& procCells = ibProcCells();
    // This requires a rewrite useng mapDistribute
    // HJ, 11/Aug/2016
    forAll (procPsi, procI)
    {
        procPsi.set
        (
            procI,
            new Field<scalar>
            (
                procCells[procI].size(),
                pTraits<scalar>::zero
            )
        );
    }
    /* 
   
    typedef List<Field<scalar> > FieldTypeList;
    if (Pstream::parRun())
    {     
        FieldTypeList procPsiList(Pstream::nProcs());  

        Field<scalar> curPsi(psi, procCells[Pstream::myProcNo()]);
        
        procPsiList[Pstream::myProcNo()] = curPsi;

        Pstream::gatherList(procPsiList); 
        Pstream::scatterList(procPsiList);

        forAll (procPsi, procI)
        {
            procPsi[procI] = procPsiList[procI];
        } 
    }
*/


    if (Pstream::parRun())
    {

        // Send
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                // Do not send empty lists
                if (!procCells[procI].empty())
                {
                    Field<scalar> curPsi(psi, procCells[procI]);

                    // Parallel data exchange
                    //OPstream toProc(Pstream::blocking,procI);
                    OPstream toProc
                    (
                        Pstream::blocking,
                        procI,
						curPsi.size()*sizeof(scalar)
                    );

                    toProc << curPsi;
                }
            }
        }


        // Receive
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
 
            if (procI != Pstream::myProcNo())
            {

                // Do not receive empty lists
                if (!procPsi[procI].empty())
                {
                    // Parallel data exchange
                    IPstream fromProc
                    (
                        Pstream::blocking,
                        procI,
						procPsi[procI].size()*sizeof(scalar)
                    );

                    fromProc >> procPsi[procI];
                }
            }
        }
    }

Pout<<procPsi[3][16]<<endl;


forAll(procPsi,id)
{
//Pout<<id<<" "<<procCells[id].size()*sizeof(scalar)<<endl;
	scalar nnn=0;
	forAll(procPsi[id],iid)
	{
		if(procPsi[id][iid]<SMALL)
		{
			nnn++;
		}
	}
    //Pout<<nnn<<endl;
}
