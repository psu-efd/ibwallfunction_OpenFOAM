/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "immersedBoundaryFvPatch.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::FieldField<Foam::Field, Type> >
Foam::immersedBoundaryFvPatch::sendAndReceive
(
    const Field<Type>& psi
) const
{

    tmp<FieldField<Field, Type> > tprocPsi
    (
        new FieldField<Field, Type>(Pstream::nProcs())
    );
    FieldField<Field, Type>& procPsi = tprocPsi();

    const labelListList& procCells = ibProcCells();
    // This requires a rewrite useng mapDistribute
    // HJ, 11/Aug/2016
 
    forAll (procPsi, procI)
    {
        procPsi.set
        (
            procI,
            new Field<Type>
            (
                ibProcCentres()[procI].size(),
                pTraits<Type>::zero
            )
        );
    }

/*
    typedef List<Field<Type> > FieldTypeList;
    if (Pstream::parRun())
    {     
        FieldTypeList procPsiList(Pstream::nProcs());  

        Field<Type> curPsi(psi, procCells[Pstream::myProcNo()]);
        
        procPsiList[Pstream::myProcNo()] = curPsi;

        Pstream::gatherList(procPsiList); 
        Pstream::scatterList(procPsiList);

        forAll (procPsi, procI)
        {
            procPsi[procI] = procPsiList[procI];
        } 
    }


*/
	forAll(procCells,i)
	{
		//Pout<<Pstream::myProcNo()<<" "<<i<<" "<<procCells[i].size()<<endl;
	}
    if (Pstream::parRun())
    {
        PstreamBuffers pBuffers(Pstream::nonBlocking);
        // Send
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                // Do not send empty lists
                if (!procCells[procI].empty())
                {
                    Field<Type> curPsi(psi, procCells[procI]);

                    // Parallel data exchange
//                    OPstream toProc(Pstream::blocking,procI);
                    OPstream toProc
                    (
                        Pstream::blocking,
                        procI,
						curPsi.size()*sizeof(Type)

                    );


                    toProc << curPsi;
                }
            }
        }

        //pBuffers.finishedSends();   // no-op for blocking

        // Receive
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
 
            if (procI != Pstream::myProcNo())
            {

                // Do not receive empty lists
                if (!procPsi[procI].empty())
                {

                    // Parallel data exchange
//                    IPstream fromProc(Pstream::blocking,procI);
                    IPstream fromProc
                    (
                        Pstream::blocking,
                        procI,
						procPsi[procI].size()*sizeof(Type)
                    );

                    fromProc >> procPsi[procI];

                }
            }
        }
    }

 //Pout<<"see if it works"<<endl;
    return tprocPsi;
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::immersedBoundaryFvPatch::toIbPoints
(
    const Field<Type>& triValues
) const
{
    if (triValues.size() != ibMesh().size())
    {
        FatalErrorIn
        (
            "template<class Type>\n"
            "Foam::tmp<Foam::Field<Type> >\n"
            "immersedBoundaryFvPatch::toIbPoints\n"
            "(\n"
            "    const Field<Type>& triValues\n"
            ") const"
        )   << "Field size does not correspond to size of immersed boundary "
            << "triangulated surface for patch " << name() << nl
            << "Field size = " << triValues.size()
            << " surface size = " << ibMesh().size()
            << abort(FatalError);
    }

    const labelList& ibc = ibCells();

    tmp<Field<Type> > tIbPsi
    (
        new Field<Type>(ibc.size(), pTraits<Type>::zero)
    );
    Field<Type>& ibPsi = tIbPsi();

    const labelList& hf = hitFaces();

    // Assuming triSurface data is on triangles
    forAll (ibPsi, cellI)
    {
        ibPsi[cellI] = triValues[hf[cellI]];
    }

//     const vectorField& p = ibPoints();
//     const List<labelledTri>& faces = ibMesh();
//     const vectorField& triPoints = ibMesh().points();

//     // Assuming triSurface data is on vertices
//     forAll (ibPsi, cellI)
//     {
//         const labelledTri& tri = faces[hf[cellI]];
//         triPointRef triPt = faces[hf[cellI]].tri(triPoints);

//         ibPsi[cellI] =
//             triValues[tri[0]]*triPt.Ni(0, p[cellI])
//           + triValues[tri[1]]*triPt.Ni(1, p[cellI])
//           + triValues[tri[2]]*triPt.Ni(2, p[cellI]);
//     }

    return tIbPsi;
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::immersedBoundaryFvPatch::toIbPoints
(
    const tmp<Field<Type> >& ttriValues
) const
{
    tmp<Field<Type> > tint = toIbPoints(ttriValues());
    ttriValues.clear();
    return tint;

}

template<class Type>
Foam::tmp<Foam::Field<Type> > 
Foam::immersedBoundaryFvPatch::toTriPoints
(
    const Field<Type>& triFaceCentresValues
) const
{
    //Pout << "triFaceCentresValues.size()= " << triFaceCentresValues.size()<< endl;
    if (triFaceCentresValues.size() != ibMesh().faceCentres().size())
    {
        FatalErrorIn
        (
            "template<class Type>\n"
            "Foam::tmp<Foam::Field<Type> >\n"
            "immersedBoundaryFvPatch::toTriPoints\n"
            "(\n"
            "    const Field<Type>& triFaceCentresValues\n"
            ") const"
        )   << "Field size does not correspond to size of face centres of "
            << "triangulated surface for patch " << name() << nl
            << "Field size = " << triFaceCentresValues.size()
            << " Triangulated surface face size = " << ibMesh().faceCentres().size()
            << abort(FatalError);
    }

//    Pout << "immersedBoundaryFvPatch::toTriPoints() " << endl;
 
	const labelListList& ctfAddr = ibMesh().pointFaces();
	const scalarListList& ctfWeights = triFacesToTriPointsWeights();

	tmp<Field<Type> > tIbPsi
	(
		new Field<Type>(ctfAddr.size(), pTraits<Type>::zero)
	);
	Field<Type>& ibPsi = tIbPsi();

	// Do interpolation
	forAll (ctfAddr, triI)
	{
		const labelList& curAddr = ctfAddr[triI];
		const scalarList& curWeights = ctfWeights[triI];

		forAll (curAddr, cellI)
		{
			ibPsi[triI] += curWeights[cellI]*triFaceCentresValues[curAddr[cellI]];
		}
	}
 

    return tIbPsi;
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::immersedBoundaryFvPatch::toTriPoints
(
    const tmp<const Field<Type>& > ttriFaceCentresValues
) const
{
    tmp<Field<Type> >  tint = toTriPoints(ttriFaceCentresValues());
    ttriFaceCentresValues.clear();
    return tint;
}

template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::immersedBoundaryFvPatch::toTriFaces
(
    const Field<Type>& ibValues
) const
{
    if (ibValues.size() != ibCells().size())
    {
        FatalErrorIn
        (
            "template<class Type>\n"
            "Foam::tmp<Foam::Field<Type> >\n"
            "immersedBoundaryFvPatch::toTriFaces\n"
            "(\n"
            "    const Field<Type>& ibValues\n"
            ") const"
        )   << "Field size does not correspond to size of IB points "
            << "triangulated surface for patch " << name() << nl
            << "Field size = " << ibValues.size()
            << " IB points size = " << ibCells().size()
            << abort(FatalError);
    }

    const labelListList& ctfAddr = cellsToTriAddr();
    const scalarListList& ctfWeights = cellsToTriWeights();

    tmp<Field<Type> > tIbPsi
    (
        new Field<Type>(ctfAddr.size(), pTraits<Type>::zero)
    );
    Field<Type>& ibPsi = tIbPsi();

    // Do interpolation
    forAll (ctfAddr, triI)
    {
        const labelList& curAddr = ctfAddr[triI];
        const scalarList& curWeights = ctfWeights[triI];

        forAll (curAddr, cellI)
        {
            ibPsi[triI] += curWeights[cellI]*ibValues[curAddr[cellI]];
        }
    }

    return tIbPsi;
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::immersedBoundaryFvPatch::toTriFaces
(
    const tmp<Field<Type> >& tibValues
) const
{
    tmp<Field<Type> > tint = toTriFaces(tibValues());
    tibValues.clear();
    return tint;

}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::immersedBoundaryFvPatch::toSamplingPoints
(
    const Field<Type>& cellValues
) const
{
    if (cellValues.size() != mesh_.nCells())
    {
        FatalErrorIn
        (
            "template<class Type>\n"
            "Foam::tmp<Foam::Field<Type> >\n"
            "immersedBoundaryFvPatch::toSamplingPoints\n"
            "(\n"
            "    const Field<Type>& cellValues\n"
            ") const"
        )   << "Field size does not correspond to cell centres "
            << "for patch " << name() << nl
            << "Field size = " << cellValues.size()
            << " nCells = " << mesh_.nCells()
            << abort(FatalError);
    }

    // Get addressing
    const labelList& ibc = ibCells();
    const labelListList& ibcc = ibCellCells();
    const List<List<labelPair> >& ibcProcC = ibCellProcCells();
    //const labelListList& procCells = ibProcCells();

    // Get weights
    const scalarListList& cellWeights = ibSamplingWeights();
    const scalarListList& cellProcWeights = ibSamplingProcWeights();
 
    tmp<Field<Type> > tIbPsi
    (
        new Field<Type>(ibc.size(), pTraits<Type>::zero)
    );
    Field<Type>& ibPsi = tIbPsi();

    // Do interpolation, local cell data
 
    forAll (ibc, cellI)
    {
        const labelList& curAddr = ibcc[cellI];
        const scalarList& curWeights = cellWeights[cellI];

        forAll (curAddr, ccI)
        {
            ibPsi[cellI] += curWeights[ccI]*cellValues[curAddr[ccI]];
        }
    }

    // Parallel communication for psi
    FieldField<Field, Type> procCellValues = sendAndReceive(cellValues);

    // Do interpolation, cell data from other processors
    forAll (ibc, cellI)
    {
        const List<labelPair>& curProcCells = ibcProcC[cellI];
        const scalarList& curProcWeights = cellProcWeights[cellI];

        forAll (curProcCells, cpcI)
        {
            ibPsi[cellI] +=
                curProcWeights[cpcI]*
                procCellValues
                [
                    curProcCells[cpcI].first()
                ]
                [
                    curProcCells[cpcI].second()
                ];
        }
    }

    return tIbPsi;
}


template<class Type>
const Foam::tmp<Foam::Field<Type> >
Foam::immersedBoundaryFvPatch::renumberField
(
    const   Field<Type>& f
) const
{
    const dynamicLabelList& triFInM = this->triFacesInMesh();

    tmp<Field<Type> > trf(new Field<Type>(triFInM.size()));
    Field<Type>& rf = trf();

    forAll(triFInM, faceI)
    {
        rf[faceI] = f[triFInM[faceI]];
    }

    return trf;
}

// added by Y. Xu   Interpolate gradient of center values at each cell face on the triSurface 
template<class Type>
Foam::tmp<Foam::vectorField > 
Foam::immersedBoundaryFvPatch::gradTriFaces
(
    const Field<Type>& triFaceCentresValues
) const
{
 
	//Field<Type> triPointValues = toTriPoints(triFaceCentresValues);

    const labelListList& ffAddress(ibMesh().faceFaces()); 
    const labelListList& feAddress(ibMesh().faceEdges()); 
    const edgeList& edgeAddress(ibMesh().edges());
	const pointField pts(ibMesh().localPoints());
    //const labelListList& efAddress(ibMesh().edgeFaces());
    //const vectorField& triCentres = ibMesh().faceCentres();
	const labelList& edgeOwner=ibMesh().edgeOwner();

	const scalarListList& triFacesInMeshType=this->triFacesInMeshType();
	const vectorListList& triFacesInMeshTypeVec=this->triFacesInMeshTypeVec();
	const labelList& triFacesInMesh=this->triFacesInMesh();

	const vectorField& triSf=this->triSf();

	const tmp<scalarField> triSfMag = mag(triSf);
	 
	tmp<Field<Type> > tedgePhi(new Field<Type>(this->triFacesToTriEdges(triFaceCentresValues)));
	Field<Type>& edgePhi = tedgePhi();	


 
    labelListList addr_c(ffAddress);
    scalarListList w_c(ffAddress.size());
    IOdictionary ibmDict
    (
        IOobject
        (
            "ibmDict",
            mesh_.time().system(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

	vector gD(ibmDict.lookup("gDirection"));

    tmp<vectorField> tIbPsi
    (
        new vectorField(addr_c.size(), pTraits<vector>::zero)
    );
    vectorField& ibPsi = tIbPsi();

    forAll(triFacesInMesh,faceI)
    {
		const label& triI = triFacesInMesh[faceI];
		const labelList& feA=feAddress[triI];
		forAll(feA, fI)
		{
			const edge& edgeA=edgeAddress[feA[fI]];
			vector vec=edgeA.vec(pts);
			vector edgeNorm=vec;
			//edgePhi[feA[fI]] = 0.5*(triPointValues[edgeA[0]]+triPointValues[edgeA[1]]);
			if(edgeOwner[feA[fI]]==triI)
			{
				ibPsi[triI] +=edgePhi[feA[fI]] *(edgeNorm^(gD/mag(gD)));				
			}
			else
			{
				ibPsi[triI] -=edgePhi[feA[fI]] *(edgeNorm^(gD/mag(gD)));
			}
		}
		ibPsi[triI] /=triSfMag()[triI];

		const scalarList& curveType = triFacesInMeshType[triI];
		forAll(curveType, triII)
		{
			if(curveType[triII]==2 )// patch zeroGradient
			{
				vector Vec = triFacesInMeshTypeVec[triI][triII];
				Vec = Vec/mag(Vec);
				ibPsi[triI] = ibPsi[triI] - (ibPsi[triI]&Vec)*Vec;
			}
		}
	}

    //Pout << "interploated values = " << ibPsi << endl;

    return tIbPsi;
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::immersedBoundaryFvPatch::gradTriFaces
(
    const tmp<const Field<Type>& > ttriFaceCentresValues
) const
{
    tmp<Field<Type> >  tint = gradTriFaces(ttriFaceCentresValues());
    ttriFaceCentresValues.clear();
    return tint;
}

// added by Y. Xu   Interpolate gradient of center values at each cell face on the triSurface 
template<class Type>
Foam::tmp<Foam::Field<Type> > 
Foam::immersedBoundaryFvPatch::triFacesToTriEdges
(
    const Field<Type>& triFaceCentresValues
) const
{
 
    tmp<Field<Type> > tresult
    (
        new Field<Type>(ibMesh().nEdges(), pTraits<Type>::zero)
    );

    Field<Type>& result = tresult();

    const edgeList& edges = ibMesh().edges();
    const labelListList& edgeFaces = ibMesh().edgeFaces();

    const scalarList& weights = triFacesToTriEdgesWeights();

    for (label edgei = 0; edgei < ibMesh().nInternalEdges(); edgei++)
    {
        result[edgei] =
            weights[edgei]*triFaceCentresValues[edgeFaces[edgei][0]]
          + (1.0 - weights[edgei])*triFaceCentresValues[edgeFaces[edgei][1]];
    }

    for (label edgei = ibMesh().nInternalEdges(); edgei < edges.size(); edgei++)
    {
        result[edgei] = triFaceCentresValues[edgeFaces[edgei][0]];
    }

    return tresult;
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::immersedBoundaryFvPatch::triFacesToTriEdges
(
    const tmp<const Field<Type>& > ttriFaceCentresValues
) const
{
    tmp<Field<Type> >  tint = triFacesToTriEdges(ttriFaceCentresValues());
    ttriFaceCentresValues.clear();
    return tint;
}

// added by Y. Xu   Interpolate gradient of center values at each cell face on the triSurface 
template<class Type>
Foam::tmp<Foam::Field<Type> > 
Foam::immersedBoundaryFvPatch::triPointsToTriFaces
(
    const Field<Type>& triFaceCentresValues
) const
{
    if (triFaceCentresValues.size() != ibMesh().nPoints())
    {
        FatalErrorIn
        (
            "immersedBoundaryFvPatch::"
            "triPointsToTriFaces(const scalarField& triFaceCentresValues)"
        )   << "given field does not correspond to patch. Patch size: "
            << ibMesh().nPoints() << " field size: " << triFaceCentresValues.size()
            << abort(FatalError);
    }

    tmp<Field<Type> > tresult
    (
        new Field<Type>
        (
            ibMesh().size(),
            pTraits<Type>::zero
        )
    );

    Field<Type>& result = tresult();

    const List<labelledTri>& localFaces = ibMesh().localFaces();

    forAll(result, facei)
    {
        const labelList& curPoints = localFaces[facei];

        forAll(curPoints, pointi)
        {
            result[facei] += triFaceCentresValues[curPoints[pointi]];
        }

        result[facei] /= curPoints.size();
    }

    return tresult;

}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::immersedBoundaryFvPatch::triPointsToTriFaces
(
    const tmp<const Field<Type>& > ttriFaceCentresValues
) const
{
    tmp<Field<Type> >  tint = triPointsToTriFaces(ttriFaceCentresValues());
    ttriFaceCentresValues.clear();
    return tint;
}

// ************************************************************************ //
