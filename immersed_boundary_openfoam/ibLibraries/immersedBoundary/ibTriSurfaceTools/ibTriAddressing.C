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

#include "ibTriSurfaceTools.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
 // address the trisurface points added by Y. Xu   return local points
void Foam::ibTriSurfaceTools::makeTriPointsToHitPoints() const
{
 

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (!triPointsToHitPointsAddrPtr_.empty() || !triPointsToHitPointsWeightsPtr_.empty())
    {
        FatalErrorIn("ibTriSurfaceTools::makeTriPointsToHitPoints() const")
            << "tri addressing already exist"
            << abort(FatalError);
    }

    const triSurface& triPatch = triSurface_;
    const vectorField& triCentres = triPatch.faceCentres();
    const vectorField& triPoints = triPatch.localPoints();
    // point-face addressing of trisurface
    const labelListList& triPointFaces = triPatch.pointFaces();

    const labelListList& triFaceFaces = triPatch.faceFaces();
//    triPointsToHitPointsAddrPtr_ = new labelListList(triPointFaces.size());
    labelListList& addr = triPointsToHitPointsAddrPtr_;
    
    addr=triPointFaces;

//    triPointsToHitPointsWeightsPtr_ = new scalarListList(triPointFaces.size());
    scalarListList& w = triPointsToHitPointsWeightsPtr_;
    
    w.setSize(addr.size());
 
    forAll (triPointFaces,ptI)
    {

		labelHashSet curAddrSet;
		
		labelList triPF = triPointFaces[ptI];
		forAll(triPF, triI)
		{

			label adjacentTriI=triPF[triI];

			if(!curAddrSet.found(adjacentTriI))
			{
				curAddrSet.insert(adjacentTriI);
			}

			forAll(triFaceFaces[adjacentTriI], triII)
			{
				if(!curAddrSet.found(triFaceFaces[adjacentTriI][triII]))
				{
					curAddrSet.insert(triFaceFaces[adjacentTriI][triII]);
				}		
			}
		}
		addr[ptI] = curAddrSet.toc();

        w[ptI].setSize(addr[ptI].size());

        labelList& curAddr = addr[ptI];
        scalarList& curW = w[ptI];
        const vector& curTriPoint = triPoints[ptI];

        scalar sumW = 0.0;


        forAll (curAddr, ibI)
        {

			vector Length=curTriPoint - triCentres[curAddr[ibI]];
			//Length.z()=0;

            curW[ibI] = 1.0/(mag(Length)+SMALL);

            sumW += curW[ibI];
        }

            // Divide weights by sum distance
        forAll (curW, ibI)
        {
            curW[ibI] /= sumW;

        }

    }

//   Pout << "addr[101] = " <<  triPointsToHitPointsWeightsPtr_[101] << endl;
//   Pout << "w[101] = " <<  w[101] << endl;


}

const Foam::labelListList&
Foam::ibTriSurfaceTools::triPointsToHitPointsAddr() const
{
    if (triPointsToHitPointsAddrPtr_.empty())
    {
        makeTriPointsToHitPoints();
    }

    return triPointsToHitPointsAddrPtr_;
}


const Foam::scalarListList&
Foam::ibTriSurfaceTools::triPointsToHitPointsWeights() const
{
    if (triPointsToHitPointsWeightsPtr_.empty())
    {
        makeTriPointsToHitPoints();
    }

    return triPointsToHitPointsWeightsPtr_;
}


template<class Type>
Foam::tmp<Foam::Field<Type> > 
Foam::ibTriSurfaceTools::toTriPoints
(
    const Field<Type>& triFaceCentresValues
) const
{
    //Pout << "triFaceCentresValues.size()= " << triFaceCentresValues.size()<< endl;
    if (triFaceCentresValues.size() != triSurface_.faceCentres().size())
    {
        FatalErrorIn
        (
            "template<class Type>\n"
            "Foam::tmp<Foam::Field<Type> >\n"
            "ibTriSurfaceTools::toTriPoints\n"
            "(\n"
            "    const Field<Type>& triFaceCentresValues\n"
            ") const"
        )   << "Field size does not correspond to size of face centres of "
            << "triangulated surface for patch "<< nl
            << "Field size = " << triFaceCentresValues.size()
            << "Triangulated surface face size = " << triSurface_.faceCentres().size()
            << abort(FatalError);
    }

//    Pout << "ibTriSurfaceTools::toTriPoints() " << endl;

    const labelListList& ctfAddr = triPointsToHitPointsAddr();
    const scalarListList& ctfWeights = triPointsToHitPointsWeights();

   //Pout << "triPointsToHitPointsAddr = " << ctfAddr.size() << endl;
   //Pout << "triPointsToHitPointsWeights = " << ctfWeights.size() << endl;

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
Foam::ibTriSurfaceTools::toTriPoints
(
    const tmp<const Field<Type>& > ttriFaceCentresValues
) const
{
    tmp<Field<Type> >  tint = toTriPoints(ttriFaceCentresValues());
    ttriFaceCentresValues.clear();
    return tint;
}

// added by Y. Xu   Interpolate gradient of center values at each cell face on the triSurface 
template<class Type>
Foam::tmp<Foam::vectorField > 
Foam::ibTriSurfaceTools::gradTriFaces
(
    const Field<Type>& triFaceCentresValues
) const
{
//    Pout << "triFaceCentresValues = " << triFaceCentresValues.size() << endl;
    if (triFaceCentresValues.size() != triSurface_.faceCentres().size())
    {
        FatalErrorIn
        (
            "template<class Type>\n"
            "Foam::tmp<Foam::Field<Type> >\n"
            "ibTriSurfaceTools::gradTriFaces\n"
            "(\n"
            "    const Field<Type>& triFaceCentresValues\n"
            ") const"
        )   << "Field size does not correspond to size of face centres of "
            << "triangulated surface for patch " << nl
            << "Field size = " << triFaceCentresValues.size()
            << "Triangulated surface face size = " << triSurface_.faceCentres().size()
            << abort(FatalError);
    }
    Field<Type> triPointValues = toTriPoints(triFaceCentresValues);

    const labelListList& ffAddress(triSurface_.faceFaces()); 
    const labelListList& feAddress(triSurface_.faceEdges()); 
    const edgeList& edgeAddress(triSurface_.edges());
	const pointField pts(triSurface_.localPoints());
    //const labelListList& efAddress(triSurface_.edgeFaces());
    const vectorField& triCentres = triSurface_.faceCentres();
    labelListList addr_c(ffAddress);
    scalarListList w_c(ffAddress.size());

    tmp<vectorField> tIbPsi
    (
        new vectorField(addr_c.size(), pTraits<vector>::zero)
    );
    vectorField& ibPsi = tIbPsi();

    forAll(feAddress,triI)
    {
		const labelList& feA=feAddress[triI];
		forAll(feA, fI)
		{
			edge edgeA=edgeAddress[feA[fI]];
			point ptC=triCentres[triI];
			point ptA=pts[edgeA[0]];
			point ptB=pts[edgeA[1]];
			vector AB=ptB-ptA;
			vector AC=ptC-ptA;
			scalar c_=AB&AC/(AB&AB);
			vector edgeNorm=c_*AB-AC;
			edgeNorm =edgeNorm/mag(edgeNorm);

			//const labelList& efA=efAddress[fI];
			scalar edgePhi = 0.5*(triPointValues[edgeA[0]]+triPointValues[edgeA[1]]);
			ibPsi[triI] +=edgePhi*edgeNorm*mag(AB);
		}
		ibPsi[triI] /=triSurface_[triI].mag(pts);
	}

    return tIbPsi;
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::ibTriSurfaceTools::gradTriFaces
(
    const tmp<const Field<Type>& > ttriFaceCentresValues
) const
{
    tmp<Field<Type> >  tint = gradTriFaces(ttriFaceCentresValues());
    ttriFaceCentresValues.clear();
    return tint;
}

// ************************************************************************* //
