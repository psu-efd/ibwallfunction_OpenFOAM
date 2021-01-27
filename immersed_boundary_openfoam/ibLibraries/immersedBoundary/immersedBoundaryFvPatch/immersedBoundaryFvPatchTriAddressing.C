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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void Foam::immersedBoundaryFvPatch::makeTriAddressing() const
{
	// Updates has include triFace with at least one pt in coupled patch in parallel, Xu Dec2017

	// it takes too much time, need to be optimized
 	const double Oldtime0=mesh_.time().elapsedCpuTime();

    if (debug)
    {
        InfoIn("void immersedBoundaryFvPatch::makeTriAddressing() const")
            << "creating tri addressing for immersed boundary " << name()
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (cellsToTriAddrPtr_ || cellsToTriWeightsPtr_)
    {
        FatalErrorIn("immersedBoundaryFvPatch::makeTriAddressing() const")
            << "tri addressing already exist"
            << "for immersed boundary" << name()
            << abort(FatalError);
    }

    // Get reference to tri patch and hit faces
    const triSurface& triPatch = ibPolyPatch_.ibMesh();
    const vectorField& triCentres = triPatch.faceCentres();

    const labelList& hf = hitFaces();
    const vectorField& ibp = ibPoints();

    // Create a markup field and mark all tris containing an ib point with its
    // index
 
    labelListList hitTris(triPatch.size());

    forAll (hf, hfI)
    {

		//hitTriSetList[hf[hfI]].insert(hfI);
        //To fix: need to consider the case that some ghost cells might don't
        //     have any hit at all (hf[hfI]=-1)
        hitTris[hf[hfI]].append(hfI);
    }
 
    // Allocate storage
    cellsToTriAddrPtr_ = new labelListList(triPatch.size());
    labelListList& addr = *cellsToTriAddrPtr_;

    cellsToTriWeightsPtr_ = new scalarListList(triPatch.size());
    scalarListList& w = *cellsToTriWeightsPtr_;

    // Algorithm:
    // For each tri face, check if it contains an IB point
    // - if so, set the addressing to the index of IB point and weight to 1
    // - if not, search the neighbouring faces of the visited faces until
    //   at least 3 IB points are found, or the neighbourhood is exhausted.
    //   When a sufficient number of points is found, calculate the weights
    //   using inverse distance weighting

    // Get addressing from the triangular patch
    //const labelListList& pf = triPatch.pointFaces();
    const labelListList& faceFaces = triPatch.faceFaces();

    const labelList& triFacesInMesh = this->triFacesInMesh();

 	const double Oldtime2=mesh_.time().elapsedCpuTime();
	Info<<"triFacesInMesh Executation Time = "<<Oldtime2-Oldtime0<< " s"<<endl;
 
    label counter = 0;
 
    Info<< "triangles: " << triPatch.size() << " triFacesInMesh: " << returnReduce(triFacesInMesh.size(), sumOp<label>()) << " hit: " << returnReduce(hf.size(), sumOp<label>()) << endl;

    // Only search for tri faces in the mesh
    forAll (triFacesInMesh, tfimI)
    {
        //const label triI = tfimI;
        const label triI = triFacesInMesh[tfimI];

        if (hitTris[triI].size() > 0) //YC Xu 2017/10/11 using surrounding cells instead of only hit points
        {

            addr[triI].setSize(hitTris[triI].size());
            w[triI].setSize(hitTris[triI].size());

            forAll(hitTris[triI], hpI)
            {
               addr[triI][hpI] = hitTris[triI][hpI];
            }

            labelList& curAddr = addr[triI];
            scalarList& curW = w[triI];

            vector curTriCentre = triCentres[triI];

            scalar sumW = 0;

            forAll (curAddr, ibI)
            {
				vector Length = curTriCentre - ibp[curAddr[ibI]];

				Length.z() = 0;
	
	            curW[ibI] = 1.0/(mag(Length)+SMALL);
		
				curW[ibI] *=curW[ibI];

                sumW += curW[ibI];
            }

            // Divide weights by sum distance
            forAll (curW, ibI)
            {
                curW[ibI] /= sumW;
            }				
        }
        else
        {
            // No direct hit.  Start a neighbourhood search

            // Record already visited faces
            labelHashSet visited;

            // Collect new faces to visit
            labelHashSet nextToVisit;

            // Collect new faces to visit
            //SLList<label> nextToVisit;

            // Collect IB points for interpolation
            labelHashSet ibPointsToUse;

            // Initialise with the original tri
            nextToVisit.insert(triI);

			do
			{
				const labelList NTV = nextToVisit.toc();
				nextToVisit.clear();
				forAll(NTV,I)
				{
					const label& next_triI=NTV[I];
					const labelList& hitPointsInNextTriI=hitTris[next_triI];
					if (!visited.found(next_triI))
					{
						if (hitPointsInNextTriI.size() > 0)
						{
							forAll (hitPointsInNextTriI, ptI)
							{
								const label& hitPtI=hitPointsInNextTriI[ptI];
								if(!ibPointsToUse.found(hitPtI))
								{
									ibPointsToUse.insert(hitPtI);
								}
							}
						}						
						visited.insert(hitPointsInNextTriI);
					}
					const labelList& ff=faceFaces[next_triI];
					forAll(ff,ffI)
					{
						const label ff_I=ff[ffI];
						if (!visited.found(ff_I) and !nextToVisit.found(ff_I))
						{
							nextToVisit.insert(ff_I);
						}
					}
				}
            } while
            (
                ibPointsToUse.size() <7
             && !nextToVisit.empty()
            );
			

			/*			
            do
            {
                const label curTri = nextToVisit.removeHead();

                // Discard tri if already visited
                if (visited[curTri]) continue;

                visited.insert(curTri);

                const triFace& curTriPoints = triPatch[curTri];

                // For all current points of face, pick up neighbouring faces
                forAll (curTriPoints, tpI)
                {
                    //the global list index of current point
                    label g_tpI = curTriPoints[tpI];

                    //the local point index of current point
                    label l_tpI = triPatch.whichPoint(g_tpI);

                    if(l_tpI == -1) //didn't find the local point
                    {
                        Pout << "Didn't find the local point index for point "
                             << g_tpI << endl;
                    }

                    // faces sharing current point
                    //const labelList curNbrs = pf[curTriPoints[tpI]];
                    const labelList& curNbrs = pf[l_tpI];

                    //loop over all faces sharing current point
                    forAll (curNbrs, nbrI)
                    {
                        //if current triangle face has not been visited and 
                        //it is inside the domain
                        if (!visited.found(curNbrs[nbrI]))
                        {
                            nextToVisit.append(curNbrs[nbrI]);

                            if (!(hitTris[curNbrs[nbrI]].empty()))
                            {
                                // Found a neighbour with a hit: use all of
                                // its hit points
                                forAll(hitTris[curNbrs[nbrI]],hpI)
                                {
                                  ibPointsToUse.insert(hitTris[curNbrs[nbrI]][hpI]);
                                }
                            }
                        }
                    }
                }
            } while
            (
                ibPointsToUse.size() < 3
             && !nextToVisit.empty()
            );
			*/
            // Found neighbourhood: collect addressing and weights
            addr[triI] = ibPointsToUse.toc();
            w[triI].setSize(addr[triI].size());

            labelList& curAddr = addr[triI];
            scalarList& curW = w[triI];

            vector curTriCentre = triCentres[triI];

            scalar sumW = 0.0;

            forAll (curAddr, ibI)
            {
				vector Length = curTriCentre - ibp[curAddr[ibI]];

				//Length.z() = 0;
	
	            curW[ibI] = 1.0/(mag(Length)+SMALL);
				curW[ibI] *=curW[ibI];
                sumW += curW[ibI];
            }

            // Divide weights by sum distance
            forAll (curW, ibI)
            {
                curW[ibI] /= sumW;
            }
        }
    }

    // Issue a warning if there are triangular faces inside the mesh without
    // neighbouring faces containing ibPoints
    if (counter > 0)
    {
        WarningIn
        (
            "immersedBoundaryFvPatch::makeTriAddressing() const"
        )   << "Not all triangular faces have neighbours with ibPoints" << nl
            << "Number of faces:" << counter << nl
            << "This might cause false force calculation," << nl
            << "consider coarsening the triangular mesh"
            << endl;
    }
 		const double Oldtime1=mesh_.time().elapsedCpuTime();
		Info<<"makeTriAddressing Executation Time = "<<Oldtime1-Oldtime0<< " s"<<endl;
}


const Foam::labelListList&
Foam::immersedBoundaryFvPatch::cellsToTriAddr() const
{

    if (!cellsToTriAddrPtr_)
    {
        makeTriAddressing();
    }

    return *cellsToTriAddrPtr_;
}


const Foam::scalarListList&
Foam::immersedBoundaryFvPatch::cellsToTriWeights() const
{
    if (!cellsToTriWeightsPtr_)
    {
        makeTriAddressing();
    }

    return *cellsToTriWeightsPtr_;
}



// address the trisurface points added by Y. Xu   return local points
void Foam::immersedBoundaryFvPatch::makeTriFacesToTriPoints() const
{
 	const double Oldtime0=mesh_.time().elapsedCpuTime();
    if (debug)
    {
        InfoIn("void immersedBoundaryFvPatch::makeTriFacesToTriPoints() const")
            << "creating tri addressing for immersed boundary " << name()
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (triFacesToTriPointsWeightsPtr_)
    {
        FatalErrorIn("immersedBoundaryFvPatch::makeTriFacesToTriPoints() const")
            << "tri addressing already exist"
            << "for immersed boundary" << name()
            << abort(FatalError);
    }

    const triSurface& triPatch = ibPolyPatch_.ibMesh();
    const vectorField& triCentres = triPatch.faceCentres();
    const vectorField& triPoints = triPatch.localPoints();

    // point-face addressing of trisurface
    const labelListList& triPointFaces = triPatch.pointFaces();

    triFacesToTriPointsWeightsPtr_ = new scalarListList(triPointFaces.size());
    scalarListList& w = *triFacesToTriPointsWeightsPtr_;

    forAll (triPointFaces,ptI)
    {	

        const labelList& curAddr = triPointFaces[ptI];
        scalarList& curW = w[ptI];
		curW.setSize(curAddr.size());

        scalar sumW = 0.0;


        forAll (curAddr, ibI)
        {
            curW[ibI] = 1.0/(mag(triPoints[ptI] - triCentres[curAddr[ibI]])+SMALL);

            sumW += curW[ibI];
        }

            // Divide weights by sum distance
        forAll (curW, ibI)
        {
            curW[ibI] /= sumW;

        }

    }
    
/*this is search two adjacent layer, which may cause problem.
    forAll (triPointFaces,ptI)
    {

		labelHashSet curAddrSet;
		
		const labelList& triPF = triPointFaces[ptI];
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
        vector curTriPoint = triPoints[ptI];

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
*/
 	const double Oldtime1=mesh_.time().elapsedCpuTime();
	Info<<"makeTriFacesToTriPoints Executation Time = "<<Oldtime1-Oldtime0<< " s"<<endl;
}



void Foam::immersedBoundaryFvPatch::makeTriFacesToTriEdges() const
{
 
    if (debug)
    {
        InfoIn("void immersedBoundaryFvPatch::makeTriFacesToTriEdges() const")
            << "creating tri addressing for immersed boundary " << name()
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (triFacesToTriPointsWeightsPtr_)
    {
        FatalErrorIn("immersedBoundaryFvPatch::makeTriFacesToTriEdges() const")
            << "tri addressing already exist"
            << "for immersed boundary" << name()
            << abort(FatalError);
    }

    const pointField& points = ibMesh().localPoints();
    const List<labelledTri>& faces = ibMesh().localFaces();
    const edgeList& edges = ibMesh().edges();
    const labelListList& edgeFaces = ibMesh().edgeFaces();

    triFacesToTriEdgesWeightsPtr_ = new scalarList(ibMesh().nInternalEdges());
    scalarList& weights = *triFacesToTriEdgesWeightsPtr_;

    forAll(weights, edgei)
    {
        vector P = faces[edgeFaces[edgei][0]].centre(points);
        vector N = faces[edgeFaces[edgei][1]].centre(points);
        vector S = points[edges[edgei].start()];
        vector e = edges[edgei].vec(points);

        scalar alpha =
            -(((N - P)^(S - P))&((N - P)^e))/(((N - P)^e )&((N - P)^e));

        vector E = S + alpha*e;

        weights[edgei] = mag(N - E)/(mag(N - E) + mag(E - P));
    }
}

 

const Foam::scalarListList&
Foam::immersedBoundaryFvPatch::triFacesToTriPointsWeights() const
{
    if (!triFacesToTriPointsWeightsPtr_)
    {
        makeTriFacesToTriPoints();
    }

    return *triFacesToTriPointsWeightsPtr_;
}

const Foam::scalarList&
Foam::immersedBoundaryFvPatch::triFacesToTriEdgesWeights() const
{
    if (!triFacesToTriEdgesWeightsPtr_)
    {
        makeTriFacesToTriEdges();
    }

    return *triFacesToTriEdgesWeightsPtr_;
}



// ************************************************************************* //
