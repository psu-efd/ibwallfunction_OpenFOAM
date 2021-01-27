  	cpuTime timer;

    boolList cutFace(mesh.nFaces(), false);

    label nCutFaces = 0;

    // Intersect mesh edges with surface (is fast) and mark all faces that
    // use edge.

    forAll(mesh.edges(), edgeI)
    {
        if (debug && (edgeI % 10000 == 0))
        {
            Pout<< "Intersecting mesh edge " << edgeI << " with surface"
                << endl;
        }

        const edge& e = mesh.edges()[edgeI];

        const point& p0 = mesh.points()[e.start()];
        const point& p1 = mesh.points()[e.end()];

        pointIndexHit pHit(querySurf.tree().findLineAny(p0, p1));

        if (pHit.hit())
        {
            const labelList& myFaces = mesh.edgeFaces()[edgeI];

            forAll(myFaces, myFaceI)
            {
                label faceI = myFaces[myFaceI];

                if (!cutFace[faceI])
                {
                    cutFace[faceI] = true;

                    nCutFaces++;
                }
            }
        }
  }

    if (debug)
    {
        Pout<< "Intersected edges of mesh with surface in = "
            << timer.cpuTimeIncrement() << " s\n" << endl << endl;
    }

    //
    // Construct octree on faces that have not yet been marked as cut
    //

    labelList allFaces(mesh.nFaces() - nCutFaces);

    label allFaceI = 0;

    forAll(cutFace, faceI)
    {
        if (!cutFace[faceI])
        {
            allFaces[allFaceI++] = faceI;
        }
    }

    if (debug)
    {
        Pout<< "Testing " << allFaceI << " faces for piercing by surface"
            << endl;
    }

    treeBoundBox allBb(mesh.points());
    // Extend domain slightly (also makes it 3D if was 2D)
    scalar tol = 1e-6 * allBb.avgDim();

    point& bbMin = allBb.min();
    bbMin.x() -= tol;
    bbMin.y() -= tol;
    bbMin.z() -= tol;

    point& bbMax = allBb.max();
    bbMax.x() += 2*tol;
    bbMax.y() += 2*tol;
    bbMax.z() += 2*tol;

    indexedOctree<treeDataFace> faceTree
    (
        treeDataFace(false, mesh, allFaces),
        allBb,      // overall search domain
        8,          // maxLevel
        10,         // leafsize
        3.0         // duplicity
    );

    const triSurface& surf = querySurf.surface();
    const edgeList& edges = surf.edges();
    const pointField& localPoints = surf.localPoints();

    label nAddFaces = 0;
    forAll(edges, edgeI)
    {
      /*  if (edgeI==2616)
        {
            Pout<< "Intersecting surface edge " << edgeI
                << " with mesh faces" << endl;
        }*/
        const edge& e = edges[edgeI];

        const point& start = localPoints[e.start()];
        const point& end = localPoints[e.end()];

        vector edgeNormal(end - start);
        const scalar edgeMag = mag(edgeNormal);
        const vector smallVec = 1e-9*edgeNormal;

        edgeNormal /= edgeMag+VSMALL;

        // Current start of pierce test
        point pt = start;
     /*   if (edgeI==2616)
        {
            Pout<< "        point pt = start; " << pt
                <<endl;
            Pout<< "        end; " << end
                <<endl;
        }
      */  while (true)
        {
            pointIndexHit pHit(faceTree.findLine(pt, end));

            if (!pHit.hit())
            {
                break;
            }
            else
            {
                label faceI = faceTree.shapes().faceLabels()[pHit.index()];

                if (!cutFace[faceI])
                {
                    cutFace[faceI] = true;

                    nAddFaces++;
                }
               // Restart from previous endpoint
                pt = pHit.hitPoint() + smallVec;
     /*  if (edgeI==2616)
        {
            Pout<< "hitPoint " << pt
                <<endl;
            Pout<< "hitPoint1111 " << ((pt-start) & edgeNormal) - edgeMag
                <<endl;
        }
       */         if (((pt-start) & edgeNormal) >= edgeMag)
                {
                    break;
                }
            }
        }
    }

 //   if (debug)
    {
        Pout<< "Detected an additional " << nAddFaces << " faces cut"
            << endl;

        Pout<< "Intersected edges of surface with mesh faces in = "
            << timer.cpuTimeIncrement() << " s\n" << endl << endl;
    }

    Pout<<cutFace.size()<<endl;;
 

