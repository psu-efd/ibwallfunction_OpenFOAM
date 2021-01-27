const boolList& piercedFace(cutFace);

List<cellInfo> cellInfoList(mesh.nCells());

    // Mark cut cells first
    forAll(piercedFace, faceI)
    {
        if (piercedFace[faceI])
        {
            cellInfoList[mesh.faceOwner()[faceI]] =
                cellInfo(cellClassification::CUT);

            if (mesh.isInternalFace(faceI))
            {
                cellInfoList[mesh.faceNeighbour()[faceI]] =
                    cellInfo(cellClassification::CUT);
            }
        }
    }

    //
    // Mark cells containing outside points as being outside
    //

    // Coarse guess number of faces
    labelHashSet outsideFacesMap(outsidePts.size() * 6 * 2);

   forAll(outsidePts, outsidePtI)
    {
        // Use linear search for points.
        label cellI = queryMesh.findCell(outsidePts[outsidePtI], -1, false);

        if (returnReduce(cellI, maxOp<label>()) == -1)
        {
            FatalErrorIn
            (
                "List<cellClassification::cType> markCells"
                "(const meshSearch&, const boolList&, const pointField&)"
            )   << "outsidePoint " << outsidePts[outsidePtI]
                << " is not inside any cell"
                << nl << "It might be on a face or outside the geometry"
                << exit(FatalError);
        }

        if (cellI >= 0)
        {
            cellInfoList[cellI] = cellInfo(cellClassification::OUTSIDE);

            // Mark faces of cellI
            const labelList& myFaces = mesh.cells()[cellI];
            forAll(myFaces, myFaceI)
            {
                outsideFacesMap.insert(myFaces[myFaceI]);
            }
        }
    }
    //
    // Mark faces to start wave from
    //

    labelList changedFaces(outsideFacesMap.toc());

    List<cellInfo> changedFacesInfo
    (
        changedFaces.size(),
        cellInfo(cellClassification::OUTSIDE)
    );
   
    int dummyTrackData_ = mesh.nCells();

    MeshWave<cellInfo> cellInfoCalc
    (
        mesh,
        changedFaces,                       // Labels of changed faces
        changedFacesInfo,                   // Information on changed faces
        cellInfoList,                       // Information on all cells
        mesh.globalData().nTotalCells()+1  // max iterations
        ,dummyTrackData_
    );

Pout<<"mesh.globalData().nTotalCells()+1 "<<mesh.globalData().nTotalCells()+1<<endl;
    // Get information out of cellInfoList
    const List<cellInfo>& allInfo = cellInfoCalc.allCellInfo();
Pout<<"cellInfoCalc.getUnsetCells().size() "<<cellInfoCalc.getUnsetCells()<<endl;
Pout<<"cellInfoCalc.getUnsetCells().data() "<<cellInfoCalc.data()<<endl;
    forAll(allInfo, cellI)
    {
        label t = allInfo[cellI].type();

        if (t == cellClassification::NOTSET)
        {
        //Pout<<t<<endl;
            t = cellClassification::INSIDE;
        }
        cellType[cellI]= t;

    }


