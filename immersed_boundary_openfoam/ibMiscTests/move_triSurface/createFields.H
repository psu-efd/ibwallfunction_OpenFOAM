//read in the parameters relevant to the immersed boundary method
    IOdictionary ibmDict
    (
        IOobject
        (
            "ibmDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    word stlFile(ibmDict.lookup("stlFile"));
    vector span(ibmDict.lookup("span"));
    scalar interpR(readScalar(ibmDict.lookup("interpR")));
    label minSurroundCells(readLabel(ibmDict.lookup("minSurroundCells")));

    word interpScheme(ibmDict.lookup("interpScheme"));


