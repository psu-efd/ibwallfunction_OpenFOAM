/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description

\*---------------------------------------------------------------------------*/

#include "triSurfaceFieldToTecplot.H"
#include "vtkSurfaceWriter.H"


namespace Foam
{

typedef DynamicList<word> dynWordList;

//scalar field stored in triangular center
void triSurfaceScalarFieldToTecplot
(const triSurface& tsf, const triSurfaceScalarField& tsfSF, const fileName& tpFilePath, const word& tpFileName)
{
    const int nVarmax = 20;
    scalar *TecVar[nVarmax];
    int nvar = 0;
    dynWordList headerNames;

    mkDir(tpFilePath);

    Pout << "tpFilePath = " << tpFilePath << endl;

    fileName tpDataFile(tpFilePath/tpFileName+".dat");
    OFstream dataOut(tpDataFile);
   
    const pointField& points = tsf.points();
    //const faceList& faces(tsf);

    for(int i=0; i<nVarmax; i++)
    {
        TecVar[i] = new scalar [points.size()];
    }

    forAll(points, pointi)
    {
        TecVar[0][pointi] = points[pointi].x();
        TecVar[1][pointi] = points[pointi].y();
        TecVar[2][pointi] = points[pointi].z();
    }
    nvar += 3;
    headerNames.append("X");
    headerNames.append("Y");
    headerNames.append("Z");

    headerNames.append("scalarVar");

    //write out to Tecplot file
    {
      dataOut
             << "TITLE = \"OpenFOAM triSurfaceScalarField \" " << endl;

                dataOut << "VARIABLES = ";
                forAll(headerNames, headerI)
                {
                        dataOut << "\"" << headerNames[headerI] << "\"";
                        if (headerI != headerNames.size()-1)
                        {
                                dataOut << ", ";
                        }
                }

           dataOut << endl
             << "ZONE N=" << points.size() << ", E="
             << tsf.size() <<", DATAPACKING=BLOCK, ZONETYPE=FETRIANGLE" 
             << ", VARLOCATION=([4]=CELLCENTERED)\n";

           //write out the coordinates (must be nodal)
           for(int i=0; i<nvar; i++)
           {
              forAll(points, pointi)
              {
                 dataOut << TecVar[i][pointi] << tab;
              }
              dataOut << endl;
           }

           //write out the scalar varialbe (cell centered)
           forAll(tsf, facei)
           {
                dataOut << tsfSF.field()[facei] << tab;
           }
           dataOut << endl;

           //write out the connectivity list
           forAll(tsf, facei)
           {
                dataOut << tsf[facei][0]+1 << " " 
                        << tsf[facei][1]+1 << " " 
                        << tsf[facei][2]+1 << endl;
           }
           dataOut << endl;
    }

    headerNames.clear();

    for(int i=0; i<nVarmax; i++)
    {
        if(!TecVar[i]) delete TecVar[i];
    }

}

//scalar field stored in triangular center
void triSurfaceScalarFieldToVTK
(const triSurface& tsf, const triSurfaceScalarField& tsfSF, const fileName& vtkFilePath, const word& varName, const word& vtkFileName)
{
    faceList faces(tsf.size());

    forAll(tsf, fI)
    {
        faces[fI] = tsf[fI].triFaceFace();
    }

    vtkSurfaceWriter().write
    (
       vtkFilePath,
       vtkFileName,    // surfaceName
       tsf.points(),
       faces,
       varName, // fieldName
       tsfSF,
       false, // isNodeValues
       true // verbose
    );

}

} // End namespace Foam

