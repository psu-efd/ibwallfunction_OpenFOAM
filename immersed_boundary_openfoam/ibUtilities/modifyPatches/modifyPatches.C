/*---------------------------------------------------------------------------*\
   II   II        II  II   Leichtweiss-Institute for Hydraulics          
  II    II  II  II   II    Dep. Hydromechanics and Coastal Eng. 
 II     IIIIIIII    II     Developed by: Hisham El Safti
IIIIII  II  II     II      Email: hsafti@gmail.com
\*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Description

This utility is modifies "constant/polyMesh/boundary" to change patch types
to empty, wall, cyclic, processor, symmetry or wedge (for use with gmshToFoam)
and to add null patches at end of file (for use with splitMesh) and remove 
existing null patches

Developed by: Hisham El Safti hsafti@gmail.com, January 2013



\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "polyMesh.H"
#include "emptyPolyPatch.H"
#include "wedgePolyPatch.H"
#include "symmetryPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "processorPolyPatch.H"
#include "wallPolyPatch.H"
#include "Time.H"
#include "polyTopoChange.H"
#include "primitiveFacePatch.H"
#include "repatchPolyTopoChanger.H"

#include "immersedBoundaryPolyPatch.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


word typeNameWithoutDash (string currentType)
{
    currentType.erase (0,1);
    return word(currentType);
}


// Add null patch
void addNullPatch(Foam::polyMesh& mesh, const word& name)
{
    const label patchI = mesh.boundaryMesh().findPatchID(name);


    if (patchI != -1 && mesh.boundaryMesh()[patchI].size())
    {
        FatalErrorIn("addNullPatch(const polyBoundaryMesh&, const word&)")
            << "Patch " << name << " is present but non-zero size"
            << exit(FatalError);
    }

    if (patchI == -1)
    {
      Info << "Adding null patch: " << name << endl;

      repatchPolyTopoChanger repatcher(mesh);
      List<polyPatch*> newPatchPtrList((mesh.boundaryMesh().size() + 1));
      forAll(mesh.boundaryMesh(), patchII)
      {
          const polyPatch& patch = mesh.boundaryMesh()[patchII];

          newPatchPtrList[patchII] = patch.clone
                  (
                      mesh.boundaryMesh(),
                      patchII,
                      patch.size(),
                      patch.start()
                  ).ptr();
      }
     
      const immersedBoundaryPolyPatch patch 
          (
              name,
              0,                  // size
              mesh.nFaces(),      // start
              mesh.boundaryMesh().size(),  // index
              mesh.boundaryMesh(),  // polyBoundaryMesh
              "immersedBoundary"
          );

      Info << patch.type() << endl;         

      newPatchPtrList[mesh.boundaryMesh().size()] = patch.clone
          (
              mesh.boundaryMesh(),
              mesh.boundaryMesh().size(),
              patch.size(),
              patch.start()
          ).ptr();

      repatcher.changePatches(newPatchPtrList);    

      repatcher.repatch();
    }
}




// remove null patch
void removeNullPatch(Foam::polyMesh& mesh, const word& name)
{
    const label patchI = mesh.boundaryMesh().findPatchID(name);


    if (patchI == -1 || mesh.boundaryMesh()[patchI].size())
    {
        FatalErrorIn("removeNullPatch(const polyBoundaryMesh&, const word&)")
            << "Patch " << name << " is not present or of non-zero size"
            << exit(FatalError);
    }

    if (patchI != -1)
    {
      Info << "Removing null patch: " << name << endl;

      repatchPolyTopoChanger repatcher(mesh);
      List<polyPatch*> newPatchPtrList((mesh.boundaryMesh().size() - 1));

      label oneLess(0);

      forAll(mesh.boundaryMesh(), patchII)
      {
      
          if (patchI == patchII)
          {
              oneLess=-1;
          } 
          else
          {
              const polyPatch& patch = mesh.boundaryMesh()[patchII];    
      
              newPatchPtrList[patchII+oneLess] = patch.clone
                      (
                          mesh.boundaryMesh(),
                          patchII+oneLess,
                          patch.size(),
                          patch.start()
                      ).ptr();
          }
      }
      
      repatcher.changePatches(newPatchPtrList);    

      repatcher.repatch();
    }
    
}



// Main program:

int main(int argc, char *argv[])
{

    argList::noParallel();
    argList::addOption("null", "patchName", "creates a new patch of zero size");
    argList::addOption("remove", "patchName", "removes a zero size patch");

    argList::removeOption("doc");
    argList::removeOption("srcDoc");


#   include "setRootCase.H"
#   include "createTime.H"
    runTime.functionObjects().off();
#   include "createPolyMesh.H"

    runTime.setTime(instant(runTime.constant()), 0);
    
    word currentType;

    if (argc == 1)
    {
        Info << "\nPlease provide command line arguments or type "
            << "\"modifyPatches -help\" for help\n" << endl; 
        return 0;
    }

    for (label i = 1; i < argc; i++)
    {
        if (word(argv[i]) == "-null" || word(argv[i]) == "-remove")  
        {
            currentType = typeNameWithoutDash((word(argv[i])).c_str());
        }
        else
        {
            if (currentType == "null")
            {
                addNullPatch(mesh, word(argv[i]));
            }
            else if (currentType == "remove")
            {
                removeNullPatch(mesh, word(argv[i]));
            }
            else
            {
                Info << "\nPlease provide command line arguments "
                     << "\n either null or remove\n" << endl;
                return 0;
            }
                
            if (!mesh.write())
            {
                FatalErrorIn(args.executable()) << "Failed writing mesh"
                    << exit(FatalError);
            }
        } 
    }


    Info << "\nModification of patches completed successfully!\n" << endl;     
        
    return 0;

}

// ************************************************************************* //
