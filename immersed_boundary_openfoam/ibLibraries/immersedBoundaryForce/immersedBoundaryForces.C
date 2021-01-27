/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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

\*---------------------------------------------------------------------------*/

#include "immersedBoundaryForces.H"
#include "immersedBoundaryFvPatch.H"
#include "immersedBoundaryFvPatchFields.H"
#include "immersedBoundaryVelocityWallFunctionFvPatchVectorField.H"

#include "volFields.H"
#include "surfaceFields.H"
#include "dictionary.H"
#include "Time.H"

#include "itoa.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(immersedBoundaryForces, 0);
}

// * * * * * * * * * * * * * * * Private member functions  * * * * * * * * * //
void Foam::immersedBoundaryForces::makeFMList() const
{
    // Only one patch is allowed in this patchSet_. However, this patch
    // could contain multiple unconnected objects.
    if(patchSet_.size()!=1)
    {
        FatalErrorIn("immersedBoundaryForces::makeFile")
            << " Only one patch is allowed in the patchSet. "
            << " patchSet_.size() = " << patchSet_.size() << nl
            << abort(FatalError);
    }

    const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);
    //const fvMesh& mesh = U.mesh();

    // The following loop will only loop once since there is only one 
    // patch in the patchSet_.
    // for each immersedBoundary patch in the domain
    forAllConstIter(labelHashSet, patchSet_, iter)
    {
         label patchI = iter.key();

         if
         (
            isA<immersedBoundaryFvPatchVectorField>
            (
                U.boundaryField()[patchI]
            )
         )
         {
             // Found immersed boundary patch and field.
             // Cast into immersed boundary type
      /*       const immersedBoundaryFvPatch& ibPatch =
                 refCast<const immersedBoundaryFvPatch>
                 (
                     mesh.boundary()[patchI]
                 );

             // Get how many unconnected objects (zones) associated
             // with this immersedBoundaryPatch
             label numZones = ibPatch.numZones();

             if(numZones <= 0)
             {
                  FatalErrorIn("immersedBoundaryForces::makeFile()")
                      << " number of unconnected objects (zones) insanity"
                      << " numZones = " << numZones << nl
                      << abort(FatalError);
             }

             // Create the forcesMoments for each object if not already created
             if (debug)
             {
                  Pout<< "Creating fmAllPtrs_." << endl;
             }*/
             label numZones = 1;
             fmAllPtrs_ = new PtrList<forcesMoments>(numZones);
             PtrList<forcesMoments>& fmIdm = *fmAllPtrs_;

             forAll(fmIdm, zoneI)
             {
                  fmIdm.set
                  (
                     zoneI,
                     new forcesMoments
                         (
                            pressureViscous(vector::zero, vector::zero),
                            pressureViscous(vector::zero, vector::zero)
                         )
                  );
             }
        }
    }
}


void Foam::immersedBoundaryForces::makeFile()
{
    // Only one patch is allowed in this patchSet_. However, this patch
    // could contain multiple unconnected objects.
    if(patchSet_.size()!=1)
    {
        FatalErrorIn("immersedBoundaryForces::makeFile")
            << " Only one patch is allowed in the patchSet. "
            << " patchSet_.size() = " << patchSet_.size() << nl
            << abort(FatalError);
    }

    const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);
    //const fvMesh& mesh = U.mesh();

    // The following loop will only loop once since there is only one 
    // patch in the patchSet_.
    // for each immersedBoundary patch in the domain
    forAllConstIter(labelHashSet, patchSet_, iter)
    {
         label patchI = iter.key();

         if
         (
            isA<immersedBoundaryFvPatchVectorField>
            (
                U.boundaryField()[patchI]
            )
         )
         {
             // Found immersed boundary patch and field.
             // Cast into immersed boundary type
  /*           const immersedBoundaryFvPatch& ibPatch =
                 refCast<const immersedBoundaryFvPatch>
                 (
                     mesh.boundary()[patchI]
                 );

             // Get how many unconnected objects (zones) associated
             // with this immersedBoundaryPatch
             //label numZones = ibPatch.numZones();

             if(numZones <= 0)
             {
                  FatalErrorIn("immersedBoundaryForces::makeFile()")
                      << " number of unconnected objects (zones) insanity"
                      << " numZones = " << numZones << nl
                      << abort(FatalError);
             }

             // Create the force file for each object
             if (debug)
             {
                      Info<< "Creating force files." << endl;
             }
*/
             label numZones = 1;
             filePtrs_ = new PtrList<OFstream>(numZones);
             PtrList<OFstream>& fileIdm = *filePtrs_;

             forAll(fileIdm, zoneI)
             {

                 // File output: only master will do this 
                 // job (if run in parallel)
                 if (Pstream::master())
                 {
                      fileName forcesDir;
                      word startTimeName =
                         obr_.time().timeName(obr_.time().startTime().value());

                      if (Pstream::parRun())
                      {
                         // Put in undecomposed case (Note: gives problems for
                         // distributed data running)
                         forcesDir = 
                            obr_.time().path()/".."/name_/startTimeName;
                      }
                      else
                      {
                         forcesDir = 
                            obr_.time().path()/name_/startTimeName;
                      }

                      // Create directory if does not exist.
                      mkDir(forcesDir);

                      // Open new file at start up
                      fileIdm.set
                      (
                        zoneI,
                        new OFstream(forcesDir/("obj" + itoa(zoneI) + ".dat"))
                      );

//                    Pout << "zoneI = " << zoneI << " " << "forcesDir = " << forcesDir << endl;

                      // Add headers to output data
                      // writeFileHeader();
                      fileIdm[zoneI]
                        << "# Time" << tab
                        << "forces(pressure, viscous) moment(pressure, viscous)"
                        << endl;
                 }//end of master
            }//end of for each file (object)
       }//end of for current ibPatch
   }
}

void Foam::immersedBoundaryForces::write()
{
    if (active_)
    {
        // Create the forces file if not already created
        if (!filePtrs_)
        {
           makeFile();
        }

        // calculate force and moment on all objects
        calcForcesMomentAll();

        forAll(*fmAllPtrs_, zoneI)
        {
           forcesMoments& fm = (*fmAllPtrs_)[zoneI];

           //Only master will do this job (if run in parallel)
           if (Pstream::master())
           {
              ((*filePtrs_)[zoneI]) << obr_.time().value() << tab << fm << endl;

              if (log_)
              {
                  Info<< "forces output:" << nl
                      << "    forces(pressure, viscous)" << fm.first() << nl
                      << "    moment(pressure, viscous)" << fm.second() << nl
                      << endl;
              }
           }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immersedBoundaryForces::immersedBoundaryForces
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    forces
    (
        name,
        obr,
        dict,
        loadFromFiles
    ),
    filePtrs_(NULL),
    fmAllPtrs_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::immersedBoundaryForces::~immersedBoundaryForces()
{
    deleteDemandDrivenData(filePtrs_);
    deleteDemandDrivenData(fmAllPtrs_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
//Calculate the list of forcesMoments for all objects in the STL file.
void Foam::immersedBoundaryForces::calcForcesMomentAll() const
{
    //Allocate memory at first entry
    if(!fmAllPtrs_)
    {
        //Pout << "forcesMoments list empty. Make a new one." << endl;
        makeFMList();      
    }

    if (directForceDensity_)
    {
          //Not relevant. 
    }
    else
    {
        volVectorField& U = 
           const_cast<volVectorField&>(obr_.lookupObject<volVectorField>(UName_));
        volScalarField& p = 
           const_cast<volScalarField&>(obr_.lookupObject<volScalarField>(pName_));

        const fvMesh& mesh = U.mesh();

//Pout << "test1" << endl;
        volSymmTensorField stress = devRhoReff();

        // Scale pRef by density for incompressible simulations
        // Why we want this?? Isn't that for incompressible simulations
        // pRef_ is already kinematic, i.e. p/rho?
        scalar pRef = pRef_/rho(p);


        // for each immersedBoundary patch in the domain
        // In fact, only one ib patch is allowed. 
        forAllConstIter(labelHashSet, patchSet_, iter)
        {
            label patchI = iter.key();

            // Check and cast into immersed boundary type
            if
            (
                isA<immersedBoundaryFvPatchVectorField>
                (
                    U.boundaryField()[patchI]
                )
            )
            {
                // Found immersed boundary patch and field.
                // Cast into immersed boundary type
                const immersedBoundaryFvPatch& ibPatch =
                    refCast<const immersedBoundaryFvPatch>
                    (
                        mesh.boundary()[patchI]
                    );

                immersedBoundaryFvPatchSymmTensorField stressPatch =
                    refCast<immersedBoundaryFvPatchSymmTensorField>
                    (
                        stress.boundaryField()[patchI]
                    );
//Pout << "stressPatch = " << stressPatch << endl;

                immersedBoundaryFvPatchScalarField pPatch =
                    refCast<immersedBoundaryFvPatchScalarField>
                    (
                        p.boundaryField()[patchI]
                    );

                stressPatch.updateCoeffs();
                pPatch.updateCoeffs();
                //stressPatch.imposeImmersedBoundaryEffect();
                //pPatch.imposeImmersedBoundaryEffect();


//                pPatch.writeTriValue("p"); 

                // Get how many unconnected objects (zones) associated
                // with this immersedBoundaryPatch
//                label numZones = ibPatch.numZones();
//Pout << "numZones = " << numZones << endl;
/*
                if(numZones <= 0)
                {
                  FatalErrorIn("immersedBoundaryForces")
                      << " number of unconnected objects (zones) insanity"
                      << " numZones = " << numZones << nl
                      << abort(FatalError);
                }
*/
                label numZones = 1;
   
                // Get ibPatch data on the whole surface mesh
                const tmp<scalarField> pTriF = pPatch.triValue();
 
                const vectorField& Sfb = ibPatch.triSf();

                // Get triangular surface area vectors
                //const tmp<vectorField> tSfbTriFInM = ibPatch.renumberField(Sfb);
                //const vectorField& SfbTriFInM = tSfbTriFInM();
                const scalarField sA = mag(Sfb);
 



   /*             // Calculate the face area for triangles
                PtrList<scalarField> sAPtr(numZones);

                forAll(sAPtr, zoneI)
                {
                     sAPtr.set
                     (
                        zoneI,
                        new scalarField
                        (
                           mag(SfbPtr[zoneI])
                        )
                     );
                }


                // Calculate distance for triangles
                PtrList<vectorField> MdPtr(numZones);

                forAll(MdPtr, zoneI)
                {
                     MdPtr.set
                     (
                        zoneI,
                        new vectorField
                        (
//                           ibPatch.triCf()[zoneI]-CofR_
                           ibPatch.triCf()-coordSys_.origin()
                        )
                     );
                }
 */
                // For each of the unconnected object, calculate the force
                // and moments
                for(label zone = 0; zone < numZones; zone++) 
                {
                   forcesMoments& fm = (*fmAllPtrs_)[zone];

                   // Reset the values of forcesMoments to zero
                   fm.first().first() = vector::zero;   
                   fm.second().first() = vector::zero;   
                   fm.first().second() = vector::zero;   
                   fm.second().first() = vector::zero;   

                   // Pressure force is an integral of interpolated pressure
                   // on triangular faces
                   vectorField pf =
                       -Sfb*(pTriF() - pRef);
                   
                   fm.first().first() += rho(p)*gSum(pf);
                   fm.second().first() += rho(p)*gSum((ibPatch.triCf()-coordSys_.origin()) ^ pf);

                if
                (
                    isA<immersedBoundaryVelocityWallFunctionFvPatchVectorField>
                    (   
                        U.boundaryField()[patchI]
                    )
                )
                {

                   // Shear force is dotted with a normal in the IB point
                   // and integrated on triangular faces
                    const
                    immersedBoundaryVelocityWallFunctionFvPatchVectorField&
                        UPatch = refCast
                        <
                            const
                            immersedBoundaryVelocityWallFunctionFvPatchVectorField
                        >
                        (   
                            U.boundaryField()[patchI]
                        );

 
                    vectorField vf =
                        -sA*ibPatch.toTriFaces(UPatch.wallShearStress());
 
                   fm.first().second() += gSum(vf);
                   fm.second().second() += gSum((ibPatch.triCf()-coordSys_.origin()) ^ vf);
                }
                else if
                (
                    isA<immersedBoundaryFvPatchVectorField>
                    (
                        U.boundaryField()[patchI]
                    )
                )
                {
                    // Laminar flow

                    // Get immersed boundary velocity
                    const immersedBoundaryFvPatchVectorField& UPatch =
                        refCast<const immersedBoundaryFvPatchVectorField>
                        (
                            U.boundaryField()[patchI]
                        );

                    // Look up the viscosity
                    if (mesh.foundObject<dictionary>("transportProperties"))
                    {
                        const dictionary& transportProperties =
                            mesh.lookupObject<dictionary>
                            (
                                "transportProperties"
                            );

                        dimensionedScalar nu(transportProperties.lookup("nu"));

                        // Integrate wall shear stress on triangular
                        // faces and get the part inside the mesh
                        const tmp<vectorField> grad = UPatch.triGrad();

                        vectorField vf =
                            -sA*rho(p)*nu.value()*grad;

                        fm.first().second() += gSum(vf);
                        fm.second().second() += gSum((ibPatch.triCf()-coordSys_.origin()) ^ vf);
                    }
                else
                    {
                        InfoIn
                        (
                            "immersedBoundaryForces::forcesMoments"
                            "immersedBoundaryForces::calcForcesMoment() const"
                        )   << "Laminar flow, but cannot find nu.  Skipping"
                            << endl;
                    }


                }


                }
            }
        }
    }
}


// ************************************************************************* //
