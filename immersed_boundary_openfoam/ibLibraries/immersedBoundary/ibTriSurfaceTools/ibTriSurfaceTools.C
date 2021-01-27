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

#include "ibTriSurfaceTools.H"
#include "addToRunTimeSelectionTable.H"
#include "vectorTools.H"
#include "SortableList.H"
#include "volFields.H"
#include "pointToPointPlanarInterpolation.H"
#include "triSurfaceMesh.H"
#include "triSurfaceFieldToTecplot.H"
#include "Time.H"
#include "fvMesh.H"
#include "simpleMatrix.H"
#include "PrimitivePatchInterpolation.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ibTriAddressing.C"
#endif 


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


//sand-slide correction algorithm added by Y. Xu
void Foam::ibTriSurfaceTools::makeSandSlide(const fvMesh& mesh, const bool& outside_ = false)
{
    triSurface& triS_(triSurface_);

    tmp<vectorField > tcf(new vectorField(triS_.faceCentres()));
    vectorField& cf = tcf();// face centres

    tmp<scalarField > tzb(new scalarField(cf.size()));
    scalarField& zb = tzb();// face centres z coordinate

    tmp<vectorField > tsf(new vectorField(cf.size()));
    vectorField& sf = tsf();// face area vector

    tmp<scalarField > tAh(new scalarField(cf.size()));
    scalarField& Ah = tAh();// projected face area

    SortableList<scalar> theta0(Ah.size());
 
    const labelListList& ffAddress(triS_.faceFaces()); // addressing cell neighbourings

    // initialize prjected face area and z cooridnate
    forAll(Ah,triI)
    {
		sf[triI] = triS_[triI].normal(triS_.points());
        Ah[triI] = mag(sf[triI]&(g_/mag(g_))); // Ah changes as zb changes
        zb[triI] = cf[triI].z();
        theta0[triI] = mag(Foam::vectorTools::radAngleBetween(sf[triI],-g_));
    }

    //sort label list according to the slope angle
    labelList newList (Ah.size());// face centres

    sortedOrder(theta0,newList);

    tmp<scalarField > tcts(new scalarField(Ah.size(),0));
    scalarField& cts = tcts();

    //iteration
    forAll(Ah,I)
    {
        label triI=newList[I];
        const labelList& pAddr = ffAddress[triI];
        scalar dz=0;// dz=dz1/(Ah-dz2)

        scalar totalMass=0; //check the mass conservation of the bathymetry      

		if(!inMesh(mesh,triI))
		{
			continue;
		}

		calc_dz
		(
			mesh,
			cf,
			pAddr,
			Ah,
    		triI,
			dz,
			totalMass,
			outside_
		);

        if(mag(totalMass)>1e-5) //check the mass conservation of the bathymetry
        {
			cts[triI]=totalMass;
            Info<<triI<<" "<<cf[triI]<<" "<<totalMass<<endl;
        }
    }
	//Info<<"sum of loss "<<cts.size()<<"/"<<zb.size()<<" "<<average(cts)<<" "<<max(cts)<<endl;

    forAll(cf,triI)
    {
        zb[triI]=cf[triI].z();
    }

    tmp<pointField > tpts(new pointField(triS_.localPoints()));
    pointField& pts = tpts(); // local points

	faceList faces(ffAddress.size());

	List<labelledTri> triFaces(triS_.localFaces());
	forAll(triFaces,triI)
	{
		faces[triI]=triFaces[triI];
	}
	PrimitivePatch<face, ::Foam::List, pointField, point>  pp(faces,pts);
	PrimitivePatchInterpolation< PrimitivePatch<face, ::Foam::List, pointField, point> > ppI(pp);


    tmp<scalarField > tptzb(new scalarField(toTriPoints(zb)));
    scalarField& ptzb = tptzb();

	//scalarField ptzb = ppI.faceToPointInterpolate(zb);
    //scalarField ptzb(toTriPoints(zb)); // interpolate from cell center to point

    forAll(pts,ptI)
    {
		if(mesh.findCell(pts[ptI]))
		{
			pts[ptI].z()=ptzb[ptI];
		}
    }

	//change local points label list to global points label list
    tmp<pointField > tglobal_pts(new pointField(triS_.points()));
    pointField& global_pts = tglobal_pts();//global points

    forAll(global_pts,gp)
    {
        label ptI=triS_.whichPoint(gp); //global points to local points
        global_pts[gp]=pts[ptI];
    }

    triS_.movePoints(global_pts); //have to use global pts which from points()

	triInMesh_.clear();
    triPointsToHitPointsAddrPtr_.clear();
    triPointsToHitPointsWeightsPtr_.clear();
	//move outside points
}

void Foam::ibTriSurfaceTools::calc_dz
(
	const fvMesh& mesh,
	vectorField& cf,
    const labelList& pAddr,
    scalarField Ah,
    label triI,
	scalar& dz,
	scalar& totalMass,
	const bool& outside_ = false
)
{
    scalar Cp=0;
    scalar Ci=0;

    tmp<scalarField > ttheta(new scalarField(pAddr.size(),0));
    scalarField& theta = ttheta();

    tmp<scalarField > tSpi(new scalarField(pAddr.size(),0));
    scalarField& Spi = tSpi();

    tmp<scalarField > tLpi(new scalarField(pAddr.size(),0));
    scalarField& Lpi = tLpi();

   	forAll(pAddr,triII)
   	{
		const label& pAddrI=pAddr[triII];
	    theta[triII]=Foam::vectorTools::radAngleBetween(cf[triI]-cf[pAddrI],-g_);

        if(radToDeg(theta[triII])>90)
        {
            theta[triII]=theta[triII]-degToRad(90);
            Spi[triII]=-Foam::tan(theta[triII]);
        }
        else
        {
            theta[triII]=degToRad(90)-theta[triII];
            Spi[triII]=Foam::tan(theta[triII]);
        }

        Lpi[triII]=mag(cf[triI]-cf[pAddrI])*mag(Foam::cos(theta[triII]));

		if(inMesh(mesh,pAddrI) and inMesh(mesh,triII))
		{
	        if(mag(Spi[triII])>repose_)
	        {
	            Spi[triII]=repose_*Spi[triII]/mag(Spi[triII]);
	        }
            Ci+=Ah[pAddrI]*(cf[triI].z()-cf[pAddrI].z()-Lpi[triII]*Spi[triII]);
	        Cp+=Ah[pAddrI];
		}
		else if(!inMesh(mesh,pAddrI) and inMesh(mesh,triII))
		{
			Spi[triII]=0;
        	Ci+=Ah[pAddrI]*(cf[triI].z()-cf[pAddrI].z()-Lpi[triII]*Spi[triII]);
	        Cp+=Ah[pAddrI];
		}
		else if(inMesh(mesh,pAddrI) and !inMesh(mesh,triII))
		{
			Spi[triII]=0;
        	Ci+=Ah[pAddrI]*(cf[triI].z()-cf[pAddrI].z()-Lpi[triII]*Spi[triII]);
	        Cp+=Ah[pAddrI];
		}
		else if(!inMesh(mesh,pAddrI) and !inMesh(mesh,triII))
		{
	        if(mag(Spi[triII])>repose_)
	        {
	            Spi[triII]=repose_*Spi[triII]/mag(Spi[triII]);
	        }
            Ci+=Ah[pAddrI]*(cf[triI].z()-cf[pAddrI].z()-Lpi[triII]*Spi[triII]);
	        Cp+=Ah[pAddrI];
		}
 	}

    if(Ah[triI]-Cp==0)
   	{
   	    dz=0;
   	}
   	else
   	{
       	dz=Ci/(Ah[triI]-Cp);
   	}
    forAll(pAddr,triII)
   	{
		const label& pAddrI=pAddr[triII];
		scalar dzi=cf[triI].z()-cf[pAddrI].z()-Lpi[triII]*Spi[triII];
        totalMass-=Ah[pAddrI]*dzi;
        cf[pAddrI].z()+=dzi;
    }

	cf[triI].z()+=dz;

    totalMass+=Ah[triI]*dz;
}


void Foam::ibTriSurfaceTools::checkSandSlide(const fvMesh& mesh, const bool& outside_ = false)
{

    triSurface& triS_(triSurface_);

    bool& ssp_=sandSlidePtr_;
    ssp_=false;

    const vectorField& cf(triS_.faceCentres());// face centres

    const labelListList& ffAddress(triS_.faceFaces()); // addressing cell neighbourings

    //iteration
    forAll(cf,triI)
    {
        const labelList& pAddr = ffAddress[triI];
		tmp<scalarField > ttheta(new scalarField(pAddr.size(),0));
		scalarField& theta = ttheta();

		tmp<scalarField > tSpi(new scalarField(pAddr.size(),0));
		scalarField& Spi = tSpi();
		bool inMesh_ = true;

		if(!inMesh(mesh,triI))
		{
			continue;
		}

        forAll(pAddr,triII)
        {
			const label& pAddrI=pAddr[triII];
			if(!inMesh(mesh,pAddrI))
			{
				inMesh_ = false; // at least one point of triangle is not inside the mesh
			}
		}
 
		if(!inMesh_)
		{
			continue;
		}
 
        forAll(pAddr,triII)
        {
			const label& pAddrI=pAddr[triII];

			if(!inMesh(mesh,pAddrI))
			{
				continue;
			}

            theta[triII]=Foam::vectorTools::radAngleBetween(cf[triI]-cf[pAddrI],-g_);

            if(radToDeg(theta[triII])>90)
            {
                theta[triII]=theta[triII]-degToRad(90);
                Spi[triII]=-Foam::tan(theta[triII]);
            }
            else
            {
                theta[triII]=degToRad(90)-theta[triII];
                Spi[triII]=Foam::tan(theta[triII]);
            }

        	if(mag(Spi[triII])>repose_)
        	{
            	ssp_=true;
        	}
        }
    }
}

bool Foam::ibTriSurfaceTools::inMesh
(
	const fvMesh& mesh,
	const label& triI
)

{
	if(triInMesh_.empty())
	{
		const vectorField& fc=triSurface_.faceCentres();
	    boolList& triInMesh = triInMesh_;
		triInMesh.setSize(fc.size());
		forAll(fc, fI)
		{
			if(mesh.findCell(fc[fI])>-1)
			{
				triInMesh[fI] = true;
			}
			else
			{
				triInMesh[fI] = false;
			}
		}
	}
	return triInMesh_[triI];
}




// * * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * //

//Haudsorff distance, distance between two data sets
Foam::scalar Foam::ibTriSurfaceTools::Hausdorff_dist
(
    const pointField& pts1,
    const pointField& pts2
)
{
    scalar max_1=0;
    forAll(pts1,I)
    {
        scalar min_=1000;
        forAll(pts2,II)
        {
            scalar dist_(mag(pts1[I]-pts2[II]));
            if(dist_<min_)
            {
                min_=dist_;
            }
        }
        if(min_>max_1)
        {
            max_1=min_;
        }
    }
    scalar max_2=0;
    forAll(pts2,I)
    {
        scalar min_=1000;
        forAll(pts1,II)
        {
            scalar dist_(mag(pts2[I]-pts1[II]));
            if(dist_<min_)
            {
                min_=dist_;
            }
        }
        if(min_>max_2)
        {
            max_2=min_;
        }
    }
	return max(max_1,max_2);
}



void Foam::ibTriSurfaceTools::sandSlide
(
    triSurface& surf,
    const fvMesh& mesh,
    vector& gD,
    scalar& reposeD
    
)
{
    triSurface_=surf;
    g_=gD;
    repose_=reposeD*0.95;
    massRatio_=1.0;
    triSurface& triS_(triSurface_);

	//typedef List<triSurface> triSurfaceList;
	//triSurfaceList procTriSurfaces(Pstream::nProcs());

	//if (Pstream::myProcNo() == Pstream::masterNo())
    {
        bool& ssp_=sandSlidePtr_;
		Info<<"Start initial sand-slide check"<<endl;
 		const double Oldtime0=mesh.time().elapsedCpuTime();

        checkSandSlide(mesh,false);
 		const double Oldtime1=mesh.time().elapsedCpuTime();
		Info<<"checkSandSlide Executation Time = "<<Oldtime1-Oldtime0<< " s"<<endl;
        label counts=0;

        const pointField& cc(triS_.faceCentres());//new face centers

		tmp<vectorField > tsf(new vectorField(cc.size()));
		vectorField& sf = tsf();// face area vector

		tmp<scalarField > tAh(new scalarField(cc.size()));
		scalarField& Ah = tAh();// projected face area

        scalar totalMass=0;

        forAll(cc,triI)
        {
			if(!inMesh(mesh,triI))
			{
				continue;
			}

			sf[triI] = triS_[triI].normal(triS_.points());
	        Ah[triI] = mag(sf[triI]&(g_/mag(g_))); // Ah changes as zb changes
            totalMass +=(cc[triI].z()+0.5)*Ah[triI];
        }
        scalar tMass0=totalMass;
        //Info <<tMass0<<" "<<totalMass/tMass0<<endl;

        while(ssp_ && counts<150)
        {
			Info<<counts<<endl;
            makeSandSlide(mesh,false);
            //makeSandSlide_outside(mesh);
            checkSandSlide(mesh,false);
            counts+=1;

            //Info<<counts<<endl;
            //word cts(Foam::Time::timeName(counts, 5));
            //fileName name_dym="collapse_"+cts+".stl";
            //Pout<<name_dym<<endl;
            //rm(name_dym);
            //triS_.write(name_dym);

            const pointField& newcc(triS_.faceCentres());//new face centers
			tmp<vectorField > tnewsf(new vectorField(newcc.size(), pTraits<vector>::zero));
			vectorField& newsf = tnewsf();// face area vector

			tmp<scalarField > tnewAh(new scalarField(newcc.size(),0));
			scalarField& newAh = tnewAh();// projected face area

			tmp<scalarField > tnewtheta(new scalarField(newsf.size(),0));
			scalarField& newtheta = tnewtheta();//slope angle

            totalMass=0;    
            forAll(newcc,triI)
            {
			
				if(!inMesh(mesh,triI))
				{				
					continue;
				} 				
				newsf[triI] = triS_[triI].normal(triS_.points());
                newAh[triI] = mag(newsf[triI]&(g_/mag(g_)));
                totalMass +=(newcc[triI].z()+0.5)*newAh[triI];
 
                newtheta[triI]=radToDeg(Foam::vectorTools::radAngleBetween(newsf[triI],-g_));
            }
            massRatio_=totalMass/tMass0;
           // Info <<totalMass<<" "<<massRatio_<<endl;
           // Info <<triS_.faceCentres().size()<<" "<<triS_.points().size()<<endl;
            //writeTriValue(newtheta,mesh,"theta",cts);
			//writeTriValue(newAh,mesh,"newAh",cts);
            Info <<"Sand-slide correction: iteration "<<counts<< " Mass conserv. "<<massRatio_<<endl;
        }
		//clearDelete();
		//triSurface& triS_(ibMeshRef());
		//fileName name_dym="constant/triSurface/sandSlide_morph"+mesh_.time().timeName()+".stl";
		//rm(name_dym);
		//triS_.write(name_dym);
		 //   Info <<"sand-slide correction: No Iterations "<<counts<<endl;
		surf=triSurface_;
		//procTriSurfaces[Pstream::myProcNo()]=triSurface_;
    }/*
	else
	{
		procTriSurfaces[Pstream::myProcNo()]=triSurface_;
	}

   	Pstream::gatherList(procTriSurfaces);
   	Pstream::scatterList(procTriSurfaces);
Info<<"dddddddddddddddddddddddd"<<endl;
	if(Pstream::myProcNo() != Pstream::masterNo())
	{
		surf=procTriSurfaces[Pstream::masterNo()];
	}*/
}




void Foam::ibTriSurfaceTools::writeTriValue
(
    const scalarField triangularValues, 
    const fvMesh& mesh,
    const word& varName, 
    const word& countName
) const
{

    //Foam::Time runTime(Foam::Time::controlDictName);
 

    triSurfaceMesh stlMesh
    (
         IOobject
         (
          "stlMesh",
          mesh.time().timeName(),
          mesh.time(),
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
         ),
         triSurface_
	);
 
    
    triSurfaceScalarField stlField
    (
         IOobject
         (
          "stlField",
          mesh.time().timeName(),
          mesh.time(),
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
         ),
         stlMesh,
         dimensionSet(0, 0,0, 0, 0, 0, 0),
         triangularValues
     );

     // fill the values to the triangle's centres
     forAll(stlField, triI)
     {
        stlField[triI] = triangularValues[triI];
     }
     
     fileName tpPath=".";   
/*    if (Pstream::parRun())
    {
        tpPath = mesh_.time().path()/".."/"triValues";
    }
    else
    {
        tpPath = mesh_.time().path()/"triValues";
    }
   if (!isDir(tpPath))
    {
        mkDir(tpPath);
    }
*/ 
     triSurfaceScalarFieldToVTK
           (stlMesh, stlField, tpPath, varName, countName);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //



