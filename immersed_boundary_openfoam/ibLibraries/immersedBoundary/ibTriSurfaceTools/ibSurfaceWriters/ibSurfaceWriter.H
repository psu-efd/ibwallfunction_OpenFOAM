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

Class
    Foam::ibSurfaceWriter

Description

SourceFiles
    ibSurfaceWriter.C

\*---------------------------------------------------------------------------*/

#ifndef ibSurfaceWriter_H
#define ibSurfaceWriter_H

#include "Field.H"
#include "typeInfo.H"
#include "autoPtr.H"
#include "pointField.H"
#include "faceList.H"
#include "fileName.H"

#include "runTimeSelectionTables.H"
#include "diagTensor.H"
#include "diagTensorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Type> class ibSurfaceWriter;
template<class Type> class nullIbSurfaceWriter;
template<class Type> class proxyIbSurfaceWriter;

/*---------------------------------------------------------------------------*\
                        Class ibSurfaceWriter Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class ibSurfaceWriter
{
    //- friendship between writer data types
    template<class Type2> friend class ibSurfaceWriter;

public:

    //- Runtime type information
    TypeName("ibSurfaceWriter");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            ibSurfaceWriter,
            word,
            (),
            ()
        );


    // Selectors

        //- Return a reference to the selected ibSurfaceWriter
        static autoPtr<ibSurfaceWriter> New(const word& writeFormat);


    // Constructors

        //- Construct null
        ibSurfaceWriter();


    // Destructor

        virtual ~ibSurfaceWriter();


    // Member Functions

        //- Return true if the surface format supports separate files
        virtual bool separateFiles()
        {
            return false;
        }


        //- Writes single surface geometry to file.
        virtual void write
        (
            const fileName& outputDir,      // <root>/<case>/surface/TIME
            const fileName& surfaceName,    // name of surface
            const pointField& points,
            const faceList& faces,
            const bool verbose = false
        ) const
        {}


        //- Writes single surface to file. Either one value per vertex or
        //  one value per face (detected by values.size()==faces.size())
        virtual void write
        (
            const fileName& outputDir,      // <root>/<case>/surface/TIME
            const fileName& surfaceName,    // name of surface
            const pointField& points,
            const faceList& faces,
            const fileName& fieldName,      // name of field
            const Field<Type>& values,
            const bool verbose = false
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ibSurfaceWriter.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
