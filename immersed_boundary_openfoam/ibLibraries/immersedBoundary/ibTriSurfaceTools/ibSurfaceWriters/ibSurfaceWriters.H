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

InClass
    Foam::ibSurfaceWriters

Description

\*---------------------------------------------------------------------------*/

#ifndef ibSurfaceWriters_H
#define ibSurfaceWriters_H

#include "ibSurfaceWriter.H"
#include "fieldTypes.H"
#include "diagTensor.H"
#include "diagTensorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Only used internally
#define makeTypeIbSurfaceWritersTypeName(typeWriter, dataType)                 \
                                                                             \
    defineNamedTemplateTypeNameAndDebug(typeWriter< dataType >, 0)

// Sometimes used externally
#define makeIbSurfaceWritersTypeName(typeWriter)                               \
                                                                             \
    makeTypeIbSurfaceWritersTypeName(typeWriter, scalar);                      \
    makeTypeIbSurfaceWritersTypeName(typeWriter, vector);                      \
    makeTypeIbSurfaceWritersTypeName(typeWriter, sphericalTensor);             \
    makeTypeIbSurfaceWritersTypeName(typeWriter, diagTensor);                  \
    makeTypeIbSurfaceWritersTypeName(typeWriter, symmTensor);                  \
    makeTypeIbSurfaceWritersTypeName(typeWriter, tensor);                      \
    //makeTypeIbSurfaceWritersTypeName(typeWriter, symmTensor4thOrder)

// Define type info for single dataType template instantiation (eg, vector)
#define makeIbSurfaceWriterType(typeWriter, dataType)                          \
                                                                             \
    defineNamedTemplateTypeNameAndDebug(typeWriter< dataType >, 0);          \
    addTemplatedToRunTimeSelectionTable                                      \
    (                                                                        \
        ibSurfaceWriter, typeWriter, dataType, word                            \
    )


// Define type info for scalar, vector etc. instantiations
#define makeIbSurfaceWriters(typeWriter)                                       \
                                                                             \
    makeIbSurfaceWriterType(typeWriter, scalar);                               \
    makeIbSurfaceWriterType(typeWriter, vector);                               \
    makeIbSurfaceWriterType(typeWriter, sphericalTensor);                      \
    makeIbSurfaceWriterType(typeWriter, diagTensor);                           \
    makeIbSurfaceWriterType(typeWriter, symmTensor);                           \
    makeIbSurfaceWriterType(typeWriter, tensor);                               \
    //makeIbSurfaceWriterType(typeWriter, symmTensor4thOrder)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
