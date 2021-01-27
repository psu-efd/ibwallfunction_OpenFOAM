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

#include "ibSurfaceWriter.H"

#include "MeshedSurfaceProxy.H"
#include "nullIbSurfaceWriter.H"
#include "proxyIbSurfaceWriter.H"

#include "HashTable.H"
#include "word.H"

template<>
Foam::ibSurfaceWriter<bool>::wordConstructorTable*
Foam::ibSurfaceWriter<bool>::wordConstructorTablePtr_;


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr< Foam::ibSurfaceWriter<Type> >
Foam::ibSurfaceWriter<Type>::New(const word& writeType)
{
    typename wordConstructorTable::iterator cstrIter =
        wordConstructorTablePtr_->find(writeType);

    if (cstrIter == wordConstructorTablePtr_->end())
    {
        // not supported for this data type, but it generally does work
        // (it handles the 'bool' specialization - ie, geometry write)
        if
        (
            Foam::ibSurfaceWriter<bool>::wordConstructorTablePtr_->found
            (
                writeType
            )
        )
        {
            // use 'null' handler instead
            return autoPtr< ibSurfaceWriter<Type> >
            (
                new nullIbSurfaceWriter<Type>()
            );
        }
        else if (MeshedSurfaceProxy<face>::canWriteType(writeType))
        {
            // generally unknown, but can be written via MeshedSurfaceProxy
            // use 'proxy' handler instead
            return autoPtr< ibSurfaceWriter<Type> >
            (
                new proxyIbSurfaceWriter<Type>(writeType)
            );
        }

        if (cstrIter == wordConstructorTablePtr_->end())
        {
            FatalErrorIn
            (
                "ibSurfaceWriter::New(const word&)"
            )   << "Unknown write type \"" << writeType << "\"\n\n"
                << "Valid write types : "
                << wordConstructorTablePtr_->sortedToc() << nl
                << "Valid proxy types : "
                << MeshedSurfaceProxy<face>::writeTypes() << endl
                << exit(FatalError);
        }
    }

    return autoPtr< ibSurfaceWriter<Type> >(cstrIter()());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::ibSurfaceWriter<Type>::ibSurfaceWriter()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::ibSurfaceWriter<Type>::~ibSurfaceWriter()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// ************************************************************************* //
