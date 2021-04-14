/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
    Copyright (C) YEAR AUTHOR, AFFILIATION
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

\*---------------------------------------------------------------------------*/

#include "codedFvOptionTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
#include "fvMatrix.H"

//{{{ begin codeInclude

//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// dynamicCode:
// SHA1 = 4be0633f558cce96a8701b85600366d18eee499a
//
// unique function name that can be checked if the correct library version
// has been loaded
extern "C" void headache_4be0633f558cce96a8701b85600366d18eee499a(bool load)
{
    if (load)
    {
        // Code that can be explicitly executed after loading
    }
    else
    {
        // Code that can be explicitly executed before unloading
    }
}


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(headacheFvOptionscalarSource, 0);
addRemovableToRunTimeSelectionTable
(
    option,
    headacheFvOptionscalarSource,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

headacheFvOptionscalarSource::
headacheFvOptionscalarSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(name, modelType, dict, mesh)
{
    if (false)
    {
        printMessage("Construct headache from components");
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

headacheFvOptionscalarSource::
~headacheFvOptionscalarSource()
{
    if (false)
    {
        printMessage("Destroy headache");
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void headacheFvOptionscalarSource::correct
(
    GeometricField<scalar, fvPatchField, volMesh>& fld
)
{
    if (false)
    {
        Info<<"headacheFvOptionscalarSource::correct()\n";
    }

//{{{ begin code
    #line 70 "/home/zstyle/OpenFOAM_files/hotheater1/constant/ellipse/fvOptions.energySource.scalarCodedSourceCoeffs"
Pout<< "**codeCorrect**" << endl;
//}}} end code
}


void headacheFvOptionscalarSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (false)
    {
        Info<<"headacheFvOptionscalarSource::addSup()\n";
    }

//{{{ begin code
    #line 75 "/home/zstyle/OpenFOAM_files/hotheater1/constant/ellipse/fvOptions.energySource.scalarCodedSourceCoeffs"
const Time& time = mesh().time();
		//const scalarField& V = mesh_.V();
                const scalar  F_0 = 82; 			//Laser to heat-Jm^-2
	 	const scalar  dt = 5e-9;			//Half-width-s
	 	const scalar  t_0 = 30e-9;		        //Peak-s
	 	const scalar  C = 4759e-18;      		//Laser absorption-m^3
	 	const scalar  volu = 2.932153143e-26;           //Volume m^3
		const scalar  F_t =(F_0/(dt*sqrt(2*3.1415)))*exp(-pow(time.value() - t_0,2)/(2*pow(dt,2)));
                scalarField& hSource = eqn.source();
		hSource = (C*F_t);
//}}} end code
}


void headacheFvOptionscalarSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (false)
    {
        Info<<"headacheFvOptionscalarSource::addSup()\n";
    }

//{{{ begin code
    #line 75 "/home/zstyle/OpenFOAM_files/hotheater1/constant/ellipse/fvOptions.energySource.scalarCodedSourceCoeffs"
const Time& time = mesh().time();
		//const scalarField& V = mesh_.V();
                const scalar  F_0 = 82; 			//Laser to heat-Jm^-2
	 	const scalar  dt = 5e-9;			//Half-width-s
	 	const scalar  t_0 = 30e-9;		        //Peak-s
	 	const scalar  C = 4759e-18;      		//Laser absorption-m^3
	 	const scalar  volu = 2.932153143e-26;           //Volume m^3
		const scalar  F_t =(F_0/(dt*sqrt(2*3.1415)))*exp(-pow(time.value() - t_0,2)/(2*pow(dt,2)));
                scalarField& hSource = eqn.source();
		hSource = (C*F_t);
//}}} end code
}


void headacheFvOptionscalarSource::constrain
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (false)
    {
        Info<<"headacheFvOptionscalarSource::constrain()\n";
    }

//{{{ begin code
    #line 89 "/home/zstyle/OpenFOAM_files/hotheater1/constant/ellipse/fvOptions.energySource.scalarCodedSourceCoeffs"
Pout<< "**codeConstrain**" << endl;
//}}} end code
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// ************************************************************************* //

