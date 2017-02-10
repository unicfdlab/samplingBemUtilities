/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/


#include "soundPressureSampler.H"
#include "volFields.H"
#include "dictionary.H"
#include "Time.H"
#include "argList.H"

//#include "wordReList.H"
//#include "PtrList.H"
//#include "PtrListIO.C"
//#include "RASModel.H"
//#include "LESModel.H"
//#include "basicThermo.H"

//sampledSurfaces stuff
#include "IOmanip.H"
#include "volPointInterpolation.H"
#include "PatchTools.H"
#include "ListListOps.H"
#include "stringListOps.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(soundPressureSampler, 0);
}

Foam::scalar Foam::soundPressureSampler::mergeTol_ = 1e-10;

// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * //


void Foam::soundPressureSampler::writeGeometry() const
{
    // Write to time directory under outputPath_
    // Skip surface without faces (eg, a failed cut-plane)


    const fileName outputDir = "surfaceGeometryData";

    Info << "write geometry \n";

    forAll(controlSurfaces_, surfI)
    {
        //Info << "in surfaces list\n";        
        
        const sampledSurface& s = controlSurfaces_.operator[](surfI);

        Info << "surfI size = " << mergeList_[surfI].faces.size() << nl;

        if (Pstream::parRun())
        {
            if (Pstream::master() && mergeList_[surfI].faces.size())
            {
                formatter_->write
                (
                    outputDir,
                    s.name(),
                    mergeList_[surfI].points,
                    mergeList_[surfI].faces
                );
            }
        }
        else if (s.faces().size())
        {
            formatter_->write
            (
                outputDir,
                s.name(),
                s.points(),
                s.faces()
            );
        }
    }
}



// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::soundPressureSampler::soundPressureSampler
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    log_(false),
    interpolationScheme_(word::null),
    controlSurfaces_(0),
    outputFormat_(word::null),
    pName_(word::null),
    mergeList_(),
    probeI_(0)
{
    read(dict);
}

// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::soundPressureSampler::~soundPressureSampler()
{}

// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

void Foam::soundPressureSampler::read(const dictionary& dict)
{
    //where should be log info 
    log_ = dict.lookupOrDefault<Switch>("log", false);
    
    if (!log_)
    {
        Info << "Direct logging to stdio disabled" << endl
        << " to enable, please insert string:" << endl
        << "log\t\t true;" << endl
        << "in dictionary" << endl;
    }

    dict.lookup("interpolationScheme") >> interpolationScheme_;

    const word writeType (dict.lookup("outputGeometryFormat"));

    //read mesh
    const fvMesh& mesh_ = refCast<const fvMesh>(obr_);

        // Define the surface formatter
        // Optionally defined extra controls for the output formats
        formatter_ = surfaceWriter::New
        (
            writeType,
            dict.subOrEmptyDict("formatOptions").subOrEmptyDict(writeType)
        );



    //read surfaces for sampling
    PtrList<sampledSurface> newList
    (
        dict.lookup("surfaces"),
        sampledSurface::iNew(mesh_)
    );

    controlSurfaces_.transfer(newList);

    //Parallel fix as it was implemented in sampledSurfaces class
    if (Pstream::parRun())
    {
        mergeList_.setSize(controlSurfaces_.size());
    }
        
    // Ensure all surfaces and merge information are expired
    expire();

    if (controlSurfaces_.size())
    {
        Info<< "Function object "<< name_<<":" << nl;

        Info<< " Reading control surface description:" << nl;
            
        forAll(controlSurfaces_, surfI)
        {
            Info<< " " << controlSurfaces_.operator[](surfI).name() << nl;
            Info<< "needsUpdate: " << controlSurfaces_.operator[](surfI).needsUpdate() << nl;
            Info<< "size: " << controlSurfaces_.operator[](surfI).points().size() << nl;
        }

        Info<< endl;
    }

    if (Pstream::master() && debug)
    {
        Pout<< "soundPressureSampler control surfaces additional info:" << nl << "(" << nl;
            
        forAll(controlSurfaces_, surfI)
        {
            Pout<< " " << controlSurfaces_.operator[](surfI) << endl;
        }
            
        Pout<< ")" << endl;
    }
    
    dict.lookup("pName") >> pName_;
}

void Foam::soundPressureSampler::execute()
{

}

bool Foam::soundPressureSampler::expire()
{
    bool justExpired = false;
    
    forAll(controlSurfaces_, surfI)                                                                                                                                                                                                                            
    {
        if (controlSurfaces_.operator[](surfI).expire())
        {
            justExpired = true;
        }
        //Clear merge information
        if (Pstream::parRun())
        {
            mergeList_[surfI].clear();
        }
    }
    
    // true if any surfaces just expired
    return justExpired;
}

bool Foam::soundPressureSampler::needsUpdate() const
{
    forAll(controlSurfaces_, surfI)
    {
        if (controlSurfaces_.operator[](surfI).needsUpdate())
        {
            return true;
        }
    }

    return false;
}

bool Foam::soundPressureSampler::update()
{
    //Actually the update() function copy-pasted from libSampling
    bool updated = false;
    
    if (!needsUpdate())
    {
        return updated;
    }
    
    // Serial: quick and easy, no merging required
    // Just like sampledSurfaces
    if (!Pstream::parRun())
    {
        forAll(controlSurfaces_, surfI)
        {
            if (controlSurfaces_.operator[](surfI).update())
            {
                updated = true;
            }
        }

        return updated;
    }

    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    // Dimension as fraction of mesh bounding box
    scalar mergeDim = mergeTol_ * mesh.bounds().mag();

    if (Pstream::master() && debug)
    {
        Pout<< nl << "Merging all points within "
        << mergeDim << " metre" << endl;
    }
    
    forAll(controlSurfaces_, surfI)
    {
        sampledSurface& s = controlSurfaces_.operator[](surfI);
        
        if (s.update())
        {
            updated = true;
        }
        else
        {
            continue;
        }
        
        PatchTools::gatherAndMerge
        (
            mergeDim,
            primitivePatch
            (
                SubList<face>(s.faces(), s.faces().size()),
                s.points()
            ),
            mergeList_[surfI].points,
            mergeList_[surfI].faces,
            mergeList_[surfI].pointsMap
        );
    }

    return updated;
}

void Foam::soundPressureSampler::end()
{
    Pout << "end function" << nl;
    // Do nothing - only valid on execute
}

void Foam::soundPressureSampler::timeSet()
{
// Do nothing - only valid on write
}

void Foam::soundPressureSampler::write()
{
    probeI_++;

    update();

    if (probeI_ == 1)
    {
        writeGeometry();
    }


    scalar cTime = obr_.time().value();
    
    fileName soundPressureSamplerDir;
    
    if (Pstream::master() && Pstream::parRun())
    {
        soundPressureSamplerDir = obr_.time().rootPath() + "/" + obr_.time().caseName().path() + "/pressureData/" + obr_.time().timeName() ;

        mkDir(soundPressureSamplerDir);
    }
    else if (!Pstream::parRun())
    {
        soundPressureSamplerDir = obr_.time().rootPath() + "/" + obr_.time().caseName() + "/pressureData/" + obr_.time().timeName();

        mkDir(soundPressureSamplerDir);
    }
    else
    {
    }


    //Write time point 
    autoPtr <OFstream> timeWriter;

    if (timeWriter.empty())
    {
        timeWriter.reset
        (
            new OFstream 
            (
                soundPressureSamplerDir + "/" + ("time.dat")
            )
        );
    };

    timeWriter() << cTime;


    forAll(controlSurfaces_, surfI)
    {
        const volScalarField& p = obr_.lookupObject<volScalarField>(pName_);

        sampledSurface& s = controlSurfaces_.operator[](surfI);

        autoPtr<OFstream> spsFilePtr;
        spsFilePtr.reset
        (
            new OFstream
            (
                soundPressureSamplerDir + "/" + (s.name() + ".txt")
            )
        );

        sampleAndWrite(p,surfI,spsFilePtr);
    }    
}

Foam::scalar Foam::soundPressureSampler::mergeTol()
{
    return mergeTol_;
}

Foam::scalar Foam::soundPressureSampler::mergeTol(const scalar tol)
{
    scalar oldTol = mergeTol_;
    mergeTol_ = tol;
    return oldTol;
}
// ************************************************************************* //
