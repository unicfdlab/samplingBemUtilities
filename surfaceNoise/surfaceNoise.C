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

Application
    surfaceNoise

Description
    Utility to perform noise analysis of pressure data in points of surface using the libFFT
    library.

    Control settings are read from the $FOAM_CASE/system/surfaceNoiseDict dictionary,
    or user-specified dictionary using the -dict option. Pressure data prepared by
    surfacePressureSampling utility and read using a simple file reader:

    \heading Usage

    \verbatim

    outputFormat gmsh;

    pressureData
    {
        fileName        "pressureData.dat"
    }
    \endverbatim


    Output:
    - FFT of the pressure data

SeeAlso
    FoamFftwDriver.H

\*---------------------------------------------------------------------------*/


#include "FoamFourierAnalysis/FoamFftwDriver.H"
#include "argList.H"
#include "Time.H"
#include "functionObjectFile.H"
#include "complex.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::List<Foam::List<scalar > > > read (const fileName& pFileName)
{
    fileName fullFileName = "pressureData/" + pFileName;

    autoPtr<IFstream> inFilePtr;

    if (inFilePtr.empty())
    {
    // Open new file at start up
        inFilePtr.reset
        (
            new IFstream
            (
                fullFileName
            )
        );
    }

    List<List<scalar> > data;

    label probeI = 0;

    scalar num = 0;

    string line;

    if (inFilePtr.valid())
    {
        //Info << inFilePtr().name() << endl;

        while(inFilePtr().getLine(line))
        {
            probeI ++;

            data.resize(probeI);

            //Info << "resize ok" << endl;

            label pointI = 0;

            IStringStream lineStream(line);

            while(lineStream >> num)
            {  
                pointI++;

                data[probeI-1].resize(pointI);

                data[probeI-1][pointI-1] = num;
            }

            //inFilePtr() >> num ;
            //Info << num << nl;

            //Info << "///////////////////////////////" << nl;
            //istr++;
        }

        //check if last string has equal size with previous strings
        label refSize = data[0].size();
        label lastStrSize = data[probeI-1].size();

        //Info << "ref size = " << refSize << ", size of last string = " << lastStrSize << endl;

        if (lastStrSize < refSize) // in case of short lasr string cut it because of lost part of data
        {
            Info << "Cut tail of damaged data"  << endl;
            data.resize(probeI-1);
        }

        Info << "OK" << nl;
    }
    else
    {
        Info << " No valid file" << endl;
    }

    return autoPtr<List<List<scalar> > >
    (
        new List<List<scalar> >
        (
            data
        )
    );

}

// just return Re and Im parts of complex pressure amplitude; frequency is out of function
Foam::autoPtr<Foam::List<Foam::complex > > fft(const List<scalar> &p, scalar lastTime)
{
    List<complex > fft_res;
    
    fft_res.resize(0);

    if ( (p.size() > 0) ) 
    {
        FoamFftwDriver fftw (p, lastTime);
        
        //Info << "lastTime = " << lastTime << endl;

        autoPtr<List<List<scalar> > > pfft = fftw.simpleForwardTransform();
        
       //Info << pfft() << endl;

        fft_res.resize(pfft().first().size());

        forAll (pfft().first(), k)
        {
            fft_res[k] = complex(pfft()[0][k],pfft()[1][k]);

            //fft_res[0][k] = pfft()[0][k]; // Re of pressure amplitude, Pa
            //fft_res[1][k] = pfft()[1][k]; // Im of pressure amplitude, Pa
        }
    }
    
    return autoPtr<List<complex > >
    (
        new List<complex >
        (
            fft_res
        )
    );
}

void writeComplexNumber(const complex& number, autoPtr<OFstream>& os)
{
    os()    << number.Re();

    label sgn = sign(number.Im());

    if (sgn == 1)
    {
        os() << "+";
    }
    else
    {
        os() << "-";
    }
    
    os() << fabs(number.Im());
    os() << "j";
}

void writeToGmsh (const fileName& outputFile, const List<complex >& data)
{

    if (Pstream::master() || !Pstream::parRun())
    {
        autoPtr<OFstream> filePtr;

        if (filePtr.empty())
        {
        // Open new file at start up
            filePtr.reset
            (
                new OFstream
                (
                    outputFile
                )
            );
        }


        filePtr()   << "$NodeData" << nl;

        //string tags
        filePtr()   << "1" << nl;   // how much
        filePtr()   << "node_data" << nl;  //name of data

        //real tags
        filePtr()   << "1" << nl    // how much
                    << "0" << nl;   // should be output time (just 0 for all files)
            
        //integer tags
        filePtr()   << "4" << nl    // how much
                    << "0" << nl    // time step index
                    << "1" << nl            // how many field components
                    << data.size() << nl  // number of entities (nodes or elements)
                    << "0" << nl;

        forAll(data, probeI)
        {
            filePtr() << probeI+1 << ' ';
            writeComplexNumber(data[probeI],filePtr);
            filePtr() << nl;
                    
        }

        filePtr() << "$EndNodeData" << nl;
    }

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();

    #include "addDictOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createFields.H"
    

    Info<< "Reading data file..." << endl;

    autoPtr<List<List<scalar> > > ptData = read(pFileName);

    //total time
    label timePoints = ptData().size();

    // Info << "timePoints = " << timePoints << endl;

    scalar lastTime = ptData()[timePoints-1][0];
        
    // Info << "lastStr = " << ptData()[0].size() << ' ' << ptData()[timePoints-1].size() << endl;

    // pressure data
    List<scalar> p; // just buffer
    p.resize(timePoints);


    //------------------------


    // list of frequencies

    scalar freqPoints = ceil( scalar(timePoints)*0.5 );
    List<scalar> freq(freqPoints);

    forAll(freq,i)
    {
        freq[i] = scalar(i)/lastTime;
    }


    //fft for all points on surface point-by-point

    autoPtr<List<complex > > fftOnePoint;

    List<List<complex > > allComplexPressures;
    allComplexPressures.resize(freqPoints);

    label points = ptData()[0].size();

    Info<< "Calculate FFT..." << endl;

    for (int i = 1; i < points; ++i) //loop starts from 1 because time value is placed on 0 position
    {            
        
        forAll(p,pI)
        {
            p[pI] = ptData()[pI][i];
        }

        //Info << "Creating  FFT" << endl;
        fftOnePoint = fft(p, lastTime);

        forAll(allComplexPressures,freqI)
        {
            allComplexPressures[freqI].resize(i);
            allComplexPressures[freqI][i-1] = fftOnePoint()[freqI];
        }

        //Info << fftOnePoint() << endl;
    }


    //Info << "frequencies: " << freq << nl;

    //Info << "last fft = " << allComplexPressures << nl;


    //------------------------

    Info<< "Write gmsh files..." << endl;

    fileName controlFile = "complexPressureData/" + outputDirectory + "Spectrum.dat";

    autoPtr<OFstream> controlFilePtr;

    if (controlFilePtr.empty())
    {
        // Open new file at start up
        controlFilePtr.reset
        (
            new OFstream
            (
                controlFile
            )
        );
    }

    forAll(freq, freqI)
    {
        word freqValue = name(freq[freqI]);

        fileName freqDirectory = "complexPressureData/" + outputDirectory + "/" + freqValue;

        //Info << freqDirectory << endl;

        mkDir(freqDirectory);

        //write appropriate frequency
        fileName freqFile  = freqDirectory + "/frequency";

        autoPtr<OFstream> freqFilePtr;

        if (freqFilePtr.empty())
        {
        // Open new file at start up
            freqFilePtr.reset
            (
                new OFstream
                (
                    freqFile
                )
            );
        }

        freqFilePtr() << freq[freqI];

        //write data to msh file
        if (outputFormat == "gmsh")
        {
            fileName dataFile = freqDirectory + "/nodeData.gmsh";
            
            writeToGmsh(dataFile,allComplexPressures[freqI]);
        }
        
        Info << "---------------\n";
        Info << "Frequency: " << freq[freqI] << nl;
        
        scalar maxAbs = 0.0;
        
        forAll(allComplexPressures[freqI],k)
        {
    	    complex currentCA = allComplexPressures[freqI][k];
    	    scalar currentAbs = Foam::sqrt( pow(currentCA.Re(),2) + pow(currentCA.Im(),2) );
    	    
    	    if (maxAbs < currentAbs)
    	    {
    	        maxAbs = currentAbs;
    	    }
    	    
        }
        
        Info << "Max abs: " << maxAbs << nl;

Info << "write..." << nl;

        controlFilePtr() << freq[freqI] << ' ' << maxAbs << nl;
    }

    return 0;
}


// ************************************************************************* //
