    word dictName("surfaceNoiseDict");
    if (args.optionFound("dict"))
    {
        dictName = args["dict"];
    }

    IOdictionary dict
    (
        IOobject
        (
            dictName,
            runTime.system(),
            runTime,
            IOobject::MUST_READ
        )
    );

    // number of samples in sampling window
//    label N = dict.lookupOrDefault("N", 65536);

    // number of sampling windows
//    label nw = dict.lookupOrDefault("nw", 100);

    // lower frequency of frequency band
//    scalar f1 = dict.lookupOrDefault("f1", 25.0);

    // upper frequency of frequency band
//    scalar fU = dict.lookupOrDefault("fU", 10000.0);


    // output format
    word outputFormat = dict.lookupOrDefault<word>("outputFormat", "gmsh");

    // output directory (similar with surface name)
    fileName outputDirectory = dict.lookup("surfaceName");

    // file name with pressure history --- will construct automatically by names of sampled surfaces
    fileName pFileName = dict.subDict("pressureData").lookup("fileName");
