    fileName rootDirectory = "pressureData/";
    fileName pFileName = inputFileName;

    autoPtr<List<List<scalar> > > ptData = read(pFileName);

    Info<< "OK" << endl;

    //number of sampled times
    label sampledTimePoints = ptData().size();

    List<scalar> tSampled;
    tSampled.resize(sampledTimePoints);


    forAll(tSampled,tI)
    {
        tSampled[tI] = ptData()[tI][0];
    }

    scalar lastTime = tSampled[sampledTimePoints-1];

    // pressure data
    List<scalar> pSampled; // just buffer for pressure data in one point
    pSampled.resize(sampledTimePoints);
