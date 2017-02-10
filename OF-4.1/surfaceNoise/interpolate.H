	// get pressure values in one point (just copy from buffer)
        forAll(pSampled,pI)
        {
            pSampled[pI] = ptData()[pI][i];
        }

        // interpolate pressure data
        scalar cTimeInt = deltaIntT;
        label tI = 1;

        pInt[0] = pSampled[0];

        for (label j = 1; j < nFftTimePoints; ++j)
        {
            while (tSampled[tI] < cTimeInt)
            {
                tI++;
                continue;
            }

            pInt[j] = linearInterpolation(cTimeInt,tSampled[tI],pSampled[tI],tSampled[tI-1],pSampled[tI-1]);
            
            cTimeInt += deltaIntT;
        }