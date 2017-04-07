# samplingBemUtilities

Utilities for sound pressure sampling in OpenFOAM, its fft transform and writing it in .gmsh format.
Using in hybrid CFD/BEM model.

History of changes:

v.2.0
-----
1) support of OpenFOAM v.4.1;
2) refactoring of surfaceNoise utility;
3) small bug (moment of writitng of control surface geometry) fixed.


v.1.2
-----

Optimisation of soundPressureSampler. Unnecessary data translating between functions fixed.


v.1.1
-----

Fix bug with wrong numeration of writing data in soundPressureSampler.
Surface geometry is written automatically in first sample step.

surfaceNoise utility has got the possibility to write spectrum of sampled data (max value of amplitude on the surface in dependency of frequency) to .csv-file.
This data can be used in gnuplot to create a plot of spectrum. 


v.1.0
-----

1 function object (soundPressureSampler), 1 utility (surfaceNoise), 1 new writer (.gmsh format).

soundPressureSampler function object just sample pressure on the given control surface and writes history of pressure changing in large file during OpenFOAM simulation.

After simulation:
1) run standard sample command without any data sampling in one time point to get the geometry data of control surface in .gmsh-format;
2) run surfaceNoise utility to make FFT of pressure data. Result is the sequence of files matching to different frequencies.


