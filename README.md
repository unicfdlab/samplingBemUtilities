# samplingBemUtilities

Utilities for sound pressure sampling in OpenFOAM, its fft transform and writing it in .gmsh format.
Using in hybrid CFD/BEM model.

History of changes:

v.1.2
-----

A bit optimisation of soundPressureSampler. Remove unnecessary data translating between functions.


v.1.1
-----

Fix bug with wrong numeration of writing data in soundPressureSampler.
Surface geometry writes automatically in first sample step.

surfaceNoise utility get a possibility to write spectrum of sampled data to .csv-file.
This data can be used in gnuplot to create a plot of spectrum. 


v.1.0
-----

1 function object (soundPressureSampler), 1 utility (surfaceNoise), 1 new writer (.gmsh format).

soundPressureSampler just sample pressure on the given control surface and writes history of data in large file. It writes only pressure data during OpenFOAM simulation.

After simulation:
1) run standard sample command without any data sampling in one time point to get the geometry data of control surface in .gmsh-format;
2) run surfaceNoise utility to make FFT of pressure data. Result is the sequence of files matching to different frequencies.


