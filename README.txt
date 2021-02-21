A1			         |
Astronomical Image Proccessing   |
Edward Jenks and Kelvin Brinham	 |
_________________________________|

This repository contains the code for our A1 lab experiment.
We have constructed and tested a galaxy counter in Python.
The code has been applied on an image from the SWIRE infrared survey taken at the Kitt Peak Observatory.


FILE BREAKDOWN
--------------

The project is sorted into 3 files

-GalaxyCounter.py
	Contains the core code of the program.
	All functions used to run the galaxy count can be found in this file.

-Run.py
	Module to run the galaxy counter and present graphs of the results.

-Validation.py
	Module containing the validation we carried out on our code to ensure that it worked as expected.


RUNNING THE CODE
----------------

To run the code open the Run.py module.
This has been set up with the configuration to gather our results and plot neccessary graphs.
The code can be edited to run different configurations.

The key function to gather results is Run(method, empty, minpix, aperture) where the keywords are:

method: Either 'fixed' or 'variable' to choose between a fixed or variable aperture.
empty: The fraction of cells that need to be empty to set the aperture (only relavent for variable aperture).
minpix: The minimum number of pixels for a source to be counted (only relevent for variable aperture).
aperture: The diameter of the fixed aperture (only relevent for fixed aperture).

Changine these variables will give different results. 
The configuration we ended up choosing proved the best in our testing.

When run the program will generate 'ObjectCatalogue.dat' and 'PlotData.dat'.

'ObjectCatalogue.dat' is the catalogue of objects detected, giving their position, pixel value, magnitude, local background and their size in pixels.
'PlotData.dat' is the objects sorted into bins for plotting, giving bin locations, the cumulative count up to that magnitude and the error on the count.


VALIDATING THE CODE
-------------------

Validation is carried out in Validation.py.

The first test ensures the aperture function works properly.
A graph should be plotted of the layout of pixels in a pixel approximation of a circle.

The second test gives the counter an empty array.
No objects should be detected.

The third gives the counter an array with a single pixel bright.
One object should be detected.

The fourth test hands the counter an array with 5 Gaussian blobs.
5 objects should be detected.

The remainder of the tests return the number of objects and detection locations on small sections of the image.
This highlights the strengths and weaknesses of the methods for edge cases such as irregular or close objects.


CORE STRUCTURE
--------------

{{THIS SECTION IS NOT NECCESSARY TO RUN THE CODE}}

1. Data Read
2. Data masked by 'maskData()'
3. 'counter()' ran
4. Start loop
5. 'findmax()' ran
6. Max unmasked pixel found
7. Total pixel count in aperture determined
8. Pixels masked
9. Background adjustment
10. Magnitude calculated and value saved
11. If cut-off not reached return to 5