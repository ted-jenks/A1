import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import GalaxyCounter as gc


#%%
'''
_________________________________________________________________________________________________

Testing the Aperture Function
_________________________________________________________________________________________________

'''
xcoords,ycoords,numpixels = gc.inAperture(7,10, 13)         #aperture test
plt.plot(xcoords,ycoords, 's', markersize = 15)        #aperture plot
plt.xlabel('x')
plt.ylabel('y')
plt.title('Aperture Plot')

#%%
'''
_________________________________________________________________________________________________

Testing the Source Detector
_________________________________________________________________________________________________

'''

#EMPTY ARRAY TEST
gc.EmptyArrayTest()     #0 objects expected

#ONE PIXEL TEST
gc.OnePixelTest()       #1 object expected

#GAUSSIAN BLOB TEST IMAGE
gc.GaussBlobTest()      #5 objects expected

#FIXED APERTURE --- IMAGE 1
bins, cumulative, popt, error, number, dataMasked = gc.RunSmallImageTest1(method = 'fixed',empty = 0.9, minpix = 1, aperture = 13)
print('Fixed:',number, 'objects detected.')
plt.figure()
dataRead = dataMasked.filled(10000)
plt.imshow(dataRead, cmap='gray',norm=LogNorm(vmin = 3455, vmax = 3500))                   #view the image
plt.gca().invert_yaxis()
plt.xlim(1847,2018)
plt.ylim(339,552)
plt.title('Fixed Aperture = 13 -- Detection Locations')
plt.colorbar()    

#VARIABLE APERTURE --- IMAGE 1
bins, cumulative, popt, error, number, dataMasked = gc.RunSmallImageTest1(method = 'variable',empty = 0.9, minpix = 40, aperture = 13)
print('Variable:',number, 'objects detected.')
plt.figure()
dataRead = dataMasked.filled(10000)
plt.imshow(dataRead, cmap='gray',norm=LogNorm(vmin = 3455, vmax = 3500))                   #view the image
plt.gca().invert_yaxis()
plt.xlim(1847,2018)
plt.ylim(339,552)
plt.title('Variable Aperture -- Detection Locations')
plt.colorbar()   

#FIXED APERTURE --- IMAGE 2
bins, cumulative, popt, error, number, dataMasked = gc.RunSmallImageTest2(method = 'fixed',empty = 0.9, minpix = 1, aperture = 13)
print('Fixed:',number, 'objects detected.')
plt.figure()
dataRead = dataMasked.filled(10000)
plt.imshow(dataRead, cmap='gray',norm=LogNorm(vmin = 3455, vmax = 3500))                   #view the image
plt.gca().invert_yaxis()
plt.xlim(709,975)
plt.ylim(1139,1261)
plt.title('Fixed Aperture = 13 -- Detection Locations')
plt.colorbar()  

#VARIABLE APERTURE --- IMAGE 2
bins, cumulative, popt, error, number, dataMasked = gc.RunSmallImageTest2(method = 'variable',empty = 0.9, minpix = 40, aperture = 13)
print('Variable:',number, 'objects detected.')
plt.figure()
dataRead = dataMasked.filled(10000)
plt.imshow(dataRead, cmap='gray',norm=LogNorm(vmin = 3455, vmax = 3500))                   #view the image
plt.gca().invert_yaxis()
plt.xlim(709,975)
plt.ylim(1139,1261)
plt.title('Variable Aperture -- Detection Locations')
plt.colorbar()    