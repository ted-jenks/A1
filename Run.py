from astropy.io import fits
from astropy.io import ascii
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.optimize import curve_fit
import numpy as np
import GalaxyCounter as gc



#%%
'''
_________________________________________________________________________________________________

Read File
_________________________________________________________________________________________________
'''

hdulist = fits.open("mosaic.fits")      #Open the image file

header = hdulist[0].header              #Open the image file header
zeropoint = hdulist[0].header["MAGZPT"] #Read the zero point value 
data = hdulist[0].data                  #Read the count data 
TotalAngle = 0.04509811                 #The solid angle (in deg^2) subtented by the data we used (ie solid angle of image - solid angle of masking)

#%%
'''
_________________________________________________________________________________________________

Background Data
_________________________________________________________________________________________________
'''

print('BACKGROUND DATA')
print('_______________')
print('\n')
flatdata = data.flatten()
snippeddata = []
for i in range(len(flatdata)):
    if flatdata[i] >= 3300 and flatdata[i] <= 3600:
        snippeddata.append(flatdata[i])

plt.figure()
plt.grid()
BINS = 100                                        #no. bins for hist
counts, bins, bars = plt.hist(snippeddata, BINS, label = 'Pixel Vales')       #flatten array and plot hist
plt.xlim(3300,3600)
plt.title('Histogram of Pixel Values')
plt.xlabel('Pixel Value')
plt.ylabel('Frequancy')

def gauss(x, mu, sigma,c):
    '''
    Gaussian Function.

    Parameters
    ----------
    x : flt
        x coord.
    mu : flt
        mean.
    sigma : flt
        std.

    Returns
    -------
    flt
        Gauss(x).

    '''
    return c*(1/(sigma*np.sqrt(2*np.pi))) * np.exp(-0.5*((x-mu)/sigma)**2)

popt, pcov = curve_fit(gauss, bins[:100]+(bins[1]-bins[0])/2,counts, maxfev=10000, p0 = np.array([3419,12,3.97e07]))
error = np.sqrt(np.diag(pcov))
mu = popt[0]
sigma = popt[1]
x = np.linspace(3300,3600,1000)
plt.plot(x, gauss(x,popt[0], popt[1], popt[2]), label = 'Gaussian Fit')        #plot gaussian function over hist
plt.legend()
plt.text(3475,0.8e6,'Mean = 3419.2 +/- 0.2',bbox={'facecolor': 'white', 'alpha': 1, 'pad': 10})
plt.text(3475,0.6e6,'Std = 12.2 +/- 0.2',bbox={'facecolor': 'white', 'alpha': 1, 'pad': 10})

print('Gaussian Mean:', mu, '+/-', error[0])
print('Gaussian Stdev:', sigma, '+/-', error[1])

#%%
'''
_________________________________________________________________________________________________

Initial Image Plots
_________________________________________________________________________________________________
'''

plt.figure()
plt.imshow(data, cmap='gray')                   #view the image
plt.title('Image')
plt.colorbar()                                  #add colour bar
plt.gca().invert_yaxis()

plt.figure()
plt.imshow(data, cmap='gray', norm=LogNorm(vmin = 3000, vmax = 4000))  #view image with logarithmic brightness scaling
plt.title('Image with Logarithmic Brightness Values')
plt.colorbar()
plt.gca().invert_yaxis()

#%%
'''
_________________________________________________________________________________________________

Data Collection
_________________________________________________________________________________________________
'''

def fit(x,a,c):
    '''
    Fit function.

    Parameters
    ----------
    x : flt
        x value.
    a : flt
        gradient.
    c : flt
        constant.

    Returns
    -------
    flt
        value of function.

    '''
    return 10**(a*x+c)

def fitStraight(x,m,c):
    '''
    Straight line fit.

    Parameters
    ----------
    x : flt
        x value.
    m : flt
        gradiant.
    c : flt
        constant.

    Returns
    -------
    flt
        function value.

    '''
    return -m*x+c

def Run(method = 'variable', empty = 0.9, minpix = 2, aperture = 13):
    '''
    Runs count on the image.

    Parameters
    ----------
    method : str
        variable or fixed.
    empty : flt
        fraction of ring empty for variabe method.
    minpix : int
        minimum object size in pixels.
    aperture : int
        aperture diameter in pixels..

    Returns
    -------
    bins : arr
        bin locations of data.
    cumulative : arr
        cumulative count.
    popt : array
        fit data.
    error : arr
        error data.
    number : int
        number of objects.
    '''
    global TotalAngle
    
    gc.maskData()
    
    Vals, Mags, Xcoords, Ycoords, localBackgrounds, numpixels = gc.counter(method, empty, minpix, aperture)

    MagPlot = []
    for i in range(len(Mags)):
        if Mags[i]>0 and Mags[i]<=17:
            MagPlot.append(Mags[i])                       #Adds magnitude data to MagPlot to be plotted
            
    BINS = 30                                             #Number of histogram bins
    values, bins = np.histogram(Mags, BINS)               #Sorts magnitudes into histogram bins
    cumulative = np.cumsum(values)                        #Cumulative count of number of objects with magnitude
    
    BINSfit = 25                                          #Number of histogram bins
    valuesfit, binsfit = np.histogram(MagPlot, BINSfit)
    cumulativefit = np.cumsum(valuesfit) 
    cumulativefit = cumulativefit[8:]
    cumulativefit = cumulativefit/TotalAngle
    binsfit = binsfit[8:]
    popt, pcov = curve_fit(fit, binsfit[:-1]+(binsfit[1]-binsfit[0])/2,cumulativefit, sigma = np.sqrt(cumulativefit), maxfev=10000, p0 = np.array([0.6, 2]))
    error = np.sqrt(np.diag(pcov))
    
    number = len(Mags)
    
    ascii.write([Xcoords, Ycoords,Vals,Mags,localBackgrounds,numpixels], 'ObjectCatalogue.dat', names=['x', 'y','Pixel Value','Magnitude','Local Background', 'Size of Source (Total Pixel Count)'],delimiter = '|', overwrite=True, format='fixed_width')
    ascii.write([bins[:-1]+(bins[1]-bins[0])/2,cumulative,np.sqrt(cumulative)], 'PlotData.dat', names=['Magnitude','Cumulative Count','Error'], overwrite=True, format='fixed_width')
    return (bins, cumulative, popt, error, number, Mags,numpixels)

bins, cumulative, popt, error, number, Mags, numpixels = Run('variable',  0.9, 40, 13)


#%%
'''
_________________________________________________________________________________________________

Data Plotting
_________________________________________________________________________________________________
'''

#Variable Aperture Galaxy Count
plt.figure()
x = np.linspace(7,26,1000)
plt.plot(x, fit(x,popt[0], popt[1]),'--', label = 'Fit',color = 'grey')
plt.errorbar(bins[:-1]+(bins[1]-bins[0])/2,cumulative/TotalAngle, yerr = np.sqrt(cumulative/TotalAngle), fmt = '.', capsize = 2, label = 'Data', color = 'black')
plt.ylim(5e1,1e5)
plt.xlim(10,20)
plt.grid()
plt.title('Variable Aperture Galaxy Count')
plt.xlabel('Magnitude')
plt.ylabel('N(m)')
plt.yscale('log',nonposy = 'clip')
plt.legend()
plt.text(15,9e2,'Gradient = 0.303 +/- 0.012',bbox={'facecolor': 'white', 'alpha': 1, 'pad': 10})
plt.text(15,3e2,'Constant = -0.55 +/- 0.13',bbox={'facecolor': 'white', 'alpha': 1, 'pad': 10})
plt.text(15,1e2,'2657 objects detected.',bbox={'facecolor': 'white', 'alpha': 1, 'pad': 10})
print('Variable Aperture Galaxy Count')
print('______________________________')
print('\n\n')
print('Gradient:',popt[0],'+/-',error[0])
print('Constant:',popt[1],'+/-',error[1])
print(number, 'objects detected')


#Size of Source Vs Magnitude
plt.figure()
plt.plot(Mags, np.log10(numpixels),'x', label = 'Data',color = 'grey')
popt, pcov = curve_fit(fitStraight, Mags, np.log10(numpixels), maxfev=10000)
error = np.sqrt(np.diag(pcov))
x = np.linspace(7,26,1000)
plt.ylim(1,5)
plt.xlim(8,22)
plt.plot(x, fitStraight(x,popt[0], popt[1]),'--', label = 'Fit',color = 'grey')
plt.title('Size of Source Vs Magnitude')
plt.xlabel('Magnitude')
plt.ylabel('Log(Number of Pixels)')
plt.text(9,2,'m = -0.195 +/- 0.002',bbox={'facecolor': 'white', 'alpha': 1, 'pad': 10})
plt.text(9,1.5,'C = 5.31 +/- 0.04',bbox={'facecolor': 'white', 'alpha': 1, 'pad': 10})
plt.grid()
plt.legend()
print('\n\n')
print('Object Size Vs Magnitude')
print('________________________')
print('\n\n')
print('Gradient:',popt[0],'+/-',error[0])
print('Constant:',popt[1],'+/-',error[1])
