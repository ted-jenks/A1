from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np



#%%
'''
_________________________________________________________________________________________________

Read File
_________________________________________________________________________________________________
'''

hdulist = fits.open("mosaic.fits")          #Open the image file

header = hdulist[0].header                  #Open the image file header
zeropoint = hdulist[0].header["MAGZPT"]     #Read the zero point value 
data = hdulist[0].data                      #Read the count data 
dataMasked = np.array([])                   #Create an array for the masked data
dataRead = np.array([])                     #Create an array for the masked data filled with background


def dataReset():
    '''
    Resets Data

    Returns
    -------
    None.

    '''
    
    global hdulist
    global header
    global zeropoint
    global data
    global dataMasked
    global dataRead
    
    hdulist = fits.open("mosaic.fits")

    header = hdulist[0].header
    zeropoint = hdulist[0].header["MAGZPT"]
    data = hdulist[0].data
    dataMasked = np.array([])
    dataRead = np.array([])

    
    return
#%%
'''
_________________________________________________________________________________________________

Background Data
_________________________________________________________________________________________________
'''

mu = 3419.1978290249785         #these values are from the gaussian
sigma = 12.231015056830977      #they have been hardcoded to avoid having to run it every time

#%%
'''
_________________________________________________________________________________________________

Masking
_________________________________________________________________________________________________
'''

def drawRect(arr,x, tl, br):
    '''
    Function to draw a rectangle of values on array.

    Parameters
    ----------
    x : flt
        Value to fill cells with.
    tl : arr
        Coords of top left corner.
    br : arr
        Coords of bottom right corner.

    Returns
    -------
    None.

    '''
    global counter
    for i in range(tl[0],br[0]):            #Iterates over x coordinate values within rectangle 
        for j in range(tl[1],br[1]):        #Iterates over y coordinate values qithin rectangle CHECK X AND Y RIGHT WAY ROUND HERE!!!
            arr[i][j] = x                   #Sets the points within the rectangle to the value x (not coordinate x)
    return

def circleMask(arr, x, i,j,apertureDiameter):

    '''
    Function to return the pixels within the circular aperture centered on i,j.
 

    Parameters
    ----------
    arr : arr
        array to edit.
    x : flt
        Value to fill cells with.
    i : int
        x coordinate of central point.
    j : int
        y coordinate of central point.
    apertureDiameter : int
        Aperture diameter.

    Returns
    -------

    NONE.
    '''
    currentx = i - (apertureDiameter-1)/2 - 1           #Sets starting x coordinate to left side of the square surrounding the circle
    currenty = j - (apertureDiameter-1)/2 - 1           #Sets starting y coordinate to the bottom side of the square surrounding the circle
                                                        
    for k in range(apertureDiameter):
        currentx += 1                                   #Iterates through the x coordinates
        currenty = j - (apertureDiameter-1)/2 - 1       #WHY IS THIS REPEATED!!!
        for h in range(apertureDiameter):
            currenty += 1                               #Iterates through the y coordinates
            if round(np.sqrt((currentx-i)**2+(currenty-j)**2)) <= (apertureDiameter-1)/2: #Checks if point is within the circle
                arr[int(currentx)][int(currenty)] = x   #Sets the values within the circlular aperture to x (not coordinate x)
                
    return

def maskData():
    '''
    Mask the data set.
 
    Returns
    -------
    None.
    '''

    global dataMasked
    global dataRead
   
    mask = np.zeros([int(np.shape(data)[0]),int(np.shape(data)[1])])   #Creates an array to hold the masked pixels (array of 0s)
   
    drawRect(mask, 1, [0,0], [150,2570])            #Crops the image edges
    drawRect(mask, 1, [0,0], [4611,150])
    drawRect(mask, 1, [0,2419], [4611,2570])
    drawRect(mask, 1, [4460,0], [4611,2570])
    
    circleMask(mask, 1, 4095, 558, 90)              #Masks circular objects (foreground stars)                 
    circleMask(mask, 1, 3315, 785, 271)
    circleMask(mask, 1, 2775, 975, 191)
    circleMask(mask, 1, 2290, 901, 181)
    circleMask(mask, 1, 2265, 703, 71)
    circleMask(mask, 1, 2304, 448, 69)

    circleMask(mask, 1, 1940,677, 61)
    circleMask(mask, 1, 1655, 978, 51)
    circleMask(mask, 1, 1498, 636, 65)
    circleMask(mask, 1, 848, 211, 65)
    circleMask(mask, 1, 667, 68, 127)
    circleMask(mask, 1, 4331, 1367, 65)
    circleMask(mask, 1, 4399, 1310, 61)
                         

    drawRect(mask, 1, [403,1000], [521,1726])       #Masks bleeding lines from bright objects
    drawRect(mask, 1, [475,993], [538,1046]) 
    drawRect(mask, 1, [295,988], [363,1749])
    drawRect(mask, 1, [209,1380], [274,1501]) 
    drawRect(mask, 1, [111,1278], [165,1557]) 

    circleMask(mask, 1, 3930, 185, 63)
    circleMask(mask, 1, 4284, 65, 65)

    drawRect(mask, 1, [0,1377], [4610,1480]) 

    circleMask(mask, 1, 4031, 1458, 77)
    circleMask(mask, 1, 3206, 1422, 575)
    circleMask(mask, 1, 3854, 2279, 65)
    circleMask(mask, 1,3760, 2095, 191)
    circleMask(mask, 1, 3420, 2462, 113)
    circleMask(mask, 1, 2206, 2277, 99)
    circleMask(mask, 1, 2987, 1415, 59)
    circleMask(mask, 1, 3035, 2519, 31)
    circleMask(mask, 1, 2964, 1495, 37)
    circleMask(mask, 1, 2313, 2132, 107)
    circleMask(mask, 1, 2306, 2334, 139)
    circleMask(mask, 1, 1960, 2490, 101)
    circleMask(mask, 1, 1433, 2096, 137)
    circleMask(mask, 1, 3297, 2260, 160)

    circleMask(mask, 1, 1125, 2500, 99)
    circleMask(mask, 1, 581, 1771, 43)
    circleMask(mask, 1, 467, 1601, 59)
    circleMask(mask, 1, 291, 1796, 65)   

    dataMasked = np.ma.masked_array(data,mask)      #Holds an array of the data with a second array indicating the masked points

    dataRead = dataMasked.filled(mu)                #Creates an array of the data with masked points containing mu
    
    return

#%%
'''
_________________________________________________________________________________________________

Source Detection
_________________________________________________________________________________________________
'''

def inAperture(i,j,apertureDiameter):
    '''
    Function to return the pixels within the circular aperture centered on i,j.

    Parameters
    ----------
    i : int
        x coordinate of central point.
    j : int
        y coordinate of central point.
    apertureDiameter : int
        Aperture diameter.

    Returns
    -------
    xcoords : arr
        array of x coords in aperture.
    ycoords : arr
        array of y coords in aperture.
    numpixels : int
        number of pixels in aperture.

    '''
    
    currentx = i - (apertureDiameter-1)/2 - 1       #start sweep in bottom left corner of square
    currenty = j - (apertureDiameter-1)/2 - 1
    xcoords = np.array([])                          #array to hold x coords
    ycoords = np.array([])                          #array to hold y coords
    for k in range(apertureDiameter):
        currentx += 1                               #iterate through x
        currenty = j - (apertureDiameter-1)/2 - 1
        for h in range(apertureDiameter):
            currenty += 1                           #iterate through y
            if round(np.sqrt((currentx-i)**2+(currenty-j)**2)) <= (apertureDiameter-1)/2: #check if cel; is within aperture
                xcoords = np.append(xcoords, currentx)
                ycoords = np.append(ycoords, currenty)
    numpixels = len(xcoords)                        #count the number of pixels
    return xcoords,ycoords,numpixels   


def drawRing(i,j,apertureDiameter):
    '''
    Function to return the pixels within the ring centered on i,j.

    Parameters
    ----------
    i : int
        x coordinate of central point.
    j : int
        y coordinate of central point.
    apertureDiameter : int
        Aperture diameter.

    Returns
    -------
    xcoords : arr
        array of x coords in aperture.
    ycoords : arr
        array of y coords in aperture.
    numpixels : int
        number of pixels in aperture

    '''
    
    currentx = i - (apertureDiameter-1)/2 - 1       #start sweep in bottom left corner of square
    currenty = j - (apertureDiameter-1)/2 - 1
    xcoords = np.array([])                          #array to hold x coords
    ycoords = np.array([])
    for k in range(apertureDiameter):
        currentx += 1                               #iterate through x
        currenty = j - (apertureDiameter-1)/2 - 1
        for h in range(apertureDiameter):
            currenty += 1                           #iterate through y
            if round(np.sqrt((currentx-i)**2+(currenty-j)**2)) == (apertureDiameter-1)/2: #check if cel; is within aperture
                xcoords = np.append(xcoords, currentx)
                ycoords = np.append(ycoords, currenty)
    numpixels = len(xcoords)
    return xcoords,ycoords,numpixels


def findMax(method, empty, minpix, aperture):
    '''
    Function to find the maximum pixel value that is not masked.
    Will return total flux within aperture centered on this pixel.
    Will mask all measured pixels.

    Parameters
    ----------
    method : str
        variable or fixed.
    empty : flt
        fraction of ring empty for variabe method.
    minpix : int
        minimum object size in pixels.
    aperture : int
        aperture diameter in pixels.

    Returns
    -------
    valueTotal : int
        total pixel value in aperture.
    maxValue : int
        maximum pixel value in aperture.
    maxI : int
        x coordinate of central pixel.
    maxJ : int
        y coordinate of central pixel.
    numpixels : int
        number of pixels in aperture.
    localBackground : flt
        background at object.

    '''
    global dataRead
    global dataMasked
    global mu
    
    apertureDiameter = -1
    valueTotal = 0          #variable to hold the value within aperture
    maxValue = np.max(dataRead)
    xarr, yarr = np.where(dataRead == maxValue)
    maxI = xarr[0]                                                                       #record x and y
    maxJ = yarr[0]
    
    if method == 'variable':
        for i in range(100):
            pixelCounter = 0
            apertureDiameter += 2
            xcoords,ycoords,numpixels = drawRing(maxI,maxJ,apertureDiameter)
            for k in range(numpixels):
                if xcoords[k] < np.shape(dataRead)[0] and ycoords[k] < np.shape(dataRead)[1]:
                    if xcoords[k] >= 0 and ycoords[k] >= 0:
                        if dataRead[int(xcoords[k])][int(ycoords[k])] <= mu + 3* sigma:
                            pixelCounter+=1
            if numpixels * empty <= pixelCounter:
                break
    
    if method == 'fixed':
        apertureDiameter = aperture
                    
    xcoords,ycoords,numpixels = inAperture(maxI,maxJ,apertureDiameter)                  #find pixels within aperture around i,j
    for k in range(numpixels):
        if xcoords[k] < np.shape(dataRead)[0] and ycoords[k] < np.shape(dataRead)[1]:
            if xcoords[k] >= 0 and ycoords[k] >= 0:
                valueTotal += dataRead[int(xcoords[k])][int(ycoords[k])]                #find total aperture value
                dataMasked[int(xcoords[k])][int(ycoords[k])] = np.ma.masked             #update mask array
                
            
    backgroundDiameter = apertureDiameter + 16
    xcoordsBack,ycoordsBack,numpixelsBack = inAperture(maxI,maxJ,backgroundDiameter)
    backgroundValue = 0
    backCount = 0
    for k in range(numpixelsBack):
        if xcoordsBack[k] < np.shape(dataRead)[0] and ycoordsBack[k] < np.shape(dataRead)[1]:
            if xcoordsBack[k] >= 0 and ycoordsBack[k] >= 0:
                if dataRead[int(xcoordsBack[k])][int(ycoordsBack[k])] <  mu + 3* sigma:
                    backgroundValue += dataRead[int(xcoordsBack[k])][int(ycoordsBack[k])]
                    backCount+=1
    if backgroundValue == 0:
        localBackground = mu
    else:
        localBackground = backgroundValue/backCount
    dataRead = dataMasked.filled(localBackground)
    
    if numpixels <= minpix:
        valueTotal = 0

    return (valueTotal,maxValue, maxI, maxJ, numpixels,localBackground)


def counter(method = 'variable', empty = 0.9, minpix = 2, aperture = 13):
    '''
    Measures all pixel values of objects.
    
    Parameters
    ----------
    method : str
        variable or fixed.
    empty : flt
        fraction of ring empty for variabe method.
    minpix : int
        minimum object size in pixels.
    aperture : int
        aperture diameter in pixels.

    Returns
    -------
    ValueCatalogue : arr
        array of total values measured.
    MagCatalogue : arr
        array of magnitude values measured.
    Xcoords : arr
        array x coords of objects.
    Ycoords : arr
        array y coords of objects.
    localBackgrounds : arr
        array of local background values.

    '''
    global mu
    global sigma
    
    ValueCatalogue = np.array([])   #array to hold pixel values of sources
    MagCatalogue = np.array([])   #array to hold magnitude values of sources
    Ycoords = np.array([])
    Xcoords = np.array([])
    numpixelsArr = np.array([])
    localBackgrounds = np.array([])

    for i in range(1000000):
        valueTotal,maxValue, maxI, maxJ, numpixels,localBackground = findMax(method,empty,minpix, aperture)      #run findMax
        Value = valueTotal - (localBackground * numpixels)                  #adjust for background
        if Value > 0:
            ValueCatalogue = np.append(ValueCatalogue, [Value])     #store value
            Ycoords = np.append(Ycoords,maxI)
            Xcoords = np.append(Xcoords,maxJ)
            numpixelsArr = np.append(numpixelsArr,numpixels)
            localBackgrounds = np.append(localBackgrounds, localBackground)
            m_uncorrected = -2.5*np.log10(Value) #the magnitudes without adding the zeropoint
            m = zeropoint + m_uncorrected #the magnitudes after adding the zeropoint
            MagCatalogue = np.append(MagCatalogue, m)
        if maxValue <= mu + 3* sigma:
            break
    return ValueCatalogue, MagCatalogue, Xcoords, Ycoords, localBackgrounds, numpixelsArr



#%%
'''
_________________________________________________________________________________________________

Testing the Source Detector
_________________________________________________________________________________________________

'''

def EmptyArrayTest():
    '''
    Tests an empty array

    Returns
    -------
    None.

    '''
    global dataMasked
    global dataRead
    
    
    data = np.zeros([1000,1000])
    mask = np.zeros([1000,1000])   
    dataMasked = np.ma.masked_array(data,mask)
    dataRead = dataMasked.filled(mu)  
     
    print('EMPTY ARRAY TEST')
    print('________________')
    print('\n')
    Vals, Mags,  Xcoords, Ycoords, localBackgrounds, numpixelsArr = counter()                                   
    print(len(Vals), 'objects Detected')
    print('\n')
    print('\n')
    
    return

def OnePixelTest():
    '''
    Tests a single pixel

    Returns
    -------
    None.

    '''
    global dataMasked
    global dataRead
    
    data = np.zeros([1000,1000])
    mask = np.zeros([1000,1000])
    data[500][500] = 1000000   
    dataMasked = np.ma.masked_array(data,mask)
    dataRead = dataMasked.filled(mu)   
    print('SINGLE PIXEL TEST')
    print('_________________')
    print('\n')
    Vals, Mags,  Xcoords, Ycoords, localBackgrounds, numpixelsArr = counter()
    print(Vals)                                       
    print(len(Vals), 'objects Detected')
    print('\n')
    print('\n')
    
    return

def GaussBlob(arr, x,y,std):
    '''
    Create Gaussian Blob Image

    Parameters
    ----------
    arr : arr
        image array.
    x : int
        x coord of center.
    y : int
        y coord of center.
    std : flt
        standard deviation of gaussian.

    Returns
    -------
    NONE.

    '''
    
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

    for i in range(-100,100):
        for j in range(-100,100):
            dist = np.sqrt((i)**2 + (j)**2)
            arr[x+i][y+j] = gauss(dist, 0, std,4e6)
            
    return

def GaussBlobTest():
    
    global dataMasked
    global dataRead
    
    data = np.zeros([500,500])
    GaussBlob(data, 250, 250, 10)
    GaussBlob(data, 300, 50, 2)
    GaussBlob(data, 100, 300, 2)
    GaussBlob(data, 400, 400, 2)
    GaussBlob(data, 150, 40, 2)
    mask = np.zeros([500,500])   
    dataMasked = np.ma.masked_array(data,mask)
    dataRead = dataMasked.filled(mu) 
                    
    plt.figure()
    plt.imshow(data, cmap='gray')                  
    plt.title('Gaussian Blob Test Image')
    plt.colorbar()                                    
     
    print('GAUSSIAN BLOB TEST IMAGE')
    print('________________________')
    print('\n')
    Vals, Mags,  Xcoords, Ycoords, localBackgrounds, numpixelsArr = counter()
    print(Vals)                                       
    print(len(Vals), 'objects Detected')
    print('\n')
    print('\n')
    
    return

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

def SmallImageMask1():
    '''
    Crops to small section of image.

    Returns
    -------
    None.

    '''
    global dataMasked
    global dataRead
    
    mask = np.zeros([int(np.shape(data)[0]),int(np.shape(data)[1])])   #array to hold masked pixels
    
    drawRect(mask, 1, [0,0], [339,2570])                               #cropping
    drawRect(mask, 1, [0,0], [4611,1847])
    drawRect(mask, 1, [0,2018], [4611,2570])
    drawRect(mask, 1, [552,0], [4611,2570])
    
    dataMasked = np.ma.masked_array(data,mask)
    
    dataRead = dataMasked.filled(mu)
    return


def RunSmallImageTest1(method = 'variable',empty = 0.9, minpix = 1, aperture = 13):
    '''
    Runs count on small section of image.

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
    
    dataReset()
    
    SmallImageMask1()
    
    Vals, Mags,  Xcoords, Ycoords, localBackgrounds, numpixelsArr = counter(method,empty,minpix, aperture)

    MagPlot = []
    for i in range(len(Mags)):
        if Mags[i]>9 and Mags[i]<=17:
            MagPlot.append(Mags[i])
            
    BINS = 30                                             #no. bins for hist
    values, bins = np.histogram(Mags, BINS)
    cumulative = np.cumsum(values) 
    
    BINSfit = 30                                          #no. bins for hist
    valuesfit, binsfit = np.histogram(MagPlot, BINSfit)
    cumulativefit = np.cumsum(valuesfit) 
    popt, pcov = curve_fit(fit, binsfit[:-1]+(binsfit[1]-binsfit[0])/2,cumulativefit, sigma = np.sqrt(cumulativefit), maxfev=10000, p0 = np.array([0.6, 2]))
    error = np.sqrt(np.diag(pcov))
    
    number = len(Mags)
    return (bins, cumulative, popt, error, number,dataMasked)    

def SmallImageMask2():
    '''
    Crops to small section of image.

    Returns
    -------
    None.

    '''
    global dataMasked
    global dataRead
    
    mask = np.zeros([int(np.shape(data)[0]),int(np.shape(data)[1])])   #array to hold masked pixels
    
    drawRect(mask, 1, [0,0], [1139,2570])                               #cropping
    drawRect(mask, 1, [0,0], [4611,709])
    drawRect(mask, 1, [1261,0], [4611,2570])
    drawRect(mask, 1, [0,975], [4611,2570])
    
    dataMasked = np.ma.masked_array(data,mask)
    
    dataRead = dataMasked.filled(mu)
    return


def RunSmallImageTest2(method = 'variable',empty = 0.9, minpix = 1, aperture = 13):
    '''
    Runs count on small section of image.

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
    
    dataReset()
    
    SmallImageMask2()
    
    Vals, Mags,  Xcoords, Ycoords, localBackgrounds, numpixelsArr = counter(method, empty, minpix, aperture)

    MagPlot = []
    for i in range(len(Mags)):
        if Mags[i]>9 and Mags[i]<=17:
            MagPlot.append(Mags[i])
            
    BINS = 30                                             #no. bins for hist
    values, bins = np.histogram(Mags, BINS)
    cumulative = np.cumsum(values) 
    
    BINSfit = 30                                          #no. bins for hist
    valuesfit, binsfit = np.histogram(MagPlot, BINSfit)
    cumulativefit = np.cumsum(valuesfit) 
    popt, pcov = curve_fit(fit, binsfit[:-1]+(binsfit[1]-binsfit[0])/2,cumulativefit, sigma = np.sqrt(cumulativefit), maxfev=10000, p0 = np.array([0.6, 2]))
    error = np.sqrt(np.diag(pcov))
    
    number = len(Mags)
    return (bins, cumulative, popt, error, number, dataMasked)
