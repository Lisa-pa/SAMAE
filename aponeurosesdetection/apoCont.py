"""Functions for the active contour model - detection of aponeuroses"""
import cv2
import math
import numpy as np

def gaussianKernel(sigma):
    """ Creates a 2D Gaussian-like Kernel with standard deviation sigma,
    centered on the middle of the kernel.
    Its size is taken as 4*sigma+1.

    Args:
        sigma (double): standard deviation

    Returns:
        numpy array: kernel as a 2D array of size (4*sigma+1)x(4*sigma+1).
    """
    sizeK = 4 * int(sigma) + 1
    a = np.ogrid[-int(sizeK/2):int(sizeK/2)+1, -int(sizeK/2):int(sizeK/2)+1][0]*np.ones((sizeK, sizeK))
    b = np.ogrid[-int(sizeK/2):int(sizeK/2)+1, -int(sizeK/2):int(sizeK/2)+1][1]*np.ones((sizeK, sizeK))
    kernel = 1/(2*math.pi*sigma*sigma) * np.exp(-(a**2 + b**2)/(2*sigma*sigma))
    return kernel

def intensities(I, previousPhi, eps, s, l1, l2):
    """Calculates global and local intensities in two regions
    of image I. These two regions are delineated by a contour
    modeled as the zero level set of a Lipschitz function Phi.
    
        - c1 is the global intensity outside the contour
        - c2 is the global intensity inside the contour
        - f1 is a matrix of local intensities outside the contour
        - f2 is a matrix of local intensities inside the contour
        - LIF is a matrix where each coefficient is the local 
        intensity force of the corresponding pixel in I
        - GIF is a matrix where each coefficient is the global 
        intensity force of the corresponding pixel in I

    Args:
        I (numpy array): 2D image
        previousPhi (numpy array): 2D matrix of a level set function that
            modelizes a closed contour - We advice you to read the references
        eps (double): parameter used in dirac and heaviside 
            functions approximation, stricly positive
        s (int): standard deviation of gaussian kernels
        l1, l2 (double): strictly positive constants 

    Returns:
        c1, c2, f1, f2, LIF, GIF where c1, c2 are doubles and
        f1, f2, LIF, GIF are 2D numpy arrays

    References:
        Based on [ Li Wang, Chunming Li, Quansen Sun, Deshen Xia, Chiu-Yen Kao,
        "Active contours driven by local and global intensity fitting energy with
        application to brain MR image segmentation" Computerized Medical Imaging 
        and Graphics 33 (2009) 520–531 ]
    """
    import scipy.signal as signal
    Id = np.ones(I.shape)
    K = gaussianKernel(s)
    H = (1 + np.arctan(previousPhi/eps)*2/math.pi)/2
    HI = H*I
    oneMH = 1-H
    oneMHI = oneMH*I
    
    #global intensities
    c1 = np.sum(HI)/np.sum(H)
    c2 = np.sum(oneMHI)/np.sum(oneMH)
    
    #local intensities
    KHI = cv2.GaussianBlur(HI, ksize = (2*int(s)+1,2*int(s)+1), sigmaX = s, sigmaY =s, borderType=cv2.BORDER_DEFAULT)
    KH = cv2.GaussianBlur(H, ksize = (2*int(s)+1,2*int(s)+1),sigmaX = s,sigmaY = s, borderType=cv2.BORDER_DEFAULT)
    K1HI = cv2.GaussianBlur(oneMHI, ksize = (2*int(s)+1,2*int(s)+1),sigmaX =s, sigmaY =s, borderType=cv2.BORDER_DEFAULT)
    K1H = cv2.GaussianBlur(oneMH, ksize = (2*int(s)+1,2*int(s)+1), sigmaX =s, sigmaY =s, borderType=cv2.BORDER_DEFAULT)
    f1 = KHI/KH
    f2 = K1HI/K1H

    #Local and global intensity forces
    #convolutions are more efficiently computed in the frequency domain
    GIF = -l1 * (I- c1)*(I-c1) + l2 * (I-c2)*(I-c2)
    LIF = (l2 - l1) * I*I * signal.fftconvolve(Id, K, mode='same')\
    -2 * I * signal.fftconvolve(-l1 * f1 + l2 * f2, K, mode='same')\
    + signal.fftconvolve(-l1 * f1*f1 + l2 * f2*f2, K, mode='same')
    
    return c1, c2, f1, f2, LIF, GIF

def initiateContour(I, typeC, setPoints = None, param = None):
    """Create an initial contour for I, as a zero level set function. The
    shape of this contour depends on typeC.

    Args:
        I (array_like): one canal image

        typeC (string): type of contour to draw. It can be 'circle',
        or 'set_of_points' or 'quardrangle_param'.
            'circle' : draws a circle centered on the first point
            of setPoints, and which diameter is 3 pixels.
            'set_of_points': draws a convex hull containing all 
            points from setPoints. If setPoints is None, an error
            is raised.
            'quadrangle_param': draws a quadrangle of width w centered on a
            line of equation x = ay+b ; x : rows of I, y : columns of I

        setPoints (array or list): list of points drawing a closed 
        curve if typeC is 'set_of_points'; or list containing one 
        single point if typeC is 'circle' (center of the circle). 
        If option is 'circle' and setPoints has more than one point, 
        the first one is considered as the center of the circle.
        If option is 'set_of_points', considering that x-axis goes
        downwards and y-axis goes right, points must be given
        in a non clockwise manner
        If option is 'quadrangle_param', this object is used to store
        the four summets' coordinates, whatever its input value.
        
        param (1D array or list): a = param[0], b = param[1] are the
        parameters for the line x = ay + b. param[2] defines the desired
        width of the quadrangle

    Returns:
        array_like of same size than I. Pixels on the contour are zeros,
        pixels inside the contour are negative, pixels outside the
        contour are positive.
    """

    C = np.ones(I.shape)
            
    if typeC =='circle':
        if setPoints is None:
            raise ValueError('Missing setPoints list.')
        for i in range(I.shape[0]):
            for j in range(I.shape[1]):
                C[i,j] = -3 + np.sqrt((setPoints[0][0]-i)**2\
                                + (setPoints[0][1]-j)**2)

    if typeC == 'set_of_points':
        if setPoints is None:
            raise ValueError('Missing setPoints list.')
        contour = cv2.convexHull(setPoints,  clockwise = False)
        for i in range(I.shape[0]):
            for j in range(I.shape[1]):
                C[i,j] = - cv2.pointPolygonTest(contour, (i,j), True)
    
    if typeC == 'quadrangle_param':
        if param is None:
            raise ValueError('Missing list containing line equation parameters.')
        a = param[0]
        b = param[1]
        w = int(param[2]/2) # half width of the quadrangle, in pixels
        
        y = np.arange(0, I.shape[1],1)
        x_m = np.int64(a*y + b - w)
        x_p = np.int64(a*y + b + w)
        setPoints = []
        
        n1 = 10
        while (x_m[n1]<0 or x_m[n1]>I.shape[0]-1) and n1 < I.shape[1]-1:
            n1 = n1 + 1
        
        n2 = I.shape[1] - 1
        while (x_m[n2]<0 or x_m[n2]>I.shape[0]-1) and n2 >= 0:
            n2 = n2 - 1
        
        n3 = 10
        while (x_p[n3]<0 or x_p[n3]>I.shape[0]-1) and n3 < I.shape[1]-1:
            n3 = n3 + 1

        n4 = I.shape[1] - 1
        while (x_p[n4]<0 or x_p[n4]>I.shape[0]-1) and n4 >= 0:
            n4 = n4 - 1

        setPoints = np.array([[x_m[n1],n1],[x_p[n3],n3],[x_p[n4],n4],[x_m[n2],n2]])
        contour = cv2.convexHull(setPoints,  clockwise = False)
        for i in range(I.shape[0]):
            for j in range(I.shape[1]):
                C[i,j] = - cv2.pointPolygonTest(contour, (i,j), True)
    return C


def activeContour(I, contourIni, thresh, l1, l2, s, eps, mu, nu, dt):
    """Determines the contour of an object in image I by recurrence. This function
    acts like a snake function. The contour evoluates from the input contourIni until 
    a stationary solution is found, such as:
    |new contour - previous contour| < thresh.
    The chosen norm is the maximum among |new contour (x,y) - previous contour (x,y)| for
    (x,y) in the image domain.
    The evolution of the contour is provided by :
            dC/dt = dirac(C) * (F1 + F2) + nu * dirac(C) * div( grad(C)/|grad(C)|)
                    + mu * (lap(C) - div( grad(C)/|grad(C)|) )
    See reference articles in References for the further detailed approach.

    Args:
        I (array_like): one-canal image
        contourIni (array_like of same size than I): contains the zero level
        set initial values. Same format as the output of function initiateContour.
        l1 (double): stricly positive constant
        l2 (double): stricly positive constant
        s (double): standard deviation for gaussian kernels 
        mu (double): stricly positive constant that weights the level set regularization term
        nu (double): stricly positive constant that weights the length term
        dt (double): time step. It must be small enough to avoid numerical issues
        but big enough to reach a reasonable computational time
        eps (double): strictly positive constant
        thresh (double): stopping condition for the recurrence

    Returns:
        array_like: stationary contour represented by pixels = 0.
        The array has the same size than I.
    
    References:
        Li Wang et al., Computerized Medical Imaging and Graphics 33 (2009) 520–531
        'Active contours driven by local and global intensity fitting energy with
        application to brain MR image segmentation'
        Shan Ling et al., IEEE JOURNAL OF BIOMEDICAL AND HEALTH INFORMATICS,
        VOL. 17, NO. 6, NOVEMBER 2013. 'Automatic Tracking of Aponeuroses and 
        Estimation of Muscle Thickness in Ultrasonography: A Feasibility Study'
    """
    
    if len(I.shape) > 2:
        I = cv2.cvtColor(I, cv2.COLOR_RGB2GRAY)
        
    'Initialization'
    previousPhi = contourIni
    
    'recurrence'
    step = 1
    stop_criterion = thresh + 1.
    while stop_criterion > thresh and step <= 10000:
        if step%10 == 0:
            print('Current tens for step : ', step)
        
        #compute global and local intensity forces at step n
        c1, c2, f1, f2, GIF, LIF = intensities(I, previousPhi, eps, s, l1, l2)

        # constant weight for local and global intensity forces
        # SHOULD BE MODIFIED IN THE FUTURE TO ADAPT EACH IMAGE
        w = 0.01

        # compute mathematical operators at step n
        dirac = eps/(math.pi*(eps*eps+previousPhi*previousPhi))        
        grad = np.gradient(previousPhi)
        normgrad = np.sqrt(grad[0]*grad[0] + grad[1]*grad[1])
        normgrad = normgrad + (normgrad == 0)*1.
        div = np.gradient(grad[0]/normgrad, axis = 0) + np.gradient(grad[1]/normgrad, axis = 1)
        lap = cv2.Laplacian(previousPhi, cv2.CV_64F)
        
        # compute derivative of contour
        dPhi = dirac * ((1-w)*LIF + w*GIF)\
                + nu * dirac * div\
                + mu * (lap - div)
        
        # compute new contour
        newPhi = previousPhi + dt * dPhi
        
        # check the stopping criterion
        stop_criterion = np.max(abs(dPhi * dt))

        # go to next step n
        previousPhi = newPhi
        step = step+1

    return previousPhi, step-1

def extractContour(levelSet, image, offSetX = 0, offSetY = 0):
    """
    This function spots the border between negative and positive values of
    the levelSet object (which represents a zero level set function).
    This border represents the contour of an object in image I.
    The function returns a list containing the points of the contour and
    an image similar to I with the contour is green-lighted. 
    
    Args:
        image (array): 3-canal image
        levelSet (array): array of size (image.shape[0], image.shape[1])
        representing a zero level set function (same format as output
        by activeContour function)
        offSetX: add an vertical offset (changes rows)
        offSetY: add an horizontal offSet (changes columns)
    Returns:
        I (array): array of same size than I, with same pixels value except
        the contour which is green
        listC (list): list of all pixels contained in the detected contour
    """
    
    listC=[]
    I = np.copy(image) 
      
    # binarization of the levelset image to separate positive from negative values
    binar = np.uint8((levelSet>=0)*0. + (levelSet<0)*255.)
    #find contours
    objects = cv2.findContours(binar, mode = cv2.RETR_EXTERNAL, method = cv2.CHAIN_APPROX_NONE)[0]
    # verify that at least one object exists, otherwise it means that there is no any detected contour in levelset
    if len(objects) == 0 :
        raise ValueError('No contour has been found. Please check: 1) that aponeuroses\
                        have been correctly located. 2) that initial contour is correct')
    # If several contours detected in levelSet, keep only the biggest
    objects_size = [(i,objects[i].size) for i in range (len(objects))]
    objects_size.sort(key=lambda x:x[1]) 
    biggest = objects_size[-1][0]
    # verify that the points in the detected contours are within I's dimensions
    for point in objects[biggest]:
        if point[0][0] < I.shape[1] and point[0][1] < I.shape[0]:
            listC.append((point[0][1] + offSetX, point[0][0] + offSetY))
            I[point[0][1],point[0][0],:] = [0,255,0]
    return I, listC

def approximateApo(p, apoType, I, typeapprox, d):
    """Function that approximates an aponeurosis shape.
    More precisely, if apoType is 'upper', meaning an upper aponeurosis
    is being processed, this function will approximate the lower boundary
    of the upper aponeurosis with an interpolated curve of degree d. If apoType 
    is 'lower', meaning a deep aponeurosis is being processed, this function 
    will approximate the upper boundary of the deep aponeurosis with an 
    interpolated curve of degree d.
    The function returns the spline that approximates the aponeurosis.

    Args:
        p (list): list containing the points of the aponeurosis' contour
        apoType (string): aponeurosis type ('upper' or 'lower')
        I (array) : image to which the processed aponeurosis belongs. It
        can either have one canal or three canals.
        d (integer): degree of the interpolated curve.
        typeapprox (string): 'Bspline' or 'polyfit', to choose which
            type of interpolation you cant between B splines or polynomial fitting.

    Returns
        spline modeling aponeurosis.
    """

    p.sort(key=lambda x:x[0]) #x-coord sort
    p.sort(key=lambda x:x[1]) #y-coord sort

    # retrieve points from the contour depending on the type of aponeurosis
    if apoType == 'lower':
        #from top to bottom, keep the first contour pixel of each column 
        comparator = -1
        line = []
        for pt in p:
            if pt[1] != comparator:
                line.append(pt)
                comparator = pt[1]
        
    elif apoType == 'upper':
        #from bottom to top, keep the first contour pixel of each column 
        comparator = -1
        line = []
        for pt in reversed(p):
            if pt[1] != comparator:
                line.append(pt)
                comparator = pt[1]

    # sort aponeurosis' points (necessary for interpolation)
    line.sort(key=lambda x:x[1]) #y-coord sort
    ycoord = [line[i][1] for i in range(len(line))]
    xcoord = [line[i][0] for i in range(len(line))]

    # approximate by a spline, depending on the chosen type of interpolation
    if typeapprox == 'Bspline':
        import scipy.interpolate as interpolate
        spline = interpolate.UnivariateSpline(ycoord, xcoord, k=d, ext = 0) 
    elif typeapprox == 'polyfit':
        polyn = np.polyfit(ycoord, xcoord, deg = d)
        spline = np.poly1d(polyn)     
    
    return spline