import cv2
import math as m
import numpy as np

def gaussianKernel(sigma):
    """ Creates a 2D Gaussian-like Kernel with standard deviation sigma,
    centered on the middle of the kernel.
    Its size is taken as 4*sigma+1.

    Args:
        sigma (double): standard deviation

    Returns:
        numpy array: 2D array of size (4*sigma+1)x(4*sigma+1).
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
        modelizes a closed contour - See references
        eps (double): parameter used in dirac and heaviside functions approximation
        s (int): standard deviation of gaussian kernels
        epsilon (double): strictly positive constant
        l1, l2 (double): strictly positive constants 

    Returns:
        tuple: (c1,c2,f1,f2, LIF, GIF) where c1, c2 are doubles and
        f1, f2, LIF, GIF are 2D numpy arrays

    References:
        Based on [ Li Wang, Chunming Li, Quansen Sun, Deshen Xia, Chiu-Yen Kao,
        "Active contours driven by local and global intensity fitting energy with
        application to brain MR image segmentation" Computerized Medical Imaging 
        and Graphics 33 (2009) 520–531 ]
    """

    Id = np.ones(I.shape)
    K = gaussianKernel(s)
    H = (1 + np.arctan(previousPhi/eps)*2/math.pi)/2
    HI = H*I
    oneMH = 1-H
    oneMHI = oneMH*I
    
    c1 = np.sum(HI)/np.sum(H)
    c2 = np.sum(oneMHI)/np.sum(oneMH)
    
    KHI = cv2.GaussianBlur(HI, ksize = (2*int(s)+1,2*int(s)+1), sigmaX = s, sigmaY =s, borderType=cv2.BORDER_DEFAULT)
    KH = cv2.GaussianBlur(H, ksize = (2*int(s)+1,2*int(s)+1),sigmaX = s,sigmaY = s, borderType=cv2.BORDER_DEFAULT)
    K1HI = cv2.GaussianBlur(oneMHI, ksize = (2*int(s)+1,2*int(s)+1),sigmaX =s, sigmaY =s, borderType=cv2.BORDER_DEFAULT)
    K1H = cv2.GaussianBlur(oneMH, ksize = (2*int(s)+1,2*int(s)+1), sigmaX =s, sigmaY =s, borderType=cv2.BORDER_DEFAULT)
    
    f1 = KHI/KH
    f2 = K1HI/K1H
    
    GIF = -l1 * (I- c1)*(I-c1) + l2 * (I-c2)*(I-c2)
    
    LIF = (l2 - l1) * I*I * signal.fftconvolve(Id, K, mode='same')\
    -2 * I * signal.fftconvolve(-l1 * f1 + l2 * f2, K, mode='same')\
    + signal.fftconvolve(-l1 * f1*f1 + l2 * f2*f2, K, mode='same')
    
    return c1, c2, f1, f2, LIF, GIF

def initiateContour(I, typeC='circle', pointIni = None, setPoints = None):
    """Create an initial contour for I as a zero level set function. The
    shape of this contour depends on typeC.

    Args:
        I (array_like): one canal image
        typeC (string): type of contour to draw. It can be 'circle',
        or 'set_of_points' (still in development).
        pointIni (tuple or list, optional): center of the circle if
        'circle' is chosen
        setPoints (array or list, optional): list of points to interpolate
        so that to draw a closed curve from the points (if 'set_of_points' is chosen)

    Returns:
        array_like ofsame size than I. The initial contour is caracterized by zero values
        in the array.
    """
    initialPhi = np.zeros(I.shape)
    if typeC=='circle':
        if pointIni is None:
            raise TypeError('Missing center point to draw circle')
        else:
            for i in range(I.shape[0]):
                for j in range(I.shape[1]):
                    initialPhi[i,j] = -3 + np.sqrt((pointIni[0]-i)**2\
                                       + (pointIni[1]-j)**2)
   
    if typeC == 'set_of_points':
        print('not ready')

        # if setPoints is None:
        #     raise TypeError('Missing list setPoints')
        # else:
            
        #     pad = 3
        #     setPoints = np.pad(setPoints, [(pad,pad), (0,0)], mode='wrap')
        #     x, y = setPoints.T
        #     i = np.arange(0, len(setPoints))
            
        #     interp_i = np.linspace(pad, i.max() - pad + 1, 5 * (i.size - 2*pad))
        #     xi = interp1d(i, x, kind='cubic')(interp_i) #ou interp2D
        #     yi = interp1d(i, y, kind='cubic')(interp_i)
        #     contour = np.array([])#np array 2D where a line is a point from the contour
        #     for i in range(I.shape[0]):
        #         for j in range(I.shape[1]):
        #             initialPhi[i,j] = cv2.pointPolygonTest(contour, [i,j], True)

    return initialPhi

def activeContour(I, contourIni, thresh, l1, l2, s, eps, mu, nu, dt):
    """Determines the contour of an object in image I by recurrence. This function
    acts like a snake function. The contour evoluates from the input contourIni until 
    a stationary solution is found, that is:
    |new contour - previous contour| < thresh.
    The chosen norm is the maximum among |new contour (x,y) - previous contour (x,y)| for
    (x,y) in the image domain.
    The evolution of the contour is provided by :
            dC/dt = dirac(C) * (F1 + F2) + nu * dirac(C) * div( grad(C)/|grad(C)|)
                    + mu * (lap(C) - div( grad(C)/|grad(C)|) )

    Args:
        I (array_like): one-canal image
        contourIni (array_like of same size than I): contains the zero level
        set initial values
        l1 (double): stricly positive constant
        l2 (double): stricly positive constant
        s (double): standard deviation for gaussian kernels 
        mu (double): stricly positive constant that weights the level set regularization term
        nu (double): stricly positive constant that weights the length term
        dt (double): time step
        eps (double): strictly positive constant
        thresh (double): stopping condition for the recurrence

    Returns:
        array_like: stationary contour. The array has the same size than I.
    
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
    while stop_criterion > thresh and step <= 1000:
        
        c1, c2, f1, f2, GIF, LIF = intensities(I, previousPhi, eps, s, l1, l2)

        w = 0.01
        dirac = eps/(math.pi*(eps*eps+previousPhi*previousPhi))        
        grad = np.gradient(previousPhi)
        normgrad = np.sqrt(grad[0]*grad[0] + grad[1]*grad[1])
        normgrad = normgrad + (normgrad == 0)*1.
        div = np.gradient(grad[0]/normgrad, axis = 0) + np.gradient(grad[1]/normgrad, axis = 1)
        lap = cv2.Laplacian(previousPhi, cv2.CV_64F)
        
        dPhi = dirac * ((1-w)*LIF + w*GIF)\
                + nu * dirac * div\
                + mu * (lap - div)
        
        newPhi = previousPhi + dt * dPhi
        
        stop_criterion = np.max(abs(dPhi * dt))
        previousPhi = newPhi
    
        step = step+1
    
    return previousPhi

def extractContour(levelSet, image):
    listC=[]
    I = image
    for x in range(1,levelSet.shape[0]-1):
        for y in range(1,levelSet.shape[1]-1):
            if levelSet[x,y]>0:
                if (levelSet[x+1,y]<0 or levelSet[x+1,y-1]<0\
                or levelSet[x+1,y+1]<0 or levelSet[x-1,y]<0\
                or levelSet[x-1,y-1]<0 or levelSet[x-1,y+1]<0\
                or levelSet[x,y-1]<0 or levelSet[x,y+1]<0):
                    listC.append((x,y))
                    I[x,y,:] = [0,255,0]
                
            elif levelSet[x,y]<0:
                if (levelSet[x+1,y]>0 or levelSet[x+1,y-1]>0\
                or levelSet[x+1,y+1]>0 or levelSet[x-1,y]>0\
                or levelSet[x-1,y-1]>0 or levelSet[x-1,y+1]>0\
                or levelSet[x,y-1]>0 or levelSet[x,y+1]>0):
                    listC.append((x,y))
                    I[x,y,:] = [0,255,0]
                    
            else:
                listC.append((x,y))
                I[x,y,:] = [0,255,0]
    
    return I, listC    
