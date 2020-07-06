import cv2
import math as m
import numpy as np


def heaviside(matrix, epsilon):
    """Calculates the smooth heaviside function for each
    coefficient of matrix.
        h : R -> [-pi/2 ; pi/2]
            x -> 1/2 * (1 + 2/pi * arctan( x /epsilon ))

    Args:
        matrix (numpy array): 2D array_like object which can be
        of any shape and contains real coefficients (doubles).
        epsilon (double): real and strictly positive constant. 

    Returns:
        numpy array: 2D array. Each coefficient is the heaviside value of
        the matrix's coefficient at same location.
    """
    if epsilon <=0:
        err = 'Error : epsilon must be strictly positive'
        print(err)
    else:
        h = np.zeros((matrix.shape))
        for u in range(matrix.shape[0]):
            h[u, :] = np.array([1/2.*(1+2.*m.atan(matrix[u,v]/epsilon)/m.pi) for v in range(matrix.shape[1])])
        return h


def gaussianKernel(ksize, point, s):
    """ Creates a 2D Gaussian-like Kernel of ksize with standard deviation s.
    The kernel is centered on point, and its coefficients decrease with distance 
    from point = (abs, ord):
            g : ZxZ -> R
                (x, y) -> 1/(2*pi*s*s) * exp( - ( (x-abs)**2 + (y-ord)**2 ) / (2*s*s) )

    Args:
        ksize (tuple or array_like): the two first coefficients of ksize give the size of the gaussian kernel.
        point (tuple or array_like): center point from which all distances are calculated
        s (double): standard deviation

    Returns:
        numpy array: 2D array. Each coefficient has a value reflecting its distance from
        point. The more distant from point, the smaller value.
    """

    absc = point[0]
    ordo = point[1]
    k = np.zeros(size)
    for a in range(size[0]):
        k[a, :] = np.array([1/(2*m.pi*s**2) * np.exp(-((a-absc)**2 + (b-ordo)**2)/(2*s*s))for b in range(size[1])])
    return k


def intensities(I, C, s, epsilon):
    """Calculates global and local intensities in two regions
    of image I. These two regions are delineated by the contour C.
        - c1 is the global intensity outside the contour
        - c2 is the global intensity inside the contour
        - f1 is a matrix of local intensities outside the contour
        - f2 is a matrix of local intensities inside the contour

    Args:
        I (numpy array): 2D image
        C (numpy array): 2D matrix of a level set function that
        modelizes a closed contour - See references
        s (int): standard deviation of gaussian kernels
        epsilon (double): strictly positive constant

    Returns:
        tuple: (c1,c2,f1,f2) where c1, c2 are doubles and
        f1, f2 are 2D numpy arrays

    References:
        Based on [ Li Wang, Chunming Li, Quansen Sun, Deshen Xia, Chiu-Yen Kao,
        "Active contours driven by local and global intensity fitting energy with
        application to brain MR image segmentation" Computerized Medical Imaging 
        and Graphics 33 (2009) 520–531 ]
    """

    heavisideMatrix = heaviside(C, epsilon)        
    OneMinusHeaviside = np.ones((I.shape[0],I.shape[1])) - heavisideMatrix
    HI = heavisideMatrix*I
    OneMinusH_I = OneMinusHeaviside*I

    #local intensities
    convolHI = cv2.GaussianBlur(HI, ksize = (2*int(s)+1,2*int(s)+1), sigmaX = s, sigmaY =s, borderType=cv2.BORDER_CONSTANT)
    convolH = cv2.GaussianBlur(heavisideMatrix, ksize = (2*int(s)+1,2*int(s)+1),sigmaX = s,sigmaY = s, borderType=cv2.BORDER_CONSTANT)
    convol1minusH_I = cv2.GaussianBlur(OneMinusH_I, ksize = (2*int(s)+1,2*int(s)+1),sigmaX =s, sigmaY =s, borderType=cv2.BORDER_CONSTANT)
    convol1minusH = cv2.GaussianBlur(OneMinusHeaviside, ksize = (2*int(s)+1,2*int(s)+1), sigmaX =s, sigmaY =s, borderType=cv2.BORDER_CONSTANT)
    
    f1 = convolHI/convolH
    f2 = convol1minusH_I/convol1minusH

    #global intensity
    c1 = np.sum(HI)/np.sum(heavisideMatrix)
    c2 = np.sum(OneMinusH_I)/np.sum(heavisideMatrix)
           
    return c1, c2, f1, f2


def apoContour(I, pointIni, l1, l2, mu, nu, dt, epsilon, omega, sigma, stop_thresh):
    """Determines the contour of an object in image I by recurrence. This function
    acts like a snake function. An initial contour is generated around the point
    pointIni. The contour then evoluates until a stationary solution is found, that is:
    |new contour - previous contour| < stop_thresh.
    The chosen norm is the maximum among |new contour (x,y) - previous contour (x,y)| for
    (x,y) in the image domain.
    The evolution of the contour is provided by :
            dC/dt = dirac(C) * (F1 + F2) + nu * dirac(C) * div( grad(C)/|grad(C)|)
                    + mu * (lap(C) - div( grad(C)/|grad(C)|) )

    Args:
        I (array_like): one-canal image
        pointIni (array_like or tuple): contains two coefficients = [x0, y0]
        l1 (double): stricly positive constant
        l2 (double): stricly positive constant
        mu (double): stricly positive constant that weights the level set regularization term
        nu (double): stricly positive constant that weights the length term
        dt (double): time step
        epsilon (double): strictly positive constant
        omega (double): constant between 0 and 1 that weights the influence of global and local intensities
        sigma (double): standard deviation for gaussian kernels 
        stop_thresh (double): stopping condition for the recurrence

    Returns:
        array_like: stationary contour
    
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
        
    'Initialization step'
    previousC = np.zeros((I.shape[0], I.shape[1])) #circle, center = point, diameter = 4
    for i in range(I.shape[0]):
        for j in range(I.shape[1]):
            previousC[i,j] = -2 + m.sqrt((i-pointIni[0])**2 + (j - pointIni[1])**2)
    test = 1000.
    step = 0
    
    'Recurrence'
    #I temporarily added a condition in while loop to limit the number of iterations
    while test> stop_thresh and step<=16:
        dirac = np.zeros((I.shape[0], I.shape[1]))
        
        #Vectorization for laplacian and divergence calculations
        #in the following names, P means Plus and M means Minus
        prevC_iP1 = np.concatenate((previousC[1:,:],np.zeros((1,previousC.shape[1]))), axis = 0)
        prevC_iM1 = np.concatenate((np.zeros((1,previousC.shape[1])),previousC[:previousC.shape[0]-1,:]), axis = 0)
        prevC_jP1 = np.concatenate((previousC[:,1:],np.zeros((previousC.shape[0],1))), axis = 1)
        prevC_jM1 = np.concatenate((np.zeros((previousC.shape[0],1)),previousC[:,:previousC.shape[1]-1]), axis = 1)
        prevC_ijM1 = np.concatenate((np.zeros((previousC.shape[0],1)), np.concatenate((np.zeros((1,previousC.shape[1]-1)), previousC[:-1,:-1]),axis = 0)),axis = 1)
        prevC_ijP1 = np.concatenate((np.concatenate((previousC[1:,1:],np.zeros((1,previousC.shape[1]-1))),axis = 0), np.zeros((previousC.shape[0],1))),axis = 1)
        prevC_iP1jM1 = np.concatenate((np.zeros((previousC.shape[0],1)), np.concatenate((previousC[1:,:-1], np.zeros((1,previousC.shape[1]-1))),axis = 0)),axis = 1)
        prevC_iM1jP1 = np.concatenate((np.concatenate((np.zeros((1,previousC.shape[1]-1)),previousC[:-1,1:]),axis = 0), np.zeros((previousC.shape[0],1))),axis = 1)

        norm = np.ones(previousC.shape)/np.sqrt(\
                      (prevC_iP1-previousC)*(prevC_iP1-previousC)\
                      +(prevC_jP1-previousC)*(prevC_jP1-previousC))
        
        divergence = (np.ones(previousC.shape)+1/4*(prevC_iP1-prevC_iM1)*(prevC_iP1-prevC_iM1))\
        *(prevC_iP1+prevC_iM1-2*previousC)\
        +(np.ones(previousC.shape)+1/4*(prevC_jP1-prevC_jM1)*(prevC_jP1-prevC_jM1))\
        *(prevC_jP1+prevC_jM1-2*previousC)\
        +1/8*(prevC_iP1-prevC_iM1)*(prevC_jP1-prevC_jM1)\
        *(prevC_ijP1-prevC_iP1jM1-prevC_iM1jP1-prevC_ijM1)
        
        div = divergence / norm
        
        lap = prevC_iP1 + prevC_iM1 + prevC_jP1 + prevC_jM1\
        -4*previousC
                
        #GLOBAL AND LOCAL INTENSITIES at current step
        c1, c2, f1, f2  = intensities(I, previousC, sigma, epsilon)
        GIF_Force = omega*(-l1*(I-c1*np.ones(I.shape))*(I-c1*np.ones(I.shape))\
                       + l2 * (I - c2*np.ones(I.shape))*(I - c2*np.ones(I.shape)))
        LIF_Force = np.zeros(I.shape)

        for i in range(previousC.shape[0]):
            #TEMPORARY : to follow what is going on
            print('step=',step)
            percent1 = i/previousC.shape[0]*100
            print('i=', i, percent1)
            #
            for j in range(previousC.shape[1]):             
                   
                K_ij = gaussianKernel(I, [i,j], sigma)
                I_f1_square_ij = ((I[i,j]*np.ones(I.shape))-f1)*((I[i,j]*np.ones(I.shape))-f1)
                I_f2_square_ij = ((I[i,j]*np.ones(I.shape))-f2)*((I[i,j]*np.ones(I.shape))-f2)
                LIF_Force[i,j] = -l1*np.sum((K_ij*I_f1_square_ij))+l2*np.sum((K_ij*I_f2_square_ij))
                LIF_Force[i,j] = (1-omega)*LIF_Force[i,j]
                
                dirac[i,j] = epsilon / (m.pi*(epsilon**2+previousC[i,j]**2))
        
        #time derivative approximation
        d_phi = dirac*((1 - omega)*LIF_Force+omega*GIF_Force)\
        +nu*dirac*div + mu*(lap - div)
        
        #what happens at borders ? temporary solution-> enforce no variations
        d_phi[0:1,:]=np.zeros((1,d_phi.shape[1]))
        d_phi[d_phi.shape[0]-1:d_phi.shape[0],:]=np.zeros((1,d_phi.shape[1]))
        d_phi[:,0:1]=np.zeros((d_phi.shape[0],1))
        d_phi[:,d_phi.shape[1]-1:d_phi.shape[1]]=np.zeros((d_phi.shape[0],1))
        
        #New contour 
        newC = d_phi*dt + previousC
                         
        #condition re-calculated for while loop
        test = np.max(np.absolute(d_phi*dt)) #is there a better test ?
        print('test=',test)
        
        #replace previous contour by new one
        previousC = newC
        step = step+1
    print(step, test)
    return previousC