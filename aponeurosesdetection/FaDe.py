"""Functions for fascicle detection - FaDe"""

import math as m
import cv2
import numpy as np
#------------------------------------------------------------------------------#
def d2_gaussianMasks(s):
    """Implements Gaussian's second derivatives masks for a given standard
    deviation s.
    The masks' size is set as 6*s rounded to the superior integer (+1 if even).

    Args:
        s (double) : standard deviation

    Returns:
        mask_xx, mask_xy, mask_yy (array_like): masks obtained from derivation of a gaussian,
        according to x and x, to x and y, to y and y respectively.

    """
    
    'Calculation of the size of the masks'
    size = m.ceil(s*6)
    if size%2==0: #if size if even, make it odd
        size+=1
        
    half =  m.floor(size/2)
    mask_xx = np.zeros((size,size))
    mask_xy = np.zeros((size,size))
    mask_yy = np.zeros((size,size))

    'Second derivatives expressions'
    for x in range (-half, half +1):
        for y in range (-half, half +1):

            mask_xx[x + half][y + half] = (x**2/s**2-1.)/(2.*m.pi*s**4)\
                                                *m.exp(-(x**2+y**2)/(2.*s**2))
                                                
            mask_xy[x + half][y + half] = x*y/(2.*m.pi*s**6)\
                                                *m.exp(-(x**2+y**2)/(2.*s**2))
                                                
            mask_yy[x + half][y + half] = (y**2/s**2-1.)/(2.*m.pi*s**4)\
                                                *m.exp(-(x**2+y**2)/(2.*s**2))
    return mask_xx, mask_xy, mask_yy

'-----------------------------------------------------------------------------'

def MVEF_2D(I, scales, thresholds):
    """Multiscale Vessel Enhancement Method for 2D images - based on Frangi's,
    Rana's, and Jalborg's works. This method searches for geometrical 
    structures which can be regarded as tubular.
    
    Args:
        I (2D array):       I is a grayscale image (otherwise it is 
                            converted to grayscale)
        thresholds (list):  thresholds is a list of two thresholds that 
                            control the sensitivity of the line filter 
                            to the measures Fr and R.  
        scales (list):      scales is a list of lengths that correspond to 
                            the diameter of the tube-like to find
    
    Outputs:
        I2 (2D array):  image I filtered by the multiscale vessel enhancement
                        filter
    References:
        Based on [ Automatic detection of skeletal muscle architecture 
        features, Frida Elen Jalborg, Master’s Thesis Spring 2016 ]
                [Frangi, 1998, Multiscale Vessel Enhancement filtering]
                [ Rana et al., 2009, Automated tracking of muscle fascicle
                orientation in B-mode ultrasound images]
                [ Rana et al., 2011, In-vivo determination of 3D muscle
                architecture of human muscle using free hand ultrasound]

    """
    from skimage.feature import hessian_matrix, hessian_matrix_eigvals
    
    if len(I.shape)>2:
        I = cv2.cvtColor(I, cv2.COLOR_RGB2GRAY)


    vesselness = np.zeros((I.shape[0], I.shape[1], len(scales)))
    I2 = np.zeros((I.shape[0], I.shape[1]))
    b = thresholds[0]
    c = thresholds[1]
    if b == 0:
        raise ValueError('first element of thresholds cannot be null')
        
    for sc in range(len(scales)): 
        H = hessian_matrix(image = I, sigma = scales[sc], order = 'rc')
        eigvals = hessian_matrix_eigvals(H) #2 eigenvalues in decreasing order; shape of ei = (2,I.shape[0], I.shape[1])
        frobeniusnorm =  np.sqrt(H[0]**2+ H[1]**2 + 2 * H[2]**2)
        if c == 0:
            c = 1 / 2 * np.amax(frobeniusnorm)

        for i in range(I.shape[0]):
            for j in range(I.shape[1]):

                'Calculation of vesselness: looking for a POSITIVE highest'
                'eigenvalue <=> looking for dark tube-like structures; looking'
                'for a NEGATIVE highest eigen value <=> looking for bright'
                'tube-like structures'
                #bright tube-like structures search did not work so the 
                #temporary solution is to inverse the ultrasound image 
                #and search for dark tube-like structures
                
                #if eigvals[0,i,j]<=0:
                    #vesselness[i,j,sc] = 0
                    #not necessary since vesselness initialized with zeros

                if eigvals[0,i,j]>0:
                    
                    'ratio - for second order ellipsoid'
                    R = eigvals[1,i,j]/eigvals[0,i,j]

                    'Frobenius norm - for second order structureness'
                    Fr = m.sqrt(eigvals[0,i,j]*eigvals[0,i,j] + eigvals[1,i,j]*eigvals[1,i,j])

                    vesselness[i,j,sc]=scales[sc]*m.exp(-R*R/(2.*b*b))*(1.-m.exp(-Fr*Fr/(2*c*c)))

    'Keep the highest value of vesselness across all scales'
    for ind1 in range(I2.shape[0]):
        for ind2 in range(I2.shape[1]):
            I2[ind1,ind2] = np.max(vesselness[ind1,ind2,:])
    
    return I2
'-----------------------------------------------------------------------------'
"""
from skimage.transform import radon, iradon
def normalizedRadon(I):


    if len(I.shape) > 2:
        I = cv2.cvtColor(I, cv2.COLOR_RGB2GRAY)
    
    #working with a square because radon function working on a circle only
    if I.shape[0] > I.shape[1]:
        I = I[int(I.shape[0]/2 - I.shape[1]/2):int(I.shape[0]/2 + I.shape[1]/2), :]
    else:
        I = I[:, int(I.shape[1]/2 - I.shape[0]/2):int(I.shape[1]/2 + I.shape[0]/2)]

    cv2.imshow('Part of the image analyzed', I)
    cv2.waitKey(0) & 0xFF
    cv2.destroyAllWindows()

    R = radon(I, circle = True)
    cv2.imshow('Radon Transform', R/np.max(R))
    cv2.waitKey(0) & 0xFF
    cv2.destroyAllWindows()
    
    X = np.arange(0,I.shape[0],1)
    Y = np.arange(0,I.shape[1],1)

    for rho in range(R.shape[0]):
        for theta in range(R.shape[1]):
            line = rho - m.cos(theta)*X + m.sin(theta)*Y
            #line = (np.abs(rho - m.cos(theta)*X + m.sin(theta)*Y) <0.6)*1. + (np.abs(rho - m.cos(theta)*X + m.sin(theta)*Y) >= 0.6)*0.
            dirac = 1./(m.pi*(1.*1.+line*line))        

            L = np.sum(dirac)
            if L != 0:
                R[rho, theta] = R[rho, theta]/L

    cv2.imshow('Radon Transform', R/np.max(R))
    cv2.waitKey(0) & 0xFF
    cv2.destroyAllWindows()
    
    R = cv2.threshold(R, 254., 255, cv2.THRESH_BINARY)[1]

    return R,I
"""
def contourAverage(inputC):
    """Function that approximate a closed contour by a line
    inputC has the same dimensions as the contours output by
    the cv2 function findContours

    ouputL : list of points [line, column]
    """
    outputL = []

    for i in range(inputC.shape[0]-1):
        check = len(outputL)
        for j in range(i+1, inputC.shape[0]):
            if inputC[i,0,0] == inputC[j,0,0]: #if two points on the same column, make average
                avera = [(inputC[i,0,1]+inputC[j,0,1])/2, (inputC[i,0,0]+inputC[j,0,0])/2]
                outputL.append(avera)
        if check == len(outputL): #if there is only one point on the column
            outputL.append([inputC[i,0,1],inputC[i,0,0]])
    return outputL


def locateFascicle(I, xcalib, ycalib, USimage):
    """
    I (array): binary image (one canal)
    """
    
    snippets = cv2.findContours(np.uint8(I), mode=cv2.RETR_EXTERNAL, method=cv2.CHAIN_APPROX_NONE)[0]
    if len(snippets) == 0:
        raise ValueError('No fascicle could be found in I.')
    
    filtered_indices = [] # list which contains the ID of the filtered snippets
    line_snip = [] #list which contains a linear approximation of the snippet filtered
    for i in range(len(snippets)):
        
    
        for u in range(snippets[i].shape[0]):
            USimage[snippets[i][u,0,1],snippets[i][u,0,0],:] = [0,0,255] #see what has be spotted
        
        #condition 1 - on length : minimum 4 mm long
        if snippets[i].shape[0]/2 > 4/m.sqrt(xcalib**2+ycalib**2):

            left_points = [snippets[i][m,0,:] for m in range(snippets[i].shape[0])\
                if snippets[i][m,0,0] == np.amin(snippets[i][:,0,0])]
            right_points = [snippets[i][m,0,:] for m in range(snippets[i].shape[0])\
                if snippets[i][m,0,0] == np.amax(snippets[i][:,0,0])]
            if len(left_points)>1:
                firstP = [int((left_points[0][0] + left_points[-1][0])/2),\
                    int((left_points[0][1] + left_points[-1][1])/2)]
            else:
                firstP = left_points[0]
            if len(right_points)>1:
                endP = [int((right_points[0][0] + right_points[-1][0])/2), int((right_points[0][1] + right_points[-1][1])/2)]
            else:
                endP = right_points[0]
            # firstP and endP are respectively the first left point
            # and the first right point of the snippet i. They will be used to
            # draw a line representative of the snippet.

            y_list = np.arange(firstP[0], endP[0],1,dtype = np.int64) #create line (x_list, y_list) between firstP and endP
            slope = (endP[1]-firstP[1])/(endP[0]-firstP[0])
            x_list = np.uint64(slope*(y_list-firstP[0])+firstP[1])      
            
            #condition 2 - on area/length
            area = cv2.contourArea(snippets[i])
            if round(area/y_list.shape[0],1) > 5. and round(area/y_list.shape[0],1) < 9.:
                
                #condition 3 - on the angle with horizontal
                if endP[0]-firstP[0] != 0:
                    angle = m.atan((endP[1]-firstP[1])/(endP[0]-firstP[0]))*180/m.pi
                    if abs(angle)>6 and abs(angle)<45:   
                    #condition 4 - on alignment
                        #count the number of white pixels on this line
                        pix_white = 0
                        for n in range(y_list.shape[0]):
                            USimage[x_list[n], y_list[n],:] = [0,0,255]
                            if I[x_list[n], y_list[n]] > 0 :
                                pix_white = pix_white + 1
                        if pix_white >0 and pix_white/y_list.shape[0] > 0.65:

                            filtered_indices.append(i)
                            line_snip.append([slope, firstP])


    #Look for snippets that are part of the same fascicles
    aligned_snip = [[j] for j in filtered_indices] # contains list of snippets indices part of same fascicle
    for s in range(len(filtered_indices) - 1):
        ylist_s = np.arange(0, I.shape[1],1, dtype = np.int64)
        slope_s = line_snip[s][0] #slope
        firstP_s = line_snip[s][1] #point on the line of the snippet s
        xlist_s = (slope_s*(ylist_s-firstP_s[0])+firstP_s[1] - 2 )
        yx_list = np.uint64(np.vstack((ylist_s, xlist_s))).T
        
        for n in range(s+1, len(filtered_indices)):
            #check if the previous line intersect a snippet n different from s
            intersection = set.intersection(set(map(tuple,yx_list)), set(map(tuple, snippets[filtered_indices[n]][:,0,:])))
            slope_n = line_snip[n][0]
            firstP_n = line_snip[n][1]
            ylist_n = np.arange(0, I.shape[1],1, dtype = np.int64)
            xlist_n = (slope_n*(ylist_n-firstP_n[0])+firstP_n[1] - 2 )
            if intersection and abs(firstP_n[1]-xlist_s[firstP_n[0]])<20 and\
            abs(firstP_s[1]-xlist_n[firstP_s[0]])<20: # if intersection not empty, then s and n are part of the same fascicle
                aligned_snip[s].append(filtered_indices[n])
    
    #example : we obtain aligned_snip = [[1,2,3],[2,5],[3],[4,5],[5],[7],[8,10],[10]]
    # we need to merge [1,2,3], [2,5], [3], [4,5] and [5]; and we need to merge
    #[8,10] with [10].
    # this is what is done below.
    fascicles = [aligned_snip[0]]
    for f in range(1, len(aligned_snip)):
        test = 0
        for ff in range(len(fascicles)):
            if set.intersection(set(aligned_snip[f]), set(fascicles[ff])):
                fascicles[ff] = list(set(fascicles[ff] + aligned_snip[f])) #conversion to set type removes duplicates
                test = 1
        if test == 0:
            fascicles.append(aligned_snip[f])
    print('final merged fascicles', fascicles)
    
    #outputs 
    F = []
    for fa in range(len(fascicles)):
        fasc = snippets[fascicles[fa][0]][:,0,:]
        for s in range(1,len(fascicles[fa])):
            fasc = np.concatenate((fasc, snippets[fascicles[fa][s]][:,0,:]), axis = 0)
        F.append(fasc)
    return F, fascicles