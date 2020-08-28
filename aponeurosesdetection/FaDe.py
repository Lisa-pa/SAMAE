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

def locateSnippets(I, xcalib, ycalib, minLength, rangeAngles, percentageAlign, offSetX = 0, offSetY = 0):
    """
    I (array): binary image (one canal)
    """
    #attention, contours points found by findContours have the structure (columns, rows)
    snippets = cv2.findContours(np.uint8(I), mode=cv2.RETR_EXTERNAL, method=cv2.CHAIN_APPROX_NONE)[0]
    
    if len(snippets) == 0:
        raise ValueError('No fascicle could be found in I.')
    
    filtered_indices = [] # list which contains the ID of the filtered snippets
    line_snip = [] #list which contains a linear approximation of the snippet filtered
    
    for i in range(len(snippets)):
        
        #condition 1 - on length :
        if snippets[i].shape[0]/2 > minLength/m.sqrt(xcalib**2+ycalib**2):
            
            
            # creation of a line representative of the snippet i
            # firstP and endP are are respectively
            # in the left part and in the right part of snippet i
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
                endP = [int((right_points[0][0] + right_points[-1][0])/2),\
                        int((right_points[0][1] + right_points[-1][1])/2)]
            else:
                endP = right_points[0]
                
            y_list = np.arange(firstP[0], endP[0], 1, dtype = np.int64) #create line (x_list, y_list) between firstP and endP


            #condition 2 - on area/length
                ###this condition is too restrictive
            """area = cv2.contourArea(snippets[i])
            if round(area/y_list.shape[0],1) > 5. and round(area/y_list.shape[0],1) < 9.:
            """   
            
            #condition 3 - on the angle with horizontal
            if endP[0]-firstP[0] != 0:
                
                slope = (endP[1]-firstP[1])/(endP[0]-firstP[0])
                angle = -m.atan(slope)*180/m.pi            

                if angle>rangeAngles[0] and abs(angle)<rangeAngles[1]:
                    
                    x_list = np.uint64(slope*(y_list-firstP[0])+firstP[1])  
                    

                #condition 4 - on alignment
                    #count the number of white pixels on this line
                    pix_white = 0
                    for n in range(y_list.shape[0]):
                        if I[x_list[n], y_list[n]] > 0 :
                            pix_white = pix_white + 1
                    if pix_white/y_list.shape[0] > percentageAlign: #0.65 for simple images, 0.8 for pano

                        filtered_indices.append(i)

                        line_snip.append([slope, [firstP[0] + offSetY, firstP[1] + offSetX]])
                            
    #outputs 
    filtered_snippets = []
    for fa in range(len(filtered_indices)):
        fasc = snippets[filtered_indices[fa]][:,0,:]
        fasc[:, [0,1]] = fasc[:,[1,0]] #change structure to (rows, columns)
        fasc[:, 0] = fasc[:, 0] + offSetX
        fasc[:,1] = fasc[:, 1] + offSetY
        filtered_snippets.append(fasc)
    
    return filtered_snippets, line_snip






def combineSnippets(I, listS, listS_paramlines, thresh_alignment = 20, thresh_length = 200):
    """listS contient les snippets en mode [x,y]"""
    #Look for snippets that are part of the same fascicles

    fascicles_points = []

    
    if len(listS) != 0:
        aligned_snip = [[j] for j in range(len(listS))] # contains list of snippets indices part of same fascicle
    
        for s in range(len(listS) - 1):
            ylist_s = np.arange(0, I.shape[1],1, dtype = np.int64)
            slope_s = listS_paramlines[s][0] #slope
            firstP_s = listS_paramlines[s][1] #point on the line of the snippet s
            xlist_s = (slope_s*(ylist_s-firstP_s[0])+firstP_s[1])
            xy_list_s = np.uint64(np.vstack((xlist_s, ylist_s))).T
                        
            for n in range(s+1, len(listS)):
                #check if the previous line intersect a snippet n different from s
                
                slope_n = listS_paramlines[n][0]
                firstP_n = listS_paramlines[n][1]
                ylist_n = np.arange(0, I.shape[1],1, dtype = np.int64)
                xlist_n = (slope_n*(ylist_n-firstP_n[0])+firstP_n[1])
                xy_list_n = np.uint64(np.vstack((xlist_n, ylist_n))).T
                
                intersection1 = set.intersection(set(map(tuple,xy_list_s)), set(map(tuple, listS[n])))
                intersection2 = set.intersection(set(map(tuple,xy_list_n)), set(map(tuple, listS[s])))
                
                #if intersection not empty and if first points of snippets roughly aligned : 
                if (intersection1 or intersection2) and\
                abs(firstP_n[1]-xlist_s[firstP_n[0]])<thresh_alignment and\
                abs(firstP_s[1]-xlist_n[firstP_s[0]])<thresh_alignment:
                    
                    # if snippets not too far from each other 
                    # ie maximum a quarter of I's length
                    if np.amin(np.array(listS[n])[:,1]) < np.amin(np.array(listS[s])[:,1]): 
                        if abs(np.amax(np.array(listS[n])[:,1])\
                               - np.amin(np.array(listS[s])[:,1])) < I.shape[1]/5:
                            aligned_snip[s].append(n)
                                                               
                    elif np.amin(np.array(listS[n])[:,1]) >= np.amin(np.array(listS[s])[:,1]):
                        if abs(np.amin(np.array(listS[n])[:,1])\
                               - np.amax(np.array(listS[s])[:,1])) < I.shape[1]/5:
                            aligned_snip[s].append(n)


        #example : we obtain aligned_snip = [[1,2,3],[2,5],[3],[4,5],[5],[7],[8,10],[10]]
        # we need to merge [1,2,3], [2,5], [3], [4,5] and [5]; and we need to merge
        #[8,10] with [10].
        # this is what is done below.
        grouped_snip = [aligned_snip[0]]
        for f in range(1, len(aligned_snip)):
            test = 0
            
            for ff in range(len(grouped_snip)):
                if set.intersection(set(aligned_snip[f]), set(grouped_snip[ff])):
                    grouped_snip[ff] = list(set(grouped_snip[ff] + aligned_snip[f])) #conversion to set type removes duplicates
                    test = 1
            
            if test == 0:
                #if more than 2 snippets in aligned_snip[f], consider it in grouped_snip:
                
                if len(aligned_snip[f]) >= 2:
                    grouped_snip.append(aligned_snip[f])
                    
                # if there is only 1 snippet in aligned_snip[f], check if it
                # is long enough
                
                if len(aligned_snip[f]) == 1 and len(listS[aligned_snip[f][0]]) > thresh_length:
                    grouped_snip.append(aligned_snip[f])
                
        print('final merged fascicles', grouped_snip)
        
        #outputs 
        for fa in range(len(grouped_snip)):
            fasc = listS[grouped_snip[fa][0]]
            for s in range(1,len(grouped_snip[fa])):
                fasc = np.concatenate((fasc, listS[grouped_snip[fa][s]]), axis = 0)
            fascicles_points.append(fasc)
    """
    couleurs = [[255,0,0], [0,255,0], [0,0,255], [255,255,0],[255,0,255], [0,255,255],\
                [255,255,255], [100,200,0],[100,200,100], [50,200,0],[50,100,50], [255,100,0],\
                [120,120,255], [255,80,80],[0,100,200], [0,100,80]]
    for t in range(len(fascicles_points)):
        for k in range(fascicles_points[t].shape[0]):
            pt = fascicles_points[t][k,:]
            I[pt[0], pt[1], :] = couleurs[t]
    """

    return fascicles_points




def contourAverage(inputC):
    """Function that approximate a contour by a line
    inputC's dimensions are (x,y); x coordinate of the image along axis 0
    ( = rows), y along axis 1 (= columns).

    ouputL : list of points [line, column]
    """
    input_list = list(inputC)
    
    outputL = []
    
    min_column = np.amin(inputC[:,1])
    max_column = np.amax(inputC[:,1])
    
    for c in range(min_column, max_column + 1):
        pt = [point for point in input_list if point[1] == c]
        if len(pt) == 2:
            outputL.append(np.uint64(np.mean(pt, axis = 0)))
    return outputL


def approximateFasc(I, typeapprox, listF, d):
    """listF is the list of fascicles.
    Each fascicle is a list of point defining its line
    
    points of each fascicle : array = [x y] where x is the line, y the column
    d : degree of the curve that will be interpolated
    """
    approx_fasc = []
    
    if typeapprox == 'Bspline':
        import scipy.interpolate as interpolate
        for n in range(len(listF)):
            if type(listF[n]) == list:
                f = listF[n]
            else:
                f = list(listF[n])
                
            f.sort(key = lambda x: x[1]) #sort according to columns
            ycoord = [f[i][1] for i in range(len(f))]
            xcoord = [f[i][0] for i in range(len(f))]               
            
            spline = interpolate.UnivariateSpline(ycoord, xcoord, k=d, ext = 0)   
            approx_fasc.append(spline)

    if typeapprox == 'polyfit':
        for n in range(len(listF)):
            if type(listF[n]) == list:
                f = listF[n]
            else:
                f = list(listF[n])
                  
            f.sort(key = lambda x: x[1]) #sort according to columns
            ycoord = [f[i][1] for i in range(len(f))]
            xcoord = [f[i][0] for i in range(len(f))]               
            
            p = np.polyfit(ycoord, xcoord, deg=d)
            spline = np.poly1d(p)
            approx_fasc.append(spline)
        
    return approx_fasc
