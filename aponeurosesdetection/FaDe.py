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
    The function is coded to enhance dark tube-like structure, so make
    sure that you do not need to inverse your image I before using this function.
    
    Args:
        I (2D array):       I is a grayscale image (otherwise it is 
                            converted to grayscale)
        thresholds (list):  thresholds is a list of two thresholds that 
                            control the sensitivity of the line filter 
                            to the measures Fr and R. It should look like
                            [b, c], where b must not be zero. If c is zero, 
                            then its value is set to half the frobenius norm
                            of the hessian matrix.
        scales (list):      scales is a list of lengths that correspond to 
                            the diameter (in pixels) of the tube-like structure to find
    
    Outputs:
        I2 (2D array):  one-canal image I, filtered by the multiscale vessel enhancement
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
         ValueError('first element of thresholds cannot be null')
    
    #calculation of Hessian matrix and its eigen values.
    #Corresponding eigen vectors are normal to and in the direction of
    #the tubular structure
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

def locateSnippets(I, xcalib, ycalib, minLength, offSetX = 0, offSetY = 0):
    """
    This function aims at detecting portions of fascicles in a binary image I,
    obtained from the binarization of an image filtered with MVEF.
    These portions of fascicles are filtered with statistical filtering
    on the following parameters: length, alignement, aspect ratio and 
    angle with the horizontal.

    Args:
        I (array): binary image (one canal)
        xcalib (float): vertical calibration factor
        ycalib (float): horizontal calibration factor
        minLength (float):  length in pixels, under which a portion of fascicle
                            will be considered or removed from the analysis.
        offSetX, offSetY (integers):  if the image I is a subimage of another bigger image
                            I2, and if you want to obtain the localization of the
                            snippets in I2, offSetX and offSetY should correspond
                            to the first line of I in I2 and the first column of I
                            in I2 respectively.
    Outputs:
        list of arrays: each array is a snippet's points coordinates
        list of parameters that caracterize the modeling line of a snippet. Each
        parameter is such as [slope, [point1_line + offSetY, pont1_column + offSetX]]
    """
    #attention, contours points found by findContours have the structure (columns, rows)
    snippets = cv2.findContours(np.uint8(I), mode=cv2.RETR_EXTERNAL, method=cv2.CHAIN_APPROX_NONE)[0]
    
    #if no contour is found, stop analysis:
    if len(snippets) == 0:
        return 'error', 'error'
    
    line_snip = [] #list which contains a linear approximation of the snippet sn
    vect = np.zeros((len(snippets), 4))
    
    for i in range(len(snippets)):
        sn = list(snippets[i][:,0,:])
        sn.sort(key=lambda x:x[0]) #sort according to columns  
    
        #######################################################################
        # construction of 2 points p1, p2 to represent sn by a line
        third = 1/3.*len(sn)/2
        left_points = [sn[i] for i in range(len(sn)) if sn[i][0] == sn[0][0] + int(third)]
        right_points = [sn[i] for i in range(len(sn)) if sn[i][0] == sn[-1][0] - int(third)]
        left_points.sort(key=lambda x:x[1]) #sort according to lines
        right_points.sort(key = lambda x:x[1])

        if len(left_points)==2:
            p1 = [int((left_points[0][0] + left_points[-1][0])/2),\
                      int((left_points[0][1] + left_points[-1][1])/2)]
        elif len(left_points) == 1:
            p1 = left_points[0]
        elif len(left_points) > 2: #take the extremal left point of the snippet
            p1 = [sn[i] for i in range(len(sn))\
                      if sn[i][0] == np.amin(np.array(sn)[:,0])][0]
        
        if len(right_points)==2:
            p2 = [int((right_points[0][0] + right_points[-1][0])/2),\
                    int((right_points[0][1] + right_points[-1][1])/2)]
        elif len(right_points)==1:
            p2 = right_points[0]
        elif len(right_points) > 2: #take right extremal point of the snippet
            p2 = [sn[i] for i in range(len(sn))\
                    if sn[i][0] == np.amax(np.array(sn)[:,0])][0]

        # create line (x_list, y_list) inside snippet's contour
        cmin = np.amin(snippets[i][:,0,0])
        cmax = np.amax(snippets[i][:,0,0])
        y_list = np.arange(cmin, cmax, 1)
        if p2[0]-p1[0] != 0:
            slope = (p2[1]-p1[1])/(p2[0]-p1[0])
            x_list = (slope*(y_list-p1[0])+p1[1])
            line_snip.append([slope, [p1[0] + offSetY, p1[1] + offSetX]])
        else:
            slope = 0
            x_list = np.arange(np.amin(snippets[i][:,0,1]), np.amax(snippets[i][:,0,1])+1, 1)
            line_snip.append([slope, [p1[0] + offSetY, p1[1] + offSetX]])
        
        
        #######################################################################
        
        # Features: length, alignement, aspect ratio and angle with horizontal
        if len(sn) > 5 : #make sure minimum of points required for fitEllipse function is reached
            (x0,y0),(mia,maa), ang = cv2.fitEllipse(snippets[i])
            angle = 90 - ang
            aspra = float(mia/maa)
            if angle == 90:
                align = 0
            else:
                pix_white = 0
                for n in range(y_list.shape[0]):
                    if int(x_list[n])<I.shape[0]:
                        if I[int(x_list[n]), int(y_list[n])] > 0 :
                            pix_white = pix_white + 1
                align = pix_white/y_list.shape[0]           
        else:
            maa = 0
            angle = 0
            aspra = 0.
            align = 0
        
        vect[i,0] = maa * xcalib
        vect[i,1] = angle
        vect[i,2] = align
        vect[i,3] = aspra
        
    # Normalization
    medianL = np.median(vect[:,0])
    medianAng = np.median(vect[:,1])
    medianAl = np.median(vect[:,2])
    medianAR = np.median(vect[:,3])
    
    madL = np.median(abs(vect[:,0] - medianL))
    madAng = np.median(abs(vect[:,1] - medianAng))
    madAl = np.median(abs(vect[:,2] - medianAl))
    madAR = np.median(abs(vect[:,3] - medianAR))
    
    vect[:,0] = (vect[:,0] - medianL) / madL
    vect[:,1] = (vect[:,1] - medianAng) / madAng
    vect[:,2] = (vect[:,2] - medianAl) / madAl
    vect[:,3] = (vect[:,3] - medianAR) / madAR
        
    norm = np.sqrt(vect[:,0]*vect[:,0] + vect[:,1]*vect[:,1]\
        + vect[:,2]*vect[:,2] + vect[:,3]*vect[:,3])
   
    #Filtering on norm
    #The value "4" has been empirically chosen by observing snippets depending on their norm value
    lines = []
    filtered_indices = [] # list which contains the ID of the filtered snippets
    for n in range(len(snippets)):
        if norm[n] >= 4 and vect[n,0] > minLength:
            filtered_indices.append(n)
            lines.append(line_snip[n])
    
    #outputs
    #retrieve snippets' points coordinates if their ID is in filtered_indices
    filtered_snippets = []
    for fa in range(len(filtered_indices)):
        fasc = snippets[filtered_indices[fa]][:,0,:]
        fasc[:, [0,1]] = fasc[:,[1,0]] #change structure to (rows, columns)
        fasc[:, 0] = fasc[:, 0] + offSetX
        fasc[:,1] = fasc[:, 1] + offSetY
        filtered_snippets.append(fasc)
    
    return filtered_snippets, lines




def combineSnippets(I, listS, listS_paramlines, min_nb_sn, thresh_alignment = 20):
    """
    This function aims at combining aligned snippets (we consider that, if two
    snippets are aligned, they correspond to the same muscle fascicle) so that to
    reconstruct muscle fascicles.


    Args:
        I : one-canal or three-canal image
        listS (list of arrays): list of arrays, each array contains one snippet's points coordinates (each
        point is given by (line, column))
        listS_paramlines (list of lists): list of sublist. Each sublists contains parameters that 
        caracterize the modeling line of a snippet. The format of this list is
        the same output from the function locateSnippets. The format of sublists:
        [slope, [point1_line + offSetY, pont1_column + offSetX]]
        min_nb_sn (integer): once aligned, you have the possibility to differentiate
        the fascicles with a number of detected snippets below min_nb_sn, and above
        min_nb_sn. This allows you to process those fascicles differently in the 
        following of your script (for exemple for the degree of modeling curves).
        thresh_alignment (interger): this is the threshold that determines if two snippets are aligned.

                Potential aligned snippets are such that the modeling line of one intercept the contour
                of the other.
                thresh_alignment corresponds to the maximum number of lines accepted between:
                        - point1 of the most right snippet
                        - the point with same column as point1, on the modeling line of the most left snippet
                In addition, the two snippets must be less than 1/4 of the image length far,
                to be considered as aligned snippets.
    
    Outputs:
        two lists of arrays. Each array contains the points coordinates of the
        aligned snippets of a same fascicle.
        The first list contains the fascicles with more than min_nb_sn snippets.
        The second list contains the fascicles with less than min_nb_sn snippets.
    """
    fascicles_points = []
    fascicles2_points = []    
    if len(listS) != 0:
        aligned_snip = [[j] for j in range(len(listS))] # contains list of snippets indices part of same fascicle
    
        #modeling line of snippet s
        for s in range(len(listS) - 1):
            ylist_s = np.arange(0, I.shape[1],1, dtype = np.int64)
            slope_s = listS_paramlines[s][0] #slope
            p1_s = listS_paramlines[s][1] #point on the line of the snippet s = [y, x]
            xlist_s = (slope_s*(ylist_s-p1_s[0])+p1_s[1])
            xy_list_s = np.uint64(np.vstack((xlist_s, ylist_s))).T

            for n in range(s+1, len(listS)):
                #check if the previous line intersect a snippet n different from s
                
                slope_n = listS_paramlines[n][0]
                p1_n = listS_paramlines[n][1]
                ylist_n = np.arange(0, I.shape[1],1, dtype = np.int64)
                xlist_n = (slope_n*(ylist_n-p1_n[0])+p1_n[1])
                xy_list_n = np.uint64(np.vstack((xlist_n, ylist_n))).T
                
                intersection1 = set.intersection(set(map(tuple,xy_list_s)), set(map(tuple, listS[n])))
                intersection2 = set.intersection(set(map(tuple,xy_list_n)), set(map(tuple, listS[s])))
                
                
                #if intersection not empty  
                if (intersection1 or intersection2):

                    if np.amin(np.array(listS[n])[:,1]) < np.amin(np.array(listS[s])[:,1]):
                        # if snippets roughly aligned :
                        if abs(p1_s[1]-xlist_n[p1_s[0]])<thresh_alignment:
                            # if snippets not too far from each other 
                            # ie maximum a quarter of I's length
                            if abs(np.amax(np.array(listS[n])[:,1])\
                               - np.amin(np.array(listS[s])[:,1])) < I.shape[1]/4:
                                aligned_snip[s].append(n)
                                                               
                    elif np.amin(np.array(listS[n])[:,1]) >= np.amin(np.array(listS[s])[:,1]):
                        # if snippets roughly aligned :
                        if abs(p1_n[1]-xlist_s[p1_n[0]])<thresh_alignment:
                            # if snippets not too far from each other 
                            # ie maximum a quarter of I's length                
                            if abs(np.amin(np.array(listS[n])[:,1])\
                                   - np.amax(np.array(listS[s])[:,1])) < I.shape[1]/4:
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
                grouped_snip.append(aligned_snip[f])
                    
        #outputs 
        for fa in range(len(grouped_snip)):
            fasc = listS[grouped_snip[fa][0]]
            if len(grouped_snip[fa]) < min_nb_sn:
                fascicles2_points.append(fasc)
            else:
                for s in range(1,len(grouped_snip[fa])):
                    fasc = np.concatenate((fasc, listS[grouped_snip[fa][s]]), axis = 0)
                fascicles_points.append(fasc)
    
    return fascicles_points, fascicles2_points




def contourAverage(inputC):
    """Function that approximate a contour by a line inside the contour.
        It realizes the mean of the lines of two points on a same column.
        If there is one point or strictly more than 2 points in a same
        column, the column is not considered and no mean is realized.
    
    Args:
        inputC (list of arrays): each array corresponds to the points coordinates
        of a fascicle. The format of inputC is the same as one of the output lists 
        of function combineSnippets. 
        Points coordinates are (rows along axis 0, columns along axis 1).

    Outputs:
        ouputL : list of arrays that contain points [row, column]
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


def approximateFasc(typeapprox, listF, d):
    """
    Realizes the interpolation of a list of fascicles, either with polynomial fitting
    or with B-splines.


    Inputs:
        listF is the list of fascicles (same format as output from function contourAverage).
        Each fascicle is an array of points defining its line. Each point = (row, column).
        d (integer) : degree of the curve that will be interpolated
        typeapprox (string): either 'Bspline' or 'polyfit', depending on the type of 
        approximation that you want to realize (spline or fitting with a polynom respectively) 
    Outputs:
        a list of splines.
    """
    approx_fasc = []
    
    if typeapprox == 'Bspline':
        import scipy.interpolate as interpolate
        
        #transform each fascicle format in list
        for n in range(len(listF)):
            if type(listF[n]) == list:
                f = listF[n]
            else:
                f = list(listF[n])
            
            #sort all fascicle points in increasing columns order
            f.sort(key = lambda x: x[1])

            #remove potential double points           
            to_remove = []
            for pt in range(len(f)-1):
                if f[pt][1] == f[pt+1][1]:
                    to_remove.append(pt+1)
            fa = [f[ind] for ind in range(len(f)) if ind not in to_remove]
            
            ycoord = [fa[i][1] for i in range(len(fa))]
            xcoord = [fa[i][0] for i in range(len(fa))]               
            
            #interpolation
            spline = interpolate.UnivariateSpline(ycoord, xcoord, k=d, ext = 0)   
            approx_fasc.append(spline)

    if typeapprox == 'polyfit':
        for n in range(len(listF)):

            #change format to list if necessary
            if type(listF[n]) == list:
                f = listF[n]
            else:
                f = list(listF[n])
                  
            #sort fascicle points in increasing columns order
            f.sort(key = lambda x: x[1])
            
            #remove potential double points
            to_remove = []
            for pt in range(len(f)-1):
                if f[pt][1] == f[pt+1][1]:
                    to_remove.append(pt+1)
            fa = [f[ind] for ind in range(len(f)) if ind not in to_remove]

            ycoord = [fa[i][1] for i in range(len(fa))]
            xcoord = [fa[i][0] for i in range(len(fa))]               
            
            #interpolation and conversion to spline
            p = np.polyfit(ycoord, xcoord, deg=d)
            spline = np.poly1d(p)
            approx_fasc.append(spline)
        
    return approx_fasc