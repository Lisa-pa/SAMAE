""" Automatic and semi-automatic cropping functions """

def autocropping(I, threshCmin, threshCmax, threshRmin, threshRmax, calibV = 0, additionalCrop1 = 0, additionalCrop2 = 0):
    """ Cropping of raw ultrasound image I to get region of interest: Removal 
    of lateral, top and bottom strips by thresholding the mean pixel value 
    of raws and columns.
    If calibV is not zero, then 'additionalCrop1' millimeters are removed from 
    the top of I (for example to remove the skin, additionalCrop1 should be around
    2), and 'additionalCrop2' millimeters are removed from the bottom.
    15 pixels are removed from the right border to make the white vertical
    line disappear.
    
    Args:
        I (array): grayscale or RGB image (the function will convert it to 
        grayscale image if it has not yet been done)
        threshCmin (double): between 0. and 255., minimum threshold value 
        for the lateral strips
        threshCmax (double): between 0. and 255., maximum threshold value 
        for the lateral strips
        threshRmin (double): between 0. and 255., minimum threshold value 
        for the bottom and top strips
        threshRmax (double): between 0. and 255., maximum threshold value 
        for the bottom and top strips
        calibV (double): vertical calibration factor. Default value is 0.
        If calibV is different from 0, then the equivalent of 'additionalCrop1'
        mm is cropped at the top of I to remove skin from the image, and 
        'additionalCrop2' mm are removed from the bottom.
        additionalCrop1, additionalCrop2 (floats): in millimeters, value to
        remove from top and bottom respectively, in addition to the first
        automatic cropping


    Returns:
        I2 (array): same dimension as I, containing only the region of interest

    References:
        Based on [ Automatic detection of skeletal muscle architecture 
        features, Frida Elen Jalborg, Master’s Thesis Spring 2016 ]

    Example:
        > RGB_image = cv2.imread('AponeurosesDetection/aponeurosesdetection/data/simple_echo.jpg', -1)
          #change this example: should use the function simpleimg() to load the image simple_echo, 
          #instead of using cv2.imread with a path
          cropped_img = autocropping(RGB_image, 10., 15., 12., 25.)
          cv2.imshow('image',cropped_img)
          cv2.waitKey(0) & 0xFF
          cv2.destroyAllWindows()
"""
    import cv2
    import numpy as np

    #Check if I is a 1-canal image
    if len(I.shape) == 3:
        grayI = cv2.cvtColor(I, cv2.COLOR_RGB2GRAY)
    elif len(I.shape) == 2:
        grayI = I
    grayI = (grayI*1.05>255)*255. + (grayI*1.05<=255)*grayI*1.1
    left_col = [0]
    right_col = [0]
    up_row = [0]
    bottom_row = [0]
    
    #vertical cropping: find columns
    while threshCmin < threshCmax and \
    (len(right_col) == 1 or len(left_col) == 1):

        for col in range(int(grayI.shape[1]/2)):
            if np.mean(grayI[:, col]) <= threshCmin:
                left_col.append(col)
            if np.mean(grayI[:, grayI.shape[1]-1-col]) <= threshCmin:
                right_col.append(grayI.shape[1]-1-col)
        threshCmin += 1.
    #add 15 pixels cropping in the right of the image (to remove white horizontal line)
    right_col[-1] = right_col[-1] - 15
  
    #horizontal cropping: find rows
    while threshRmin < threshRmax and \
    (len(up_row) == 1 or len(bottom_row) == 1):
        
        for lig in range(int(grayI.shape[0]/2)):
            if np.mean(grayI[lig, :]) <= threshRmin:
                up_row.append(lig)
            if np.mean(grayI[grayI.shape[0]-1-lig, :]) <= threshRmin:
                bottom_row.append(grayI.shape[0]-1-lig)
        threshRmin += 2.

    #cropping image
    if len(I.shape) == 3:
        I2 = I[up_row[-1]+5:bottom_row[-1]-5, left_col[-1]+5:right_col[-1]-5, :]
    elif len(I.shape) == 2:
        I2 = I[up_row[-1]:bottom_row[-1], left_col[-1]:right_col[-1]]
            
    if calibV != 0:
        addcrop1 = int(additionalCrop1 / calibV)
        addcrop2 = int(additionalCrop2 / calibV)
        if addcrop1<I2.shape[0]:
            up_row[-1] = up_row[-1] + addcrop1
            I2 = I2[addcrop1:, :]
        if addcrop2 !=0:
            I2 = I2[: - addcrop2, :]
            bottom_row[-1] = bottom_row[-1] - addcrop2
        
        
    return I2, up_row[-1], bottom_row[-1], left_col[-1], right_col[-1]

def manualcropping(I, pointsfile):
    """This function crops a copy of image I according to points stored 
    in a text file (pointsfile) and corresponding to aponeuroses (see 
    Args section).

    Args:
        I (array): 3-canal image
        pointsfile (text file): contains points' coordinates. Pointsfile must be 
        organized such that:
            - column 0 is the ID of each point
            - column 1 is the X coordinate of each point, ie the corresponding 
            column in I
            - column 2 is the Y coordinate, ie the raw in I
            - raw 0 is for columns' name
            - raws 1 and 2 are for two points of the scale
            - raws 3 to 13 are aponeuroses' points in panoramic images // raws 3 
            to 10 in simple images
            - following raws are for muscle fascicles (and are optional for this 
            function)
            pointsfile's name must 1) include extension 2) indicates whether I
            is panoramic or simple by having 'p' or 's' just before the point
            of the extension.

        Returns:
            I2 (array) : array of same type than I. It is the cropped image of I according
            to the aponeuroses' points manually picked and stored in pointsfile.
            point_of_intersect (tuple) : point which is the most at the right 
            of the image and which should correspond to the point of intersection of deep 
            and upper aponeuroses.
    """

    import numpy as np

    data = open(pointsfile, 'r')
    
    #finds whether the image is panoramic or simple
    search_point = -1
    while (pointsfile[search_point] != '.') and (search_point > (-len(pointsfile))):
        search_point = search_point-1
    if (search_point == -len(pointsfile)):
        raise TypeError("Input pointsfile's name is not correct. Check extension.")
    else:
        imagetype = pointsfile[search_point-1]
    
    #extract points from the input file
    picked_points = []
    for line in data:
        line = line.strip('\n')
        x = line.split('\t')
        picked_points.append((x[1], x[2]))
       
    #keep aponeuroses points according to image type
    if imagetype == 'p': #keep points 3 to 13 included
        apos = np.asarray(picked_points[3:14], dtype=np.float64, order='C')
    elif imagetype == 's': #keep points 3 to 10 included
        apos = np.asarray(picked_points[3:11], dtype=np.float64, order='C')
    else:
        raise ValueError("pointsfile's name does not fulfill conditions. See docstrings")

    #find max and min indexes for columns and raws to crop image I
    #with a margin of 10 pixels (5 pixels for min_raw).
    #Coordinates are inverted in apos
    min_raw = max(0, np.min(apos[:, 1])-10)
    max_raw = min(I.shape[0], np.max(apos[:, 1])+20)
    min_col = max(0, np.min(apos[:, 0])-10)
    max_col = min(I.shape[1], np.max(apos[:, 0])+10)

    i_cropped = np.copy(I[int(min_raw):int(max_raw), int(min_col):int(max_col), :])
    
    index = np.argmax(apos[:, 0])

    point_of_intersect = (apos[index][1] - min_raw, apos[index][0] - min_col)

    #close file
    data.close()
    
    return i_cropped, point_of_intersect, int(min_raw), int(max_raw), int(min_col), int(max_col)