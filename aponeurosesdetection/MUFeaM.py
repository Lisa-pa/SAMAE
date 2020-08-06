def muscleThickness(I, points1, points2, start, end, calibV, calibH):
    """ Function that calculates muscle thickness  in mm on the 
    interval of columns [start, end] of a US muscle image I.
    
    Args:
        I (array): three-canal image
        points1 (array): dimensions = (n,2)
                        contains the interpolated points of one 
                        aponeurosis in I.
        points2 (array): dimensions = (n,2)
                        contains the interpolated points of the 
                        other aponeurosis in I.
        start (int) : starting column to calculate muscle thickness. Its
                        value should be >= 0 and <I.shape[1]
        end (int): ending column to calculate muscle thickness. Its value
                    should be > start, > 0 and <I.shape[1]
        calibV (float) : vertical calibration factor
        calibH (float) : horizontal calibration factor
    
    Outputs:
            absc (list) : list of abscissa representing distances the beginning
                        of the image I
            mt (list) : muscle thickness in mm at the respecting abscissa
            spl (tuple): spline that interpolates the previous points (x,y)
    """
    if start < 0 or end < 0  or end < start\
        or start >= I.shape[1] or end >= I.shape[1]:
        raise ValueError("start or end values go out of bounds of I's dimensions")

    mt = []
    absc = []
    for col in range(start, end + 1):

        #search if there exist points in each aponeurosis with abscissa col
        search1 = [pt for pt in points1 if pt[1]==col]
        search2 = [pt for pt in points2 if pt[1]==col]

        if search1 and search2: #if there exist a point whith abscissa col 
                                #in each list points1, points2, then calculates MT
            mt.append(abs(search1[0][0]-search2[0][0])*calibV)
            absc.append(col*calibH)
    
    import scipy.interpolate as interpolate
    spl = interpolate.UnivariateSpline(absc, mt, k=5, ext = 0) 
    return absc, mt, spl