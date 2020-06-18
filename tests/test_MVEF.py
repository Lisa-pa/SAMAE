def test_MVEF(chemin):
    """Test Multiscale Vessel Enhancement Method on an image located at chemin."""
     
    'Opening of the image and color transform.'
    RGB_image = cv2.imread(chemin, -1);
    gray_I = cv2.cvtColor(RGB_image, cv2.COLOR_RGB2GRAY);
    gray_I2 = 255 - gray_I; #inverse black and white

    'Multiscale Vessel Enhancement Method'
    seg, vesselness,hessXY, hessXX, hessYY = MVEF_2D(gray_I2, [4.], [0.5, 0.5]);

    'Rescale the pixels values if max. value is too low'
    maxi = np.max(seg);
    if maxi <=0.5 and maxi>0:
        seg = (seg/maxi)*255;

    'Visualization'
    cv2.imshow('MVEF',seg);
    cv2.imshow('Initial image', gray_I);
    cv2.imshow('Inverse image', gray_I2);
    cv2.waitKey(0) & 0xFF;
    cv2.destroyAllWindows();