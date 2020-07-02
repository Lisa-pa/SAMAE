"*****************************************************************************"
"**********************************PACKAGES***********************************"
"*****************************************************************************"
import cv2;
import numpy as np;
import math as m
"*****************************************************************************"
"**********************************FUNCTIONS**********************************"
"*****************************************************************************"

'-----------------------------------------------------------------------------'
def heaviside(matrix, epsilon):
    h = np.zeros((matrix.shape));
    for u in range(matrix.shape[0]):
        h[u,:] = np.array([1/2.*(1+2.*m.atan(matrix[u,v]/epsilon)/m.pi) for v in range(matrix.shape[1])]);
    return h
'-----------------------------------------------------------------------------'
def gaussianKernel(I, point, sigma):
    absc = point[0];
    ordo = point[1];
    k = np.zeros((I.shape));
    for a in range(I.shape[0]):
        k[a,:] = np.array([1/(2*m.pi*sigma**2) * np.exp(-((a-absc)**2 + (b-ordo)**2)/(2*sigma*sigma))for b in range(I.shape[1])]);
    return k
'-----------------------------------------------------------------------------'
def intensities(I, previousC, sigma, epsilon):

    #local intensities
    heavisideMatrix = heaviside(previousC, epsilon);        
    OneMinusHeaviside = np.ones((I.shape[0],I.shape[1])) - heavisideMatrix;
    HI = heavisideMatrix*I;
    OneMinusH_I = OneMinusHeaviside*I;

    convolHI = cv2.GaussianBlur(HI, ksize = (2*int(sigma)+1,2*int(sigma)+1), sigmaX = sigma, sigmaY =sigma, borderType=cv2.BORDER_CONSTANT);
    convolH = cv2.GaussianBlur(heavisideMatrix, ksize = (2*int(sigma)+1,2*int(sigma)+1),sigmaX = sigma,sigmaY = sigma, borderType=cv2.BORDER_CONSTANT);
    convol1minusH_I = cv2.GaussianBlur(OneMinusH_I, ksize = (2*int(sigma)+1,2*int(sigma)+1),sigmaX =sigma, sigmaY =sigma, borderType=cv2.BORDER_CONSTANT);
    convol1minusH = cv2.GaussianBlur(OneMinusHeaviside, ksize = (2*int(sigma)+1,2*int(sigma)+1), sigmaX =sigma, sigmaY =sigma, borderType=cv2.BORDER_CONSTANT);
    
    f1 = convolHI/convolH;
    f2 = convol1minusH_I/convol1minusH;

    #global intensity
    c1 = np.sum(HI)/np.sum(heavisideMatrix);
    c2 = np.sum(OneMinusH_I)/np.sum(heavisideMatrix);
           
    return c1, c2, f1, f2
'-----------------------------------------------------------------------------'
def apoContour(I, pointIni, l1, l2, mu, nu, dt, epsilon, omega, sigma, stop_thresh):
    
    if len(I.shape) > 2:
        I = cv2.cvtColor(I, cv2.COLOR_RGB2GRAY);
        
    'Initialization step'
    previousC = np.zeros((I.shape[0], I.shape[1])); #circle, center = point, diameter = 4
    for i in range(I.shape[0]):
        for j in range(I.shape[1]):
            previousC[i,j] = -2 + m.sqrt((i-pointIni[0])**2 + (j - pointIni[1])**2);
    test = 1000.;
    step = 0;
    
    'Recurrence'
    #I added a condition in while loop to limit the number of iterations
    while test> stop_thresh and step<=16:
        dirac = np.zeros((I.shape[0], I.shape[1]));

        #GLOBAL AND LOCAL INTENSITIES at current step

        c1, c2, f1, f2  = intensities(I, previousC, sigma, epsilon);
        GIF_Force = omega*(-l1*(I-c1*np.ones(I.shape))*(I-c1*np.ones(I.shape))\
                       + l2 * (I - c2*np.ones(I.shape))*(I - c2*np.ones(I.shape)));
        LIF_Force = np.zeros(I.shape);
        
        #in the following names, P means Plus and M means Minus
        prevC_iP1 = np.concatenate((previousC[1:,:],np.zeros((1,previousC.shape[1]))), axis = 0);
        prevC_iM1 = np.concatenate((np.zeros((1,previousC.shape[1])),previousC[:previousC.shape[0]-1,:]), axis = 0)
        prevC_jP1 = np.concatenate((previousC[:,1:],np.zeros((previousC.shape[0],1))), axis = 1);
        prevC_jM1 = np.concatenate((np.zeros((previousC.shape[0],1)),previousC[:,:previousC.shape[1]-1]), axis = 1);
        prevC_ijM1 = np.concatenate((np.zeros((previousC.shape[0],1)), np.concatenate((np.zeros((1,previousC.shape[1]-1)), previousC[:-1,:-1]),axis = 0)),axis = 1);
        prevC_ijP1 = np.concatenate((np.concatenate((previousC[1:,1:],np.zeros((1,previousC.shape[1]-1))),axis = 0), np.zeros((previousC.shape[0],1))),axis = 1);
        prevC_iP1jM1 = np.concatenate((np.zeros((previousC.shape[0],1)), np.concatenate((previousC[1:,:-1], np.zeros((1,previousC.shape[1]-1))),axis = 0)),axis = 1);
        prevC_iM1jP1 = np.concatenate((np.concatenate((np.zeros((1,previousC.shape[1]-1)),previousC[:-1,1:]),axis = 0), np.zeros((previousC.shape[0],1))),axis = 1);

        norm = np.ones(previousC.shape)/np.sqrt(\
                      (prevC_iP1-previousC)*(prevC_iP1-previousC)\
                      +(prevC_jP1-previousC)*(prevC_jP1-previousC));
        
        divergence = (np.ones(previousC.shape)+1/4*(prevC_iP1-prevC_iM1)*(prevC_iP1-prevC_iM1))\
        *(prevC_iP1+prevC_iM1-2*previousC)\
        +(np.ones(previousC.shape)+1/4*(prevC_jP1-prevC_jM1)*(prevC_jP1-prevC_jM1))\
        *(prevC_jP1+prevC_jM1-2*previousC)\
        +1/8*(prevC_iP1-prevC_iM1)*(prevC_jP1-prevC_jM1)\
        *(prevC_ijP1-prevC_iP1jM1-prevC_iM1jP1-prevC_ijM1);
        
        div = divergence / norm;
        
        lap = prevC_iP1 + prevC_iM1 + prevC_jP1 + prevC_jM1\
        -4*previousC;
        
        for i in range(previousC.shape[0]):
            
            #to follow what is going on
            print('step=',step);
            percent1 = i/previousC.shape[0]*100;
            print('i=', i, percent1);
                
            for j in range(previousC.shape[1]):             
                   
                K_ij = gaussianKernel(I, [i,j], sigma);
                I_f1_square_ij = ((I[i,j]*np.ones(I.shape))-f1)*((I[i,j]*np.ones(I.shape))-f1);
                I_f2_square_ij = ((I[i,j]*np.ones(I.shape))-f2)*((I[i,j]*np.ones(I.shape))-f2);
                LIF_Force[i,j] = -l1*np.sum((K_ij*I_f1_square_ij))+l2*np.sum((K_ij*I_f2_square_ij));
                LIF_Force[i,j] = (1-omega)*LIF_Force[i,j];
                
                dirac[i,j] = epsilon / (m.pi*(epsilon**2+previousC[i,j]**2));
        
        d_phi = dirac*((1 - omega)*LIF_Force+omega*GIF_Force)\
        +nu*dirac*div + mu*(lap - div);
        
        #what happens at borders ? temporary solution-> enforce no variations
        d_phi[0:1,:]=np.zeros((1,d_phi.shape[1]));
        d_phi[d_phi.shape[0]-1:d_phi.shape[0],:]=np.zeros((1,d_phi.shape[1]));
        d_phi[:,0:1]=np.zeros((d_phi.shape[0],1));
        d_phi[:,d_phi.shape[1]-1:d_phi.shape[1]]=np.zeros((d_phi.shape[0],1));
        
        #calculation of new value 
        newC = d_phi*dt + previousC;
                         
        #condition re-calculated for while loop
        test = np.max(np.absolute(d_phi*dt)); #is there a better test ?
        print('test=',test);
        
        #replace previous contour by new one
        previousC = newC;
        step = step+1; 
    
    return previousC,step, test, LIF_Force, GIF_Force

"*****************************************************************************"
"************************************TESTS************************************"
"*****************************************************************************"
image = cv2.imread('Ini2.bmp', -1);
imageG = cv2.cvtColor(image, cv2.COLOR_RGB2GRAY);  
pt1 = [76,85];
contour, n, maxiDPHI, LIF, GIF = apoContour(image, pt1, 1.0, 1.0, 1.0, 65.025, 0.15, 1, 0.01, 3.0, 0.1);
cv2.imshow('Initial I',image);
cv2.imshow('LGIF',contour);
cv2.waitKey(0) & 0xFF;
cv2.destroyAllWindows();