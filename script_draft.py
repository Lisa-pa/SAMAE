'''#what is this script?
from skimage.color import rgb2gray
import cv2

from SAMAE.apoCont import apoContour
import SAMAE.data as dt

# Import sample images 
cirimg = rgb2gray(dt.circle()) # circle sample img
skimg = rgb2gray(dt.skmuscimg()) # ultrasound sample img

# Let's test apoContour
pt1 = [76,85]
contour, n, maxiDPHI, LIF, GIF = apoContour(cirimg, pt1, 1.0, 1.0, 1.0, 65.025, 0.15, 1, 0.01, 3.0, 0.1)

cv2.imshow('Initial I',cirimg)
cv2.imshow('LGIF',contour)
cv2.waitKey(0) & 0xFF
cv2.destroyAllWindows()'''