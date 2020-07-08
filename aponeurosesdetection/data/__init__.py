  
"""Standard test images.


"""

import os
from skimage.io import imread

data_dir = os.path.abspath(os.path.dirname(__file__))


__all__ = ['data_dir', 'circle', 'skmuscimg']


def _load(f, as_gray=False):
    """Load an image file located in the data directory.
    Parameters
    ----------
    f : string
        File name.
    as_gray : bool, optional
        Whether to convert the image to grayscale.
    Returns
    -------
    img : ndarray
        Image loaded from ``data_dir``.
    """
    # importing io is quite slow since it scans all the backends
    # we lazy import it here
    return imread(f, as_gray=as_gray)


def circle():
    """Synthetic image of a circle
    Returns
    -------
    circle : (xdim, ydim) bool ndarray
        Circle image.
    """
    return _load(os.path.join(data_dir, "circle.bmp"))


def skmuscimg():
    """[summary]
    """

    return _load(os.path.join(data_dir, "skmuscle.jpg"))


def downloadFromDropbox(tok, path2file):
    """Downloads a file from a dropbox account

    Args:
        tok (string): access token that connects to the wanted
                      app in Dropbox account
        path2file (string): Path of the file to download, in the
                            app corresponding to the upward token.

    Output:
            image (numpy.ndarray): 3-channel color image, with 
                                    coefficients' type == uint8

    Example (not working yet, I have to find a way to generate tokens that
    to not give total access to my dropbox account but only to the desired
    file):
        > token = 'randomstring_887_85sf84654f5dfdfggkf'
        > path = '/cropped_20181002_153426_image.jpg'
        > dt = downloadFromDropbox(token, path);
    """
    import dropbox
    import numpy as np
    import cv2

    dbx  = dropbox.Dropbox(tok)
    try:
        metadata, file = dbx.files_download(path2file)
    except dropbox.exceptions.HttpError as err:
        print('*** HTTP error', err)
        return None
       
    data = np.frombuffer(file.content, np.uint8)
    image = cv2.imdecode(data, 1) # cv2.IMREAD_COLOR in OpenCV 3.1

    return image