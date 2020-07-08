  
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
    """Download an image from a Dropbox account.

    Args:
        tok (string): access token that connects to the wanted
                      app in Dropbox account
        path2file (string): Path of the file to download, in the
                            app corresponding to the above token.

    Output:
            image (numpy.ndarray): 3-channel color image, with 
                                    coefficients' type == uint8

    Example:
            1) Register a new app in the App Console of your Dropbox
            account. Set up parameters as you want.
            2) In Dropbox>Applications>MyApp, import your data.
            3) In the settings page of MyApp, generate a token and copy it.
                It should look like a random string of letters and figures,
                as below. (!!!This access token can be used to access your
                account via the API. Don’t share your access token with anyone!!!) 
        > token = 'Q8yhHQ4wquAAAAAAAAABRPb9LYdKAr2WGcmhhJ8egiX4_Qak6YZwBw4GUpX9DVeb'
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