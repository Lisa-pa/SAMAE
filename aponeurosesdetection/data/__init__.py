  
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
