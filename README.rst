====================
AponeurosesDetection
====================


.. image:: https://img.shields.io/pypi/v/aponeurosesdetection.svg
        :target: https://pypi.python.org/pypi/aponeurosesdetection

.. image:: https://img.shields.io/travis/Lisa-pa/aponeurosesdetection.svg
        :target: https://travis-ci.com/Lisa-pa/aponeurosesdetection

.. image:: https://readthedocs.org/projects/aponeurosesdetection/badge/?version=latest
        :target: https://aponeurosesdetection.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status




This package contains functions that will enable automatic aponeuroses and tubular structures detection in ultrasound images.


* Free software: MIT license
* Documentation: https://aponeurosesdetection.readthedocs.io.


Features
--------
* Structure of the folder containing your US images to analyze
        TO DO

* How to use this package?
        To analyse US images with this package, open the filemanager.py file and adapt the 
                - path_to_folder variable (path to the folder that contains the images);
                - part variable (list of participants);
                - colors variable (list of colors for the visualization);
        then run the file !

* What happens when you run the filemanager.py file?

        - Simple images
                1) The image appears with a window asking you to validate the start of the analysis
                2) Scale is automatically detected
                3) The image is automatically cropped and a window asks you to validate the cropping
                        Initial thresholds are set to (10, 15, 12, 25, 2, 6). In some cases, no image
                        is output, so a second automatic try is launched with thresholds (6, 15, 6, 25, 0, 0).
                        (to check which thresholds have been used, have a look at the command window).
                         If the cropping is not validated, new thresholds are asked (maximum five tries for new thresholds are asked per image).
                          - Case a) The cropping is too aggressive (often happens in case of
                           darker images). Try thresholds (6, 15, 6, 25, 0, 0) (that is, reduce 1st and 3rd thresholds).
                          - Case b) you want to crop more at the top or at the bottom of the image: adapt
                           the last two thresholds, that correspond to additional cropping in mm at the top
                           and at the bottom of the image respectively. This is an interesting option in case
                           of an image that presents large white bands different from aponeuroses: by eliminating
                           these bands, you increase your chances to obtain a correct processing.
                4) The search for aponeuroses begins (1st: superficial aponeurosis; 2nd: deep aponeurosis)
                        A window asks you to validate the contour of the aponeurosis if it has been found.
                        If not satisfied, a second try is launched.
                        If the aponeurosis is not found or the contour does not satisfy the user, then
                        a linear approximation of the location of the aponeurosis is used.
                5) Fascicles are automatically looked for.
                6) The final image with the detected aponeuroses (in blue)
                   and the fascicles (in green) appears. Close it to move on to the following image.


        - Panoramic images
                1) The image appears with a window asking you to validate the start of the analysis
                2) Scale is automatically detected
                3) The image is cropped according to data labelled manually and stored in txt file.
                   The cropped image appears in a window. Close to move on.
                4) The image is divided vertically into sub-images to look for aponeuroses.
                   The size of sub-images depends on the width of the cropped image.
                   Per sub-images, the search for aponeuroses follows the same process as
                   for simple images. Superficial aponeurosis is searched all along the 
                   cropped image. Deep aponeurosis is searched only in the first half of the
                   cropped image, according to observations. Once all subimages have been processed,
                   aponeuroses are fitted with 2-degree polynomes (meaning, if aponeuroses are not found
                   in a sub-image, this is not a problem; to ensure the fitting of aponeuroses, there 
                   should be at least one portion of each aponeurosis found among all sub-images. However
                   the more portions of aponeuroses are detected, the better the estimation would be).
                5) The search for fascicles is realized automatically sub-image by sub-image before reunification.
                6) Two windows show 
                        - the original image with the intersection points between fascicles
                          and aponeuroses
                        - the cropped image with fascicles (in green) and aponeuroses (in blue)
                   Close them to move on to the following image.

        - When the data set has been fully processed, 10 plots are created to visualize results
                
                5 plots for simple images, 5 plots for panoramic images:
                - 1 plot for the comparison of calibration factors (auto vs. manual);
                - 1 plot for the comparison of fascicle lengths estimation (auto vs. manual)
                - 1 plot for the comparison of pennation angles estimation with superficial aponeurosis (auto vs. manual)
                - 1 plot for the comparison of pennation angles estimation with deep aponeurosis (auto vs. manual)
                - 1 plot for the comparison of muscle thickness estimation at discrete points (auto vs. manual)



Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
