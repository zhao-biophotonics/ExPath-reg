# ExPath-reg
The code is designed for rigid image registration between two images, specifically for pre-expansion image and post-expansion image. It can handle z stacks and generate 2D maximum intensity projection for image registration. This code doesnot handle 3D image registration.
This code has been tested on MatLab 2014b.
# Dependency
All the dependencies have been included.
# Installation guide
Download the code in a folder and put the folder path in the search path of Matlab. Unzip the helper.zip under the same folder and make sure the folder is in the search path of Matlab. Installation time is less than 1 minute. The code doesnot require any non-standard hardware to run.
# Instruction for use
Please see the annotation of the code. It is straightforward, first change the parameters section of the code, execute the code in matlab, select the 'map' (pre-expansion) and 'query' (post-expansion) images, and the code will automatically output the results in a folder.
Example parameters are included in the annotation of the code.
# Expected output
A folder contains the registered images (multicolor channels are splited into individual single-channel images), all keypoints for registrations, matlab variables and a note with information such as expansion factor. 
# Example data can be downloaded via https://drive.google.com/open?id=1lG8wp_UIl_ki0XPKiCKwT6DCZP_-3wWQ
