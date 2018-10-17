# LFMotionECCV
Code for the ECCV 2018 paper:

>Sizhuo Ma, Brandon M. Smith, Mohit Gupta. _3D Scene Flow from 4D Light Field Gradients_. European Conference on Computer Vision (ECCV 2018)

Website: http://wisionlab.cs.wisc.edu/project/lfsceneflow

GitHub: https://github.com/sizhuom/LFMotionECCV

## How to Use
Run initLFMotion.m to set up the necessary MATLAB paths. Then take a look at runTest.m.

Dataset used in the paper can be downloaded here: https://uwmadison.box.com/s/abmn5baub87b4i5tqxz6r76cwyi7003t
Please extract the dataset under the folder ./data

## Dependencies
__Light Field Toolbox for MATLAB__: http://dgd.vision/Tools/LFToolbox/

__LightField_GeoCalibration_ver2__: https://sites.google.com/site/yunsubok/lf_geo_calib (if you want to use your own light field dataset)

__cocolib__: http://cocolib.net/ (if you want to generate depth maps by yourself)
