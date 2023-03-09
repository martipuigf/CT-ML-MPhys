#!/usr/bin/env python3

# Import packages
import os
import numpy as np # Who does not use Numpy?


from tifffile import imwrite # Write TIFF files

from gvxrPython3 import gvxr # Simulate X-ray images
import json2gvxr
from tifffile import imwrite
from tifffile import imread


number_of_projections = 3000
angular_step = 360 / number_of_projections
theta_deg = np.linspace(0.0, angular_step * number_of_projections, num=number_of_projections, endpoint=False)

print("Number of projections:", theta_deg.shape[0])
print("Angle between successive projections:", angular_step)
print("First angle:", theta_deg[0])
print("Last angle:", theta_deg[-1])

#fname = 'thick_test.stl'

json2gvxr.initGVXR('notebook.json', "OPENGL")

json2gvxr.initSamples(verbose=2)

json2gvxr.initSourceGeometry()

spectrum, unit, k, f= json2gvxr.initSpectrum(verbose=0)
energy_set = sorted(spectrum.keys())
count_set = []
for energy in energy_set:
    count_set.append(spectrum[energy])
    
json2gvxr.initDetector("notebook.json")
gvxr.moveToCentre("cube")
detector_response = np.loadtxt("energyResponseDetector.txt")


x_ray_projections = np.array(gvxr.computeProjectionSet(0,0,0, "cm", number_of_projections, 
                                                 angular_step)).astype(np.single)
imwrite('3000_proj.tif', x_ray_projections) #Saves the projections as a stack o TIFF files
