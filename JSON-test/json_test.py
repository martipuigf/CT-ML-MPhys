#!/usr/bin/env python3

# Import packages
import os
import numpy as np # Who does not use Numpy?

has_mpl = True
try:
    import matplotlib # To plot images
    import matplotlib.pyplot as plt # Plotting
    from matplotlib.colors import LogNorm # Look up table
    from matplotlib.colors import PowerNorm # Look up table

    font = {'family' : 'serif',
             'size'   : 15
           }
    matplotlib.rc('font', **font)

    # Uncomment the line below to use LaTeX fonts
    # matplotlib.rc('text', usetex=True)
except:
    has_mpl = False

from tifffile import imwrite # Write TIFF files

from gvxrPython3 import gvxr # Simulate X-ray images
import json2gvxr
from tifffile import imwrite
import tomopy


number_of_projections = 50
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
plt.plot(detector_response[:,0] * 1000, detector_response[:,1] * 1000)
plt.xlabel('Incident energy: E (in keV)')
plt.ylabel('Detector energy response: $\\delta$(E) (in keV)')

plt.tight_layout()

x_ray_projections = np.array(gvxr.computeProjectionSet(0,0,0, "cm", number_of_projections, 
                                                 angular_step)).astype(np.single)
imwrite('proj.tif', x_ray_projections) #Saves the projections as a stack o TIFF files


dark = np.zeros(np.shape(x_ray_projections))
energy_bins = gvxr.getEnergyBins('keV')
photon_per_bin = gvxr.getPhotonCountEnergyBins()

total_energy = 0.0
for energy, count in zip(energy_bins, photon_per_bin):
    total_energy += energy * count
flat = np.ones(x_ray_projections.shape) * total_energy

projections = tomopy.normalize(x_ray_projections, flat, dark)
projections = tomopy.minus_log(projections)
    
theta_rad = np.array(theta_deg) * np.pi / 180
rot_centre = x_ray_projections.shape[2] / 2
#rot_centre = tomopy.find_center(projections, theta_rad, init=rot_centre, ind=0, tol=0.01)
print("Projection sizes:", x_ray_projections.shape)
print("Rotation centre:", rot_centre)



# recon = tomopy.recon(projections,
#                           theta_rad,
#                           center=rot_centre,
#                           algorithm='fbp')

# plt.imshow(recon[int(projections.shape[1]/2), :, :])

    

for n in range(int(number_of_projections)):
    plt.figure(figsize=(10, 5))
    plt.title("Image simulated using gVirtualXray\nusing a logarithmic colour scale")
    plt.imshow(x_ray_projections[n]/1e6, cmap="gray", norm=PowerNorm(gamma=1./2.))
    plt.colorbar(orientation='vertical');
#    plt.savefig('square'+str(n)+'.png')
#    np.savetxt(fname=('array_'+str(n)+'.txt'), X=x_ray_image[n])