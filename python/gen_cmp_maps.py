import healpy as hp
import numpy as np
from tqdm import tqdm  # Import tqdm for progress tracking

# Define folder to store CMB maps
folder1 = f"/scratch/shambhavij.sps.iitmandi/fortran_code/genmaps/"  

# Set Healpy parameters
nside = 1024  # Healpix resolution parameter
nlmax = 2507  # Maximum multipole moment for power spectra
num_map = 1000  # Number of CMB maps to generate

# Number of pixels in the map based on NSIDE
npix = hp.nside2npix(nside)

# Set a seed value for reproducibility of random numbers
SEEDVALUE = 1234
np.random.seed(SEEDVALUE)

# Create an empty array to store all maps (each map in a column)
all_maps = np.zeros((npix, num_map))

# Function to generate CMB maps from a given power spectrum (cl file)
def generate_maps(cl_file):
    cls = np.loadtxt(cl_file)  # Load the power spectrum from the file
    Map = hp.synfast(cls, nside, lmax=nlmax)  # Generate a synthetic map using Healpy
    return Map

# Generate multiple CMB maps with a progress bar and store them in the array
print("Generating CMB maps and storing them in an array:")
for i in tqdm(range(num_map), desc="CMB Map Generation"):
    all_maps[:, i] = generate_maps('cl.dat')  # Store each map in a column

# Save the entire array as a single .npy file
np.save(folder1 + "all_cmb_maps.npy", all_maps)

print("All CMB maps generated and saved successfully!")

