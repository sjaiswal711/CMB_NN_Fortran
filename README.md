# Beam Convolution Simulations

This repository contains Fortran and Python codes for simulating satellite scanning of HEALPix sky maps with different beam models:

- **Elliptical Gaussian Beam** (idealized, parametric)
- **Realistic Beam** (based on pre-computed response grid)

The outputs can be used for forward modeling of time-ordered data (TOD), response matrices, and validation of scanning strategies.

---


---

## ⚙️ Code Descriptions

### 🔹 Elliptical Beam Convolution (`elliptical_convolution.f90`)
- **Input:**
  - `map.fits` → HEALPix sky map to be scanned.
- **Process:**
  - Simulates scanning with an **elliptical Gaussian beam** (`fwhm_x`, `fwhm_y`).
  - At each time step, computes:
    - Pointing direction (`pix_ring`).
    - Convolved temperature (detector signal).
- **Output:**
  - `convolved_map.dat` containing:
    ```
    time_step   pixel   convolved_temperature
    ```

---

### 🔹 Real Beam Convolution (`RealBeam_convolution.f90`)
- **Input:**
  - `grid.txt` → Pre-computed real beam response grid.
  - (Optionally) `map.fits` for testing.
- **Process:**
  - Uses MPI to parallelize over sky pixels.
  - Computes a **response matrix**:
    - Each row corresponds to a satellite pointing.
    - Each entry gives beam weights for neighboring HEALPix pixels.
  - Stores intermediate results per rank (`results_0.dat … results_47.dat`).
- **Output:**
  - `results_rank.dat` files (one per MPI rank), containing:
    ```
    node_id   pixel   count   weight1   weight2   ...
    ```

---

### 🔹 Response Matrix Neighbors (`response_matrix_neighbors.f90`)
- **Input:**
  - `results_rank.dat` files from Real Beam simulation.
- **Process:**
  - Extracts the **48 neighboring HEALPix pixel indices** for each response matrix pixel.
  - Facilitates linking response weights back to sky pixels.
- **Output:**
  - `neighbors_rank.dat` files with neighbor pixel indices.

---

### 🔹 Yearly Scan Frequency
- Generated as part of both simulations.
- Stores **hit counts**: how many times each HEALPix pixel was visited during the one-year scan.
- Useful for coverage maps and validation.

---

### 🔹 Python Scripts
- **`gen_cmb_maps.py`**: Generate CMB realizations from CAMB power spectra (`Cl`).
- **`convolve_maps.py`**: Convolve CMB + foreground maps with the beam model.
- Additional scripts for pre/post-processing.

---

## 🚀 Typical Workflow

1. **Generate or provide input maps:**
   - Place sky maps in `map.fits`.
   - Place beam grid in `grid.txt` (for real beam mode).

2. **Run elliptical beam scan:**
   ```bash
   gfortran elliptical_convolution.f90 -o elliptical_convolution
   ./elliptical_convolution
