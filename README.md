# AA Tau ALMA Analysis

This repository contains the final data reduction and analysis pipeline for AA Tau ALMA observations. The imaging parameters were systematically optimized through extensive testing to achieve the best image quality and sensitivity.

## Overview

AA Tau is a young T Tauri star at a distance of 145 pc. This repository processes ALMA observations of multiple molecular lines from three observing campaigns (2013 SG1, 2013 SG2, and 2015 SG1).

## Source Parameters

- **Distance**: 145 pc
- **Stellar mass**: 0.62 M☉ (updated from previous value of 0.85 M☉)
- **Inclination**: 59.1°
- **Position angle**: 93.0° (measured East of North to receding side)
- **Systemic velocity**: 6.445 km/s
- **Disk center offset**:
  - RA offset (dx0): 0.0065"
  - Dec offset (dy0): -0.2573"

## Molecular Lines

### 2013 SG2 Dataset
- **N2H+** (spw 0): 279.512 GHz, 0.3 km/s channels
- **DCO+** (spw 2): 288.144 GHz, 0.3 km/s channels
- **H2CO** (spw 6): 290.623 GHz, 0.3 km/s channels

### 2013 SG1 Dataset
- **HCN** (spw 0): 265.886 GHz, 0.2 km/s channels (multiple hyperfine components)
- **HCO+** (spw 1): 267.558 GHz, 0.2 km/s channels

### 2015 SG1 Dataset
- **CN** (spw 1): 340.248 GHz, 0.2 km/s channels (multiple hyperfine components)
- **CN** (spw 3): 340.032 GHz, 0.2 km/s channels (multiple hyperfine components)
- **C18O** (spw 5): 329.331 GHz, 0.2 km/s channels
- **13CO** (spw 6): 330.588 GHz, 0.2 km/s channels

## Imaging Parameters

### Optimized Parameters
The following parameters were determined through systematic testing (see `/cleaning_tests/` in parent directory):

- **Weighting**: Briggs with robust = 0.5
- **UV taper**: ['0.05arcsec', '0.1483arcsec', '26deg']
  - Elliptical taper applied to circularize the synthesized beam
  - Major axis: 0.05"
  - Minor axis: 0.1483"
  - Position angle: 26°
- **Multiscale clean**: scales = [0, 5, 10, 20] pixels
- **Image size**: 500 × 500 pixels
- **Pixel size**: 0.03 arcsec/pixel
- **Spectral setup**: 100 channels starting at -3.0 km/s

### Keplerian Mask Parameters
- **Beam smoothing**: 2.0 × synthesized beam
- **Line width model**: dV(r) = 500 m/s × (r/1")^(-0.5)
- **Center offsets**: Included in mask generation for accurate spatial alignment

## Cleaning Procedure

The pipeline follows a standard two-step cleaning process:

### Step 1: Initial Dirty Clean
- Run `tclean` with `niter=0` to create a dirty cube
- High threshold (5 mJy) to prevent any cleaning
- Purpose: Generate cube for mask creation

### Step 2: Keplerian Mask Creation
- Use `keplerian_mask.py` (by Richard Teague)
- Create velocity mask based on disk geometry
- Mask includes all expected line emission based on Keplerian rotation
- Accounts for:
  - Disk inclination and position angle
  - Stellar mass and distance
  - Radially varying line width
  - Center position offsets

### Step 3: RMS Calculation
- Calculate RMS from line-free channels (channels 0-30)
- Set cleaning threshold to **2 × RMS**

### Step 4: Final Clean
- Run `tclean` with Keplerian mask
- Use calculated 2×RMS threshold
- Maximum 50,000 iterations
- Produces final cleaned data cube

## Usage

### Requirements
- CASA (Common Astronomy Software Applications)
- Access to continuum-subtracted measurement sets

### Running the Pipeline

From within CASA:

```python
execfile('clean_all_lines.py')
```

The script will:
1. Process all 10 spectral lines sequentially
2. Create CASA image files in `casa_images/`
3. Export FITS files to `fits_products/`
4. Print detailed progress and diagnostics

### Output Files

For each spectral line, the following files are created:

**CASA Images** (in `casa_images/`):
- `AATau_{molecule}_contsub_clean0.*` - Initial dirty cube and mask
- `AATau_{molecule}_contsub_clean1.*` - Final cleaned cube

**FITS Files** (in `fits_products/`):
- `AATau_{molecule}_contsub_clean1.fits` - Final cleaned cube
- `AATau_{molecule}_contsub_clean0.mask.fits` - Keplerian mask

## Files in This Repository

- **`clean_all_lines.py`**: Main cleaning script for all molecular lines
- **`keplerian_mask.py`**: Keplerian mask generation script (by Richard Teague)
- **`README.md`**: This file

## Testing History

The parameters used in this pipeline were determined through systematic testing:

1. **Weighting tests**: Compared natural, Briggs (robust 0.5), and uniform weighting
2. **UV taper tests**: Tested various taper sizes to circularize the beam
3. **Threshold tests**: Compared 2×, 3×, and 5× RMS thresholds
4. **Stellar mass update**: Updated from 0.85 to 0.62 M☉ based on latest measurements

All tests showed that:
- Robust 0.5 provides the best balance of sensitivity and resolution
- Elliptical UV taper successfully circularizes the beam without excessive resolution loss
- 2× RMS threshold provides good sensitivity while avoiding over-cleaning

## Notes

- The 2015 dataset uses a different directory naming convention (`contsub_2015` instead of `contsub_2015SG1`)
- HCN and CN lines include multiple hyperfine components (all specified in `restfreqs`)
- All frequency values in the script are in GHz and are converted to Hz for `make_mask()`
- The pipeline is robust to failures - if one line fails, processing continues with the next line

## Citation

If you use this pipeline or the Keplerian mask script, please cite:
- Keplerian mask script: Teague, R. (2020)

## Contact

For questions about this pipeline, please contact Jea Adams.

## Version History

- **v1.0** (November 2024): Initial final version with optimized parameters
  - Updated stellar mass to 0.62 M☉
  - Implemented elliptical UV taper for beam circularization
  - Added center position offsets to mask generation
  - Set cleaning threshold to 2× RMS
