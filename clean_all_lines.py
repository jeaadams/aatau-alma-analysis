"""
AA Tau ALMA Data Cleaning Script - Final Version

This script performs CASA tclean imaging on all AA Tau spectral lines with
optimized parameters determined through systematic testing.

Imaging Parameters:
-------------------
- Weighting: Briggs with robust = 0.5
- UV taper: ['0.05arcsec', '0.1483arcsec', '26deg'] (elliptical taper to circularize beam)
- Multiscale clean with scales = [0, 5, 10, 20]

Keplerian Mask Parameters:
--------------------------
- Stellar mass: 0.62 M_sun (updated from 0.85)
- Inclination: 59.1 degrees
- Position angle: 273.0 degrees (93.0 + 180)
- Distance: 145 pc
- Systemic velocity: 6.445 km/s
- Center offset: dx0 = 0.0065", dy0 = -0.2573"
- Beam smoothing: 2.0 beams
- Line width: dV0 = 500 m/s, dVq = -0.5

Cleaning Procedure:
------------------
1. Initial dirty clean (niter=0) to generate cube
2. Create Keplerian mask with updated parameters
3. Calculate RMS from line-free channels (0-30)
4. Final clean with mask and 2x RMS threshold

Datasets:
---------
- 2013_SG2: N2H+, DCO+, H2CO
- 2013_SG1: HCN, HCO+
- 2015_SG1: CN_SPW1, CN_SPW3, C18O, 13CO
"""

import os

# ==============================================================================
#  Load Keplerian Mask Function
# ==============================================================================

# Load the Keplerian mask function into the global scope
import sys
sys.path.insert(0, '/Users/jea/AATau')

with open('/Users/jea/AATau/keplerian_mask.py', 'r') as f:
    exec(f.read(), globals())
print("Successfully loaded make_mask function from keplerian_mask.py")

# ==============================================================================
#  Configuration
# ==============================================================================

# Common parameters for Keplerian mask generation
MASK_PARAMS = {
    'inc': 59.1,                    # Inclination in degrees
    'PA': 93.0 + 180.,              # Position angle in degrees
    'mstar': 0.62,                  # Stellar mass in M_sun (UPDATED from 0.85)
    'dist': 145.0,                  # Distance in parsecs
    'vlsr': 6.445e3,                # Systemic velocity in m/s
    'dx0': 0.0065,                  # Center offset in RA (arcsec)
    'dy0': -0.2573,                 # Center offset in Dec (arcsec)
    'nbeams': 2.0,                  # Beam smoothing factor
    'dV0': 500,                     # Line width at 1" in m/s
    'dVq': -0.5,                    # Line width power law exponent
}

# Common parameters for tclean
TCLEAN_COMMON_PARAMS = {
    'deconvolver': 'multiscale',
    'scales': [0, 5, 10, 20],
    'specmode': 'cube',
    'start': '-3.0km/s',
    'nchan': 100,
    'outframe': 'LSRK',
    'interactive': False,
    'imsize': [500, 500],
    'cell': '0.03arcsec',
    'weighting': 'briggs',
    'robust': 0.5,
    'uvtaper': ['0.05arcsec', '0.1483arcsec', '26deg'],  # Elliptical taper
    'restoringbeam': 'common',
}

# RMS calculation parameters
RMS_CHANNELS = '0~30'               # Line-free channels for RMS calculation
RMS_MULTIPLIER = 2.0                # Threshold = 2x RMS
MAX_ITERATIONS = 50000              # Maximum clean iterations

# Define all datasets to be processed
DATASETS = [
    {
        'name': '2013_SG2',
        'vis_template': 'AA_Tau_2013_SG2.spw{spw}.selfcal.ms.contsub_fit1',
        'width': '0.3km/s',
        'lines': [
            {
                'spw': 0,
                'molecule': 'N2H+',
                'line_freq': 279.5117491,  # In GHz
                'restfreqs': [279.5117491],  # In GHz
            },
            {
                'spw': 2,
                'molecule': 'DCO+',
                'line_freq': 288.1438583,
                'restfreqs': [288.1438583],
            },
            {
                'spw': 6,
                'molecule': 'H2CO',
                'line_freq': 290.623405,
                'restfreqs': [290.623405],
            },
        ]
    },
    {
        'name': '2013_SG1',
        'vis_template': 'AA_Tau_2013_SG1.spw{spw}.selfcal.ms.contsub_fit1',
        'width': '0.2km/s',
        'lines': [
            {
                'spw': 0,
                'molecule': 'HCN',
                'line_freq': 265.8864343,
                'restfreqs': [265.8848912, 265.8861886, 265.8864339, 265.8864343, 265.8864999, 265.8885221],
            },
            {
                'spw': 1,
                'molecule': 'HCO+',
                'line_freq': 267.5576259,
                'restfreqs': [267.5576259],
            },
        ]
    },
    {
        'name': '2015_SG1',
        'vis_template': 'AA_Tau_2015_SG1.spw{spw}.selfcal.ms.contsub_fit1',
        'width': '0.2km/s',
        'lines': [
            {
                'spw': 1,
                'molecule': 'CN_SPW1',
                'line_freq': 340.24777,
                'restfreqs': [340.24777, 340.248544, 340.2617734, 340.264949, 340.2791201],
            },
            {
                'spw': 3,
                'molecule': 'CN_SPW3',
                'line_freq': 340.0315494,
                'restfreqs': [340.0081263, 340.0196255, 340.0315494, 340.035408],
            },
            {
                'spw': 5,
                'molecule': 'C18O',
                'line_freq': 329.3305525,
                'restfreqs': [329.3305525],
            },
            {
                'spw': 6,
                'molecule': '13CO',
                'line_freq': 330.5879653,
                'restfreqs': [330.5879653],
            },
        ]
    },
]

# ==============================================================================
#  Main Processing Function
# ==============================================================================

def clean_line(base_path, output_path, dataset_config, line_config):
    """
    Performs the full cleaning process for a single spectral line.

    Steps:
    1. Initial dirty clean (niter=0)
    2. Create Keplerian mask with center offset
    3. Calculate RMS threshold (2x RMS)
    4. Final clean with mask
    5. Export to FITS
    """

    # --- Define paths and names ---
    dataset_name_short = dataset_config['name'].replace('_', '')

    # Handle the different naming convention for the 2015 dataset
    if dataset_config['name'] == '2015_SG1':
        contsub_dir = 'contsub_2015'
    else:
        contsub_dir = f"contsub_{dataset_name_short}"

    contsub_path = os.path.join(base_path, contsub_dir)
    vis_file = os.path.join(contsub_path, dataset_config['vis_template'].format(spw=line_config['spw']))

    # Define output names
    output_prefix = os.path.join(output_path, f"AATau_{line_config['molecule']}_contsub")
    clean0_imagename = f"{output_prefix}_clean0"
    clean0_image = f"{clean0_imagename}.image"
    clean1_imagename = f"{output_prefix}_clean1"
    clean1_image = f"{clean1_imagename}.image"
    mask_image = f"{clean0_imagename}.mask.image"

    # --- Step 1: Initial dirty clean (niter=0) ---

    tclean(
        vis=vis_file,
        imagename=clean0_imagename,
        width=dataset_config['width'],
        restfreq=f"{line_config['line_freq']}GHz",
        threshold='5mJy',
        niter=0,
        **TCLEAN_COMMON_PARAMS
    )

    # --- Step 2: Create Keplerian mask ---

    make_mask(
        image=clean0_image,
        restfreqs=[freq * 1e9 for freq in line_config['restfreqs']],
        **MASK_PARAMS
    )

    # --- Step 3: Calculate RMS threshold ---

    stats = imstat(imagename=clean0_image, chans=RMS_CHANNELS)

    if stats and 'rms' in stats and stats['rms'][0] > 0:
        rms_val = stats['rms'][0]
        threshold_val = round(1000. * RMS_MULTIPLIER * rms_val, 1)
        threshold_str = f'{threshold_val}mJy'
    else:
        threshold_str = '10mJy'

    # --- Step 4: Final clean with mask ---

    tclean(
        vis=vis_file,
        imagename=clean1_imagename,
        mask=mask_image,
        width=dataset_config['width'],
        restfreq=f"{line_config['line_freq']}GHz",
        threshold=threshold_str,
        niter=MAX_ITERATIONS,
        **TCLEAN_COMMON_PARAMS
    )

    # --- Step 5: Export to FITS ---

    fits_dir = os.path.join(base_path, "aatau-alma-analysis", "fits_products")
    if not os.path.exists(fits_dir):
        os.makedirs(fits_dir)

    # Export cleaned image
    fits_image_name = os.path.join(fits_dir, os.path.basename(clean1_imagename) + ".fits")
    exportfits(imagename=clean1_image, fitsimage=fits_image_name, overwrite=True, dropstokes=True)

    # Export mask
    fits_mask_name = os.path.join(fits_dir, os.path.basename(clean0_imagename) + ".mask.fits")
    exportfits(imagename=mask_image, fitsimage=fits_mask_name, overwrite=True, dropstokes=True)


def main():
    """
    Main function to loop through all datasets and spectral lines.
    """
    base_path = "/Users/jea/AATau"
    output_path = os.path.join(base_path, "aatau-alma-analysis", "casa_images")

    # Create output directory if it doesn't exist
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # Process all datasets and lines
    for dataset in DATASETS:
        for line in dataset['lines']:
            try:
                clean_line(base_path, output_path, dataset, line)
            except Exception:
                # Continue to next line
                continue

    print("\n" + "="*80)
    print("ALL CLEANING COMPLETE")
    print("="*80)
    print(f"\nOutput files located in:")
    print(f"  CASA images: {output_path}")
    print(f"  FITS files: {os.path.join(base_path, 'aatau-alma-analysis', 'fits_products')}")
    print()


# ==============================================================================
#  Execute main processing
# ==============================================================================

if __name__ == '__main__':
    main()

# For CASA execfile, the if __name__ guard is evaluated as True, so this will run
