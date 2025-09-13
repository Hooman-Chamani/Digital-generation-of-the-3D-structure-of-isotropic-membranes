README – 3D Porous Structure Reconstruction from 2D SEM Images

This MATLAB project is designed to generate realistic 3D porous structures from a single 2D Scanning Electron Microscope (SEM) image. 

How It Works

The user provides a 2D SEM image (DN.tif) placed in the Raw/ directory.

Each reconstruction's central slice is saved side-by-side with the original binary image.

The optimized parameters for each reconstruction are also saved.

Getting Started

Place your SEM image in a folder called Raw/ and name the file DN.tif.

Launch MATLAB and run the Main function.

When prompted, enter the pixel length in microns of your input image.

Confirm by typing y to proceed with the reconstruction process.

By default, the code will generate 32 different 3D structures.

Output

The following folders will be used or created during execution:

Raw/: Should contain the input SEM image named DN.tif.

Rec/: Will contain output .tif images showing the reconstructed 2D slice next to the original.

X/: Will contain .mat files with the optimized parameters used for each reconstruction.

Configuration Options

Within the code, you can adjust the following:

Para: Set to 1 for parallel processing (parfor), or 0 for serial execution.

Repeats: Number of 3D structures to generate (default is 32).

Requirements

This code requires the following MATLAB toolboxes:

Image Processing Toolbox

Statistics and Machine Learning Toolbox (for Bayesian optimization)

Parallel Computing Toolbox (if running in parallel)



Contact

For questions, improvements, or collaboration inquiries, please contact:

[Hooman Chamani and Arash Rabbani]
Email: [Chamani.hooman@gmail.com 
A.Rabbani@leeds.ac.uk]
