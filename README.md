# IR.Cai

(c) Jonathan Goodman, Benji Rowlands  
2021-2024  
jmg11@cam.ac.uk  
University of Cambridge  

## Overview

This repository contains an implementation of the IR.Cai algorithm described in
the paper "Towards automatically verifying chemical structures: the powerful
combination of <sup>1</sup>H NMR and IR spectroscopy". 

The script is contained in `cai.py`. An example input file, with parameters set
to the values used to obtain the results in the paper, is contained in
`input.txt`.

## Requirements
- Operating system: this software has been tested on macOS Sonoma 14.5. It
  should work on other operating systems, but has not been tested on them.
- Python version: this software was developed using Python 3.12.

## Installation

To run this software, ensure you have Python 3.12 installed along with the
following package:
- `matplotlib` (version 3.9.0)

To configure the environment correctly, follow the following steps:
1. **Install Python**: Ensure you have Python 3.12 installed. You can
   download it from the [official Python
   website](https://www.python.org/downloads/release/python-3124/).

2. **Create a virtual environment**:
    ```
    python3.12 -m venv .venv
    source .venv/bin/activate 
    ```

3. **Activate the virtual environment and install the required package**:
    ```
    source .venv/bin/activate
    pip install matplotlib==3.9.0
    ```

## Usage

1. **Prepare your input file**: The input file should contain sections for
   settings, experimental data files and calculation data files. Each section
   should be denoted by tags `<Settings>`, `<Experiments>` and
   `<Calculations>`, and each section should be closed with a tag of the form
   `</Settings>`.

2. **Run the script**: Execute the script from the command line with your input file as an argument.
```bash
python /path/to/cai.py your_input_file.txt
```

### Example Input File
```
<Settings>
print_csv spectra scaling_factor summary
minimum_wavenumber 1250
maximum_wavenumber 1600
defined_scaling_factor 0.98
min_scale_factor 0.95
temperature 300
max_scale_factor 1.0
optimise_scaling_factor False
boltzmann_cutoff 20.0
broadening 12
save_path /path/to/your/save/file.csv (optional)
</Settings>

<Experiments>
/path/to/your/experimental/data/file.txt
</Experiments>

<Calculations>
/path/to/your/dir/containing/calculation/files
</Calculations>

```

## Output

The software generates the following outputs:
- A `.log` file containing detailed logs of the process.
- A `.csv` file with the experimental and calculated spectra data.
- A graph in `.pdf` format showing the experimental and calculated spectra.
- A summary output appended to a specified path.

## Example data
Usage of the software can be demonstrated using the molecules from the test set
in the paper. All of the necessary files are contained in the Apollo repository:
`https://doi.org/10.17863/CAM.110235`. This folder contains its own README detailing how to reproduce
the IR.Cai scores from the paper in conjunction with the IR.Cai script.

## Contact

For any questions or issues, please contact Jonathan Goodman at jmg11@cam.ac.uk.