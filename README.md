# Flipper: a tool for automatically detecting conformational changes from Ringer measurements

Author: Tim Stachowski, PhD and Marcus Fischer, PhD | email: tim.stachowski@stjude.org
#### [paper using Flipper](https://)
If you use Flipper in your work, please cite:

### About Flipper
Flipper is a tool used for finding side-chain conformational changes between two sets of Ringer electron density (sigma) measurements. Flipper uses a peak finding algorithm to find heighted areas of electron density that correspond to side chain conformers. First, based on this, Flipper reports on whether a residue gains or loses a conformation between two samples being compared. Second, Flipper detects if the major conformer changes (FLIP!) between the same residue in two samples through peak integration.

The core functions are stored in the `ringer_tools.py` module. The actual script to run Flipper is `flipper.py`. 

## Installation
Flipper can be installed by simply cloning the git repository. 

## Requirements
Flipper requires Python >= 3.7, NumPy, SciPy, Pandas, and matplotlib are installed. These packages are often installed by default but can installed using package managers such as [Anaconda](https://continuum.io/downloads). The current program was tested and built on Mac OS.

## Input files
Flipper uses the standard output files from Ringer (files ending in `_ringer.csv`). Any modifications to the column positioning or addition of whitespace will most likely break Flipper. Flipper works best when the Ringer output files are broken down into separate files for each chain. For example, if your protein contains chains A and B, you should feed Flipper separate files:   `prot_chainA_ringer.csv` and `prot_chainB_ringer.csv`.

## Usage
Flipper is simply run as:
`python flipper.py -f1 file1_ringer.csv -f2 file2_ringer.csv`

The `examples` folder contains two example files from the manuscript.


Options you may want to see are:
```
  -h, --help            show this help message and exit
  -i1 [FILENAME1], --filename1 [FILENAME1]
                        input file #1 - must be a standard Ringer output CSV
                        with only one chain
  -i2 [FILENAME2], --filename2 [FILENAME2]
                        input file #2 - must be a standard Ringer output CSV
                        with only one chain
  -t SIGMATHRESHOLD, --sigmathreshold SIGMATHRESHOLD
                        sigma threshold for peak finding. default = 0.3
  -plot PLOT, --plot PLOT
                        Save individual plots showing peak finding results?
                        This is slow. default = False
```
The `-t` flag controls the sigma threshold used for peakfinding. The default is 0.3. 
The `-plot` flag controls printing plots. These plots show the data range used for integration (based on the sigma threshold) and          shows detected peaks by marking them with an `x`.


## Results
As the program runs, the current status will be printed to the screen.
Flips are printed like:
```
Identifying FLIPS....

FLIP! at residue 85 chi1
FLIP! at residue 127 chi1
FLIP! at residue 201 chi1
FLIP! at residue 85 chi2
```
Output files include:
```
peak_finder_file1_ringer.csv - includes: (1) number of peaks, (2) peak angles, and (3) peak area for each residue and dihedral for each residue in file1
peak_finder_file2_ringer.csv - dito for file2
file1_ringer_file2_ringer_gain_loss_peaks.csv - a list of residues and the change in number of peaks between the two files (file1-file2)
 
 
 
 
 
 
 
 
 

