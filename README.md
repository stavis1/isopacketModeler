# IsopacketModeler
## A tool for identifying and modeling stable isotope label incorporation in peptides.

IsopacketModeler is a command line tool written in Python for identifying peptides that have incorporated stable isotopes and modeling isotope uptake levels.

## Installation

1. Make sure you have a working conda installation. If conda is not installed please follow the instructions at this page: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html
2. Download the repository as a .zip file from the green "<> Code" drop-down menu in the upper right of the GitHub page and extract the contents. Alternatively, if you have git installed run `git clone https://github.com/stavis1/isopacketModeler/`.
3. On the command line navigate to the root of the downloaded directory and run `conda env create -n isopacketModeler -f env/run.yml`. This will create a conda environment and install the necessary dependencies for this tool.
4. Run `conda activate isopacketModeler` to activate this environment.

## Usage

IsopacketModeler requires four text files to run:
1. A [TOML](https://toml.io/) formatted options file. 
2. A tab-separated file specifying the elemental composition of each amino acid.
3. A tab-separated file describing the isotopic label added to each raw file.
4. A tab-separated file describing the peptide-spectrum matches to use when searching for isotopically enriched peptides.

Additionally, the raw data must be provided as centroid mode mzML files.

Default versions of files 1 and 2, and templates for files 3 and 4 can be found in the `templates` folder.

File 3, the experimental design file, must have at least two columns 'file' and 'label' which list the mzML files and their heavy label isotopes, respectively. The allowable isotopes are 'H[2]', 'C[13]', 'N[15]', 'O[17]', 'O[18]', 'S[33]', 'S[34]', and 'S[36]'. Control samples should leave this field blank. Only one label is allowed per file. Other columns can be included. The information in these columns will be kept associated with the peptide data and will be included in the output but will otherwise be ignored. This can be used to, for example, keep track of time points or conditions in an experiment. 

File 4, the PSM file, must include columns for peptide sequence, mzml file, scan #, charge, and protein names. These columns do not need to have specific names, rather, the names of the columns corresponding to these pieces of information will be specified in file 1, the options file. This means that tab-separated PSM files from multiple tools can be used with minimal processing. Scan # can be either the MS2 or the parent MS1 scan for the PSM. The charge column must have only integer values, e.g. 2 not +2 and cannot be blank or unknown. The mzML files must be either in the same directory that the command is run from or must be specified as paths relative to that directory. Flanking residues in peptide sequences will be ignored if they are separated by '.', e.g. 'P.EPTID.E' will be interpreted as 'EPTID'. Post-translational modifications can be specified as any ASCII character that is neither a whitespace nor an open bracket or they can be specified as bracket bounded strings of characters following a character. For example either 'm' or 'M[15.9949]' could represent a modified residue. Either of these options must be added to File 2, the amino acid formula file. PTMs that use the bracket notation must include the preceding character in File 2, so an example entry line would look like `M[15.9949]	5	9	1	2	1	0`.

To run the tool ensure that the `isopacketModeler` conda environment is active then run `python -m isopacketModeler --options options.toml`. 
