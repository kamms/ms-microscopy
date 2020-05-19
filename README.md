# ms-microscopy
MS-microscopy uses a comprehensive subcellular location specific reference data to predict the localization of a given bait protein from interactome experiment (e.g. BioID) data.

## System requirements
MS-microscopy is not operating-system specific, but has only been tested on windows 10 and requires python 3 and R 3.5.

python 3.7:
- with matplotlib and numpy
R 3.5:
- with shiny and ggplot2

## Included files
msmicroscopy -directory contains all files needed to run the msmic.py script within, as well as example output. This directory The output directory also contains sample output made with the included input files.

msmic_online -directory contains files needed to set up ms-microscopy in R shiny. reference_data.csv is generated with the python script by default each time it is run. 

## Input files and parameters for msmic.py 
1. training_set_identifiers: This file contains information on what bait corresponds to what cellular compartment. Important columns are Bait, Bait ID, and Organelle. The rest are unused. 
2. training_set_datafile: This file contains spectral count values for the baits identified as localization markers in file 1. The format is three columns: Bait (this should be identical to Bait column in identifier file), Prey (Uniprot IDs preferred. Can be other, but these should be from the same database as in query file), AvgSpec (spectral count values).
3. QUERY_file: Same format as training set datafile, but should contain the data for query experiments. The values will be normalized, but filtering will need to be done beforehand. 
4. abbreviations: Three columns: Organelle, abbreviation (for output), and color. The color column is currently unused.

Fifth parameter should be the output directory. If it doesn't exist, it will be created. 

## Usage: 
msmic.py from command line: 
- python3 msmic.py training_set_identifiers training_set_datafile QUERY_file abbreviations outputDirectory
msmic_online:
- Can be run and tested via RStudio (version 1.2.5 tested)
- Can be uploaded to e.g. shinyapps.io 

## Output
In the specified output directory: 
- results.tsv contains scores for each bait for each subcellular localization.
- Figures.pdf contains the drawn polar plots for all baits, one plot per page.
- Loc_info is a text file containing the size of each reference dataset (as a sum of psm values)
- Localization_preys is a text file with four columns: Bait, localization, prey, and avgspec. It is a full description of the QUERY dataset. It shows which preys drive what localization for each bait.
- unique_set -files contain unique preysets for each reference localization.
- In for_R_server -directory, reference_data.csv can be used as a reference database for R server.

## License

This project is covered under the [Apache 2.0 License](https://github.com/kamms/ms-microscopy/blob/master/LICENSE)
