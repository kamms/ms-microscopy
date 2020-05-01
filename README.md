# ms-microscopy
msmicroscopy -directory contains all files needed to run the msmic.py script within, as well as example output. 
Input files: 
1. training_set_identifiers: This file contains information on what bait corresponds to what cellular compartment. Important columns are Bait, Bait ID, and Organelle. The rest are unused. 
2. training_set_datafile: This file contains spectral count values for the baits identified as localization markers in file 1. The format is three columns: Bait (this should be identical to Bait column in identifier file), Prey (Uniprot IDs preferred. Can be other, but these should be from the same database as in query file), AvgSpec (spectral count values).
3. QUERY_file: Same format as training set datafile, but should contain the data for query experiments. The values will be normalized, but filtering will need to be done beforehand. 
4. abbreviations: Three columns: ORganelle, abbreviation (for output), and color. The color column is currently unused.

This directory also contains generateLocalizationReferenceForMSMIC.ipynb -notebook, which can be used to convert data from output directory into data usable as reference set for the shiny server.

msmic_online -directory contains files needed to set up ms-microscopy in R shiny. reference_data.csv can be generated with the notebook in msmicroscopy -directory. 
