# tims_tof_data_extraction

## Step 1

Install packages or create enviroment to facilitate packages, for example conda environment.

Clone the repository

  git clone https://github.com/Micholms/tims_tof_data_extraction

**To enable REST API for MetaboScape:**

Install the swagger_client using README information in python_client_generated the .zip file. 
Make sure the computer you are on can talk to the computer where MetaboScape is installed. 

**Install packages for processing of data:**
pandas, numpy, matplotlib.pyplot and seaborn.

**Optional: Install fiora:**

From here: https://github.com/BAMeScience/fiora.git

## Step 2
Run REST API for MetaboScape to analyze the experiments. "Rest_API_metaboscape.ipynb"
More information located in the file.

The output can then be further analyzed by "data_handling_REST_API.ipynb". 
For example, counting matches, pre-filtering and preparing the data for modelling.

Alternatively, the matching and filtering can be done through

  ! python target_match.py -i "test_unprocessed.csv" -o "test_filtered.csv" -t "plate_targets.csv"

  !python filter_data.py -i "test_filtered.csv" -o "test_preprocessed.csv"

Editing the filter_data script, parameters can be changed. Note: to filter data, Fiora in required! 

## Step 3
Optional: If fiora is installed, training using the pre-processed data could be done.


