# tims_tof_data_extraction

___Step 1___

Install packages or create enviroment to facilitate packages, for example conda environment.

**To enable REST API for MetaboScape:**

Install the swagger_client using README information in python_client_generated the .zip file. 
Make sure the computer you are one can talk to the computer where MetaboScape is installed. 

**Install packages for processing of data:**
pandas, numpy, matplotlib.pyplot and seaborn.

**Optional: Install fiora:**
From here: 


___Step 2___
Run REST API for MetaboScape to analyze the experiments. "Rest_API_metaboscape.ipynb"
More information located in the file.

The output can then be further analyzed by "data_handling_REST_API.ipynb". 
For example, counting matches, pre-filtering and preparing the data for modelling.

__Step 3__
Optional: If fiora is installed, training using the pre-processed data could be done.


