# Mapping the relationship between tuberculosis burden and case notifications in Uganda

This repository contains code used to produce the results described in Chapter 3 of the doctoral thesis, "Assessing local health outcomes using spatially-resolved health surveillance data."


## Downloading data

All input data for this project was extracted from publicly-available reports on the 2014-15 National TB Prevalence Survey and the annual reports of the Uganda National TB and Leprosy Control Programme (NTLP), available from the [data portal](http://library.health.go.ug/) of the Uganda Ministry of Health. Extracted and formatted data can be shared upon request.


## Code organization

- `01_prepare_shapefile.R`: Load, clean, and save the district-level shapefile of Uganda which is used throughout later data preparation and modeling steps.
- `02_spatial_data_prep.R`: Prepare TB case notifications, prevalence survey data, and covariates for the joint spatial model.
- `03_run_tb_joint_model.R`: Run the notifications and prevalence survey joint spatial model. This script executes functions defined in `model_functions.R`. The core statistical model, written in Template Model Builder, is defined in `joint_completeness_tmb_model.cpp`
- `04a_visualize_model_data.R`: Visualize the TB prevalence survey data and case notifications prepared in `02_spatial_data_prep.R`.
- `04b_visualize_model_results.R`: Visualize model results, including figure generation scripts for all figures included in this chapter.
- `05_manuscript_number_plugging.R`: Organize results in order to reproduce findings and statistics included in the final chapter.


## Code execution

All input data, model results, and visualizations are stored in subdirectories of a central folder, defined as  `work_dir` at the top of each script:
- `raw_data/`: Stores input shapefile and CSV tables extracted from reports.
- `prepped_data/`: Prepared data used as inputs for the joint small area model.
- `model_results/`: Configuration and output objects from the joint small area model.
- `viz/`: Visualizations and output from other secondary analyses.

Once the `work_dir` directory is organized in this format, the scripts can be run in order to prepare and run the tuberculosis multi-source spatial model.