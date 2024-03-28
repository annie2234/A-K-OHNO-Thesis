This repository contains the orignal DEFRA measurement data for London Marylebone for the years 2005, 2010 and 2020.

**The LM data formatting notebook**
Processes the orignal data and calculates the OH concentrations for each temperature, year and temperature dependant rate constant type which includes:
  - literature
  - unpublished rate constant, all data: all k
  - unpublished rate constant, all data with parameter Fcent unfixed: float f k
  - unpublished rate constant, all data with parameter m unfixed: float m k
  - unpublished rate constant, only data in 298-258 K included: 298-258 k

Outputs a .feather file for each temperature, rate constant type and is organised in folders named 'rate constant type-year', these are automatically made in the current directorty when the code is run.

**Please note:** The code will only run if the DEFRA measurement data is contained within a folder labelled 'DEFRA' and follow the naming convention of the CSV files uploaded to this repository.

**The Analysis notebook**
Carries out all data analysis and graph production seen in the thesis.

**Please note:** Please leave the previously generated folders and filenames as produced and in the working directory otherwise the code won't run.

**The Analysis Package**
Contains all the functions useed in the LM data formatting and Analysis notebook, all functions have a doc string descriptor accesible using help(analysis.function_name) for extra information on the use, inputs and outputs.

**Please Note:** I have also uploaded a zipped folder containing all the notebooks, packages, DEFRA data and formatted data in the required folders to the Blackboard submission point if that is easier.
  - named: 2113291-A,K-complete code and data