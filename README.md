This repository contains the orignal DEFRA measurement data for London Marylebone for the years 2005, 2010 and 2020.

**The LM data formatting Notebook**
Processes the orignal data and calculates the OH concentrations for each temperature, year and temperature dependant rate constant type which includes:
  - literature
  - unpublished rate constant, all data: all k
  - unpublished rate constant, all data with parameter Fcent unfixed: float f k
  - unpublished rate constant, all data with parameter m unfixed: float m k
  - unpublished rate constant, only data in 298-258 K included: 298-258 k

Outputs a .feather file for each temperature, rate constant type and is organised in folders named 'rate constant type-year', these are automatically made in the current directorty when the code is run.

**Please note:** The code will only run if the DEFRA measurement data is contained within a folder labelled 'DEFRA' and follow the naming convention of the CSV files uploaded to this repository ('LMyear.csv')

**The Analysis Notebook**
Carries out all data analysis and graph production seen in the thesis.

**Please note:** Please leave the previously generated folders and filenames as produced and in the working directory to allow the code to run correctly.

**The Temperature Dependant Rate Constant EXCEL Spreadsheet**
Contains all the temperature dependant rate constant calculations for the OH-NO reaction, including the literature and unpublished parameters. The 'Final data' sheet is read into the Jupyter Notebooks. The analysis of the NIST experimental and review OH-NO rate constants is also contained in a separate sheet.

**The Rate Constant EXCEL Spreadsheet**
Contains all the literature OH-species rate constant values and species molecular masses used in the calculations (Sheet 1). Sheet 2 contains the calculation of the O3 photolysis rate from the MCM parameters.

**The Analysis Package**
Contains all the user-defined functions useed in the LM data formatting and Analysis notebook, all functions have a doc string descriptor accesible using help(analysis.function_name) for extra information on the use, inputs and outputs.

**Please Note:** I have also uploaded a zipped folder containing all the notebooks, packages, DEFRA data and formatted data in the required folders to the Blackboard submission point if that is easier.
  - named: 2113291-A,K-complete code and data
