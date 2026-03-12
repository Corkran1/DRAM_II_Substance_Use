# DRAM_II_Substance_Use
MATLAB code to run the DRAM algorithm used to generate our mean transmission rates per substance use nursing unit pair. 

This GitHub repository provides the materials needed to implement the Delayed Rejection Adaptive Metropolis (DRAM) algorithm for estimating parameters in a heterogeneous infectious disease transmission model.

The model distinguishes between substance-using and non-substance-using patient populations within hospital nursing units, allowing transmission dynamics to differ across these groups.

The repository is organized into three main sections:

Section 1.0: How this Github Page is Organized

Section 2.0: Paremter Dictonary for DRAM 

This section provides definitions and descriptions of the parameters estimated by the DRAM algorithm. These parameters represent key components of the transmission process, including:

[I] Transmission rates

[II] Admission rates

[III] Recovery rates

[IV] Discharge rates

Because the model is heterogeneous, some parameters correspond specifically to substance-using patients, while others correspond to non-substance-using patients.

Refer to this section when interpreting the parameter estimates produced by the MATLAB code.

Section 3.0: Instructions for Running  DRAM 

Step 1

Navigate to the folder, within section_3.0_codes, then find corresponding to the specific nursing unit (substance-using amd non-substance-using populations are incldued in file ) that you wish to analyze.

Step 2

Copy all files from that folder and place them into a new working folder inside your MATLAB directory.

Step 3

Locate the MATLAB .m file titled:

“NursingUnitName_Main_Code.m”

Open this file in MATLAB.

Step 4

Update the following weekly totals in the script:

Weekly susceptible population

Weekly infected admissions

Weekly susceptible admissions

Weekly discharges

Weekly infected individuals

Weekly recovered individuals

For the results presented in our paper, these values were derived from electronic health record (EHR) data. However, users are not required to use the same data source. If you are applying this model to a different hospital system or transmission setting, you may substitute your own observed or simulated values.

Step 5

Run the MATLAB script (e.g., NURSING_UNIT_NAME1). After execution, examine the MATLAB workspace, where the DRAM algorithm will output the estimated parameter values.

The workspace will contain posterior mean estimates for each parameter.

Step 6

Match each estimated parameter with the corresponding description listed in Section 2.0 to interpret the results. These parameters include estimates for:

Transmission rates

Discharge rates

Recovery rates

Admission rates

These estimates quantify how transmission dynamics vary across nursing units and between substance-using and non-substance-using patient populations.
