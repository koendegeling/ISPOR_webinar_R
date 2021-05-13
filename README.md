# ISPOR Student Network Webinar: Health Economic Modelling in R

*Repository for the ISPOR Student Network Webinar on using R for Health Economic Modelling by Petros Pechlivanoglou and Koen Degeling*

## Information

This repository includes the files that are used in the webinar to demonstrate the implementation of a three-state Markov model in R to compare a standard of care (SoC) treatment to an experimental treatment (Exp).
The health states used are *ProgressionFree* (F), *Progression* (P) and *Dead* (D), which is a model structure that has been frequently used to evaluate systemic treatment for advanced cancers.
Therefore, this model is often referred to as the three-state oncology model.
Compared to the standard of care, the experimental treatment reduces the risk of progression, but it is also more expensive.
A cost-utility analysis is performed with the Markov model to assess whether the increase in quality-adjusted life years justify the increase in costs.
A stpe-by-step implementation of the Markov is provided in the `oncologyMarkov_illustration.R` script, which also includes a demonstration of deterministic sensitivity analyses and a probabilistic analysis of the model.
For these analayses, the Markov has also been implemented in a function that can easily be used to evaluate the model for a certain set of parameter values.
This function is defined in the `oncologyMarkov_function.R` script, which is sourced (i.e. run) from within the `oncologyMarkov_illustration.R` script.

Besides the Markov model, an illustration is provided for a semi-Markov model with time-dependent transition probabilities for the transition from the *ProgressionFree* to the *Progression* state.
The `oncologySemiMarkov_illustration.R` script demonstrates a survival analysis for hypothetical time-to-progression (TTP) data in the `df_TTP.RData` data object and how this can be implemented in the semi-Markov model.
As for the Markov model, the `oncologySemiMarkov_function.R` script defines a functions that can be used to evaluate the semi-Markov model for deterministic sensitivity analyses and a probabilistic analysis.


## Using the files

In general, there are two ways to get the files from this repository onto your computer:
  1)  Downloading the repository as a .zip file and unpacking the compressed folder somewhere on your computer 
  2)  Cloning the repository on git on your computer

Regardles of how you get the files onto your computes, it is recommended to open the 'ISPOR_webinar_R.Rproj' project file rather than individual scripts to assure the working directory is automatically set to the right folder.


The following files are included in this repository:

| File                                | Description                                                                  |
| ----------------------------------- | ---------------------------------------------------------------------------- |
| `df_TTP.RData`                      | Hypothetical data used in the illustration of the semi-Markov model          |
| `ISPOR_webinar_R.Rproj`             | R Studio project file (use to ensure working directory is correct)           |
| `oncologyMarkov_function.R`         | Function wrapping the Markov model used in the Markov illustration           |
| `oncologyMarkov_illustration.R`     | Script illustrating the implementation and analysis of the Markov model      |
| `oncologySemiMarkov_function.R`     | Function wrapping the semi-Markov model used in the semi-Markov illustration |
| `oncologySemiMarkov_illustration.R` | Script illustrating the implementation and analysis of the semi-Markov model |


## Acknowledgement of the DARTH workgroup

The scripts are heavily based on the scripts that have been developed by the Decision Analysis in R for Technologies in Health (DARTH) workgroup, which consists of:
  - Fernando Alarid-Escudero, PhD (1)
  - Eva A. Enns, MS, PhD (2)
  - M.G. Myriam Hunink, MD, PhD (3,4)
  - Hawre J. Jalal, MD, PhD (5)
  - Eline M. Krijkamp, MSc (3)
  - Petros Pechlivanoglou, PhD (6,7)
  - Alan Yang, MSc (7)

  1. Division of Public Administration, Center for Research and Teaching in Economics (CIDE), Aguascalientes, Mexico
  2. University of Minnesota School of Public Health, Minneapolis, MN, USA
  3. Erasmus MC, Rotterdam, The Netherlands
  4. Harvard T.H. Chan School of Public Health, Boston, USA
  5. University of Pittsburgh Graduate School of Public Health, Pittsburgh, PA, USA
  6. University of Toronto, Toronto ON, Canada
  7. The Hospital for Sick Children, Toronto ON, Canada

Please cite the applicable DARTH publications when using any of their code:
  - Jalal H, Pechlivanoglou P, Krijkamp E, Alarid-Escudero F, Enns E, Hunink MG. An Overview of R in Health Decision Sciences. Med Decis Making. 2017; 37(3): 735-746. https://journals.sagepub.com/doi/abs/10.1177/0272989X16686559
  - Alarid-Escudero F, Krijkamp EM, Enns EA, Yang A, Hunink MGM Pechlivanoglou P, Jalal H. Cohort State-Transition Models in R: A Tutorial. arXiv:200107824v2. 2020:1-48. http://arxiv.org/abs/2001.07824
  - Krijkamp EM, Alarid-Escudero F, Enns EA, Jalal HJ, Hunink MGM, Pechlivanoglou P. Microsimulation modeling for health decision sciences using R: A tutorial. Med Decis Making. 2018;38(3):400–22. https://journals.sagepub.com/doi/abs/10.1177/0272989X18754513
  - Krijkamp EM, Alarid-Escudero F, Enns E, Pechlivanoglou P, Hunink MM, Jalal H. A Multidimensional Array Representation of State-Transition Model Dynamics. Med Decis Making. 2020 Feb;40(2):242-248. https://journals.sagepub.com/doi/10.1177/0272989X19893973
