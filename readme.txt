REMEDID: Retrospective Methodology to Estimate Daily Infections from Deaths 

The REMEDID algorithm allows reconstructing daily infections of COVID-19 (or any infectious disease) from:
(1) Daily deaths from that disease
(2) The probability density functions (PDF) discretised in days of the incubation period (IP) and the period from onset of symptoms to death (IOD)
(3) Case fatality rate (CFR).

To calculate the REMEDID infections for a PI PDF and an IOD PDF, use the Matlab file: remedid.m

If the coefficients defining the PDFs are given by confidence intervals (assumed here to be uniformly distributed), REMEDID infections can be calculated with the associated errors. For this purpose, a Monte Carlo simulation is performed with the Matlab file: remedid_monte_carlo_error.m. This file uses two auxiliary files: remedid_pdf.m and pdf_infection2death_monte_carlo_error.m.



Chile and Santiago de Chile Metropolitan Region (CRM):
COVID data for Chile and Santiago de Chile Metropolitan Region used and calculated in Márquez et al. (2024) are found in the Excel files:
- Datos_COVID_Santiago_Region_Metropolitana_2020.xlsx


Citation
If the remedid.m function is used, it should be cited:
García-García, D., I. Vigo, E. S. Fonfría, Z. Herrador, M. Navarro, and

If the rest of the functions or data are used, García-García et al. (2021) should be cited and: 