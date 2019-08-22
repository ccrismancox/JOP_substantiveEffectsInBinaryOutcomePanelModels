# Replication instructions for "Estimating Substantive Effects in Binary Outcome Panel Models: A Comparison"
# Casey Crisman-Cox 
# August 22, 2019

## A note for replicators
There is one file missing from this archive.  The data file `goldmanData.rdata` contains the 2008 National Annenberg Election Survey (NAES) Internet Panel data used to replicate Goldman (2016). However, the terms and conditions of the NAES prohibit me from publicly distributing, displaying, or posting this data to any website. I am allowed to use and disclose the data for academic research and will share it with university-affiliated academic researcher via email (`ccrismancox@gmail.com`). The data can also be requested here: ` https://www.annenbergpublicpolicycenter.org/data-access/`

## Replication package contents
Contents:

- Background information
    - `readme.md`: plain text readme
- Data and helper functions
    - `sparseFirth.R`: `R` code for a sparse version of the penalized method.  Useful for producing the results in  Appendix B.
	- `disputes1.dta`: Replication data for Green, Kim, and Yoon (2001).
	- `temp-micro.dta`: Replication data for Escriba-Folch, Meseguer, and Wright (2018).
- Monte Carlos 
    - `mainMonteCarlos.r`: Replication code for the main Monte Carlos and those in Appendices A and B.  Outputs: `MC_main.rdata`, `Figure1.pdf`, `Figure2.pdf`, `Figure3.pdf`, `Figure4.pdf`, `FigureA1.pdf`, `FigureA2.pdf`, and `FigureA3.pdf`. Prints Tables B1-2 to the console.
	- `appendixC.r`: Replication code for Appendix C. Outputs: `MC_expanded.rdata`, `FigureC1.pdf`, `FigureC2.pdf`, and `FigureC3.pdf`.
	- `appendixD1.r`: Replication code for the first part of Appendix D. Outputs: `MCset_relevant0.rdata`, `FigureD1.pdf`, `FigureD2.pdf`, `FigureD3.pdf`, and `FigureD4.pdf`.
	- `appendixD2.r`: Replication code for the second part of Appendix D. Outputs: `MCset_irrelevant0s.rdata`, `FigureD5.pdf`, `FigureD6.pdf`, and `FigureD7.pdf`.
	- `appendixE.r`: Replication code for Appendix E. Outputs: `MC_highlyCensored.rdata`. Prints Table E.1 to the console.
    - `appendixF.r`: Replication code for Appendix F. Outputs: `MC_nonlinear.rdata`, `FigureF1.pdf`, `FigureF2.pdf`, and `FigureF3.pdf`.
    - `appendixG.r`: Replication code for Appendix G. Outputs: `MC_within.rdata`, `FigureG1.pdf`, `FigureG2.pdf`, and `FigureG3.pdf`.
	- `appendixH.r`: Replication code for Appendix H. Outputs: `MC_pml.rdata`. Prints Tables H1-5.
- Estimation and analysis
    - `emw_JOP_Revision.r`: Replication code for the EMW example (main text and Appendix J). Outputs: `Figure5.pdf` and `FigureJ1.pdf`. Prints Table 1 to console
	- `genderExample_JOP_Revision.R`: Replication for the Goldman (2018) example. Outputs: `Figure6.pdf`. Prints Table 2 to console.
	- `dirtypool_JOP_Revision.r`: Replication for the GKY example. Outputs: `Figure7.pdf`. Prints Table 3 to console.

