
# From Stomach to Striatum: Ghrelin Increases Work for Rewards

This project investigates the influence of intravenous infusions of the stomach-derived peptide acyl ghrelin on human motivation and its role as a neuromodulator - affecting functional connectivity and dopamine responses, using simultanous [C11]raclopride PET-fMRI during infusions of ghrelin vs. saline. 

This project code accompanies [this paper]
If you use code or data from this repository please cite the original publication. 


## Authors
- Corinna Schulz 
- Franziska Peglow
- Christian la Fougère
- Benjamin Bender
- Johannes Klaus
- Sabine Ellinger
- Martin Walter
- Gerald Reischl 
- Matthias Reimold 
- Nils B. Kroemer

## Table of Contents

- [Abstract](#abstract)
- [Getting Started](#getting_started)
- [Installing](#installing)
- [Code Structure](#code_structure)

## Abstract: 

Preclinical evidence demonstrates that gut signals influence motivated behavior through the dopaminergic system. However, translational research in humans is scarce, and it is not known if surges of the stomach-derived hormone ghrelin acutely boost motivation via dopamine transmission. To close this gap, we investigated the effects of acyl ghrelin infusions (vs. saline) on instrumental work for rewards, brain responses, and dopamine transmission in healthy adults in a double-blind, randomized and placebo-controlled crossover study with [11C]raclopride PET-fMRI. In line with preclinical findings, ghrelin enhanced hunger sensations (p = .007), and effort for food and small rewards. Crucially, ghrelin enhanced functional connectivity between the hypothalamus and the striatum (NAcc: pSVC = .016), as well as within the striatum (NAcc<>Putamen: pSVC = .001), pointing to an enhanced coupling between homeostatic and motivational circuits. However, ghrelin did not alter dopamine binding potential (reflecting changes in endogenous tone), but boosted task-free phasic pulses of the NAcc BOLD signal (p = .020), which were negatively associated with binding potential (p = .034). During the task, ghrelin increased cue- (precentral gyrus: pFWE = .003) and effort-evoked brain responses (e.g., VTA/SN: psvc = .002) in accordance with increases in motivation. Together, these findings recapitulate preclinical work that ghrelin primarily boosts phasic dopamine release, reflecting enhanced incentive salience of food and smaller rewards. This highlights the potential of targeting gut-brain interactions to improve motivational symptoms, such as loss of appetite or anhedonia.
 
 ![Graphical Abstract](/publication_figures/Abstract_FIN.jpg)

## Getting_Started

### Clone the project using SSH from github 

`$ git clone git@github.com:neuromadlab/Schulz_Stomach_Striatum.git`

Read more about cloning a repository [here](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository) 

### Installing 

Code is written in `R V.4.2.2.`, within `VSCode` using Radian. 
To set up R in VSCode follow these [instructions] (https://code.visualstudio.com/docs/languages/r). [Better R in VSCode](https://schiff.co.nz/blog/r-and-vscode/). 
If using Radian from windows, set the R home env in the terminal.integrated.env

`RStudio` should work equally fine. 

All required libraries are installed and loaded at the beginning of each script using [librarian] (https://cran.r-project.org/web/packages/librarian/vignettes/intro-to-librarian.html) 

## Code_Structure 

### Abbreviations: 
* TUE008 NIMG: lab internal study identifier, Neuroimaging part

### Input Folder
* all relevant data that is used for plotting (except imaging plots created with Mango) and for statistical analysis on the group level (except voxelwise imaging analysis and FC analyses (CONN, SPM).) 

#### Code Folder: 

* `0.Sideeffects.r` Inspection of subjective ratings of infusion (word clouds), or debriefing at the end of the session (sideeffects; Chi-Square)

* `blood_preprocessing.R` Preprocessing of the blood data collected during the study (3 Timepoints) + quality checks + residualization for further analyses 

* `1.VAS_blood.r` Visualization and analysis of visual analog scale data (e.g. reports of hunger and mood), and plotting of ghrelin values 

* `2.IMT_PET.r` Visualization and analysis of the instrumental motivation task (IMT). Model comparison of ghrelins effect on motivation. Several sensitivity analyses/checks. Stats and plots for ghrelin and dopamine binding potential. Model comparisons. 

* `3.FC.r` ROI-based FC analyses (checking voxel based results) and then explored association with behavior, PET. 

* `4.PulseDetection.r` Visualization and analysis ghrelin-induced changes in pulse count and pulse matches and their association with dopamine binding potential. Sensitivity checks for different parameters. 


#### Output Folder: 

Figures are saved here. Output of LMEs as reported in the papers is saved as .csv or .html.   


## Citation

If you use this code or data, please cite the original publication.


## Contact

For questions or issues, please contact the corresponding authors or open an issue in this repository.
