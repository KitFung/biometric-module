%% Data %%

% In this folder, you will find two databases (matlab workspaces): 

% data1.mat is the in silico data from DREAM5 challenge 4. DataTest is simulated data to % get a hold of the algorithm.
% For more information and to download further datasets, please visit the DREAM website: 
% http://wiki.c2b2.columbia.edu/dream/index.php/The_DREAM_Project

% Also used in the paper is the e.Coli data from Faith et al., 2007. 
% Below are the versions that can be downloaded. 
% The expression data can be downloaded from M^3D (http://m3d.bu.edu/) version 
% 6 build 4.
% The interactions can be downloaded from RegulonDB 
% (http://regulondb.ccg.unam.mx/), version 7.2

% All of these are structures will the following fields: 
% - expdata: the p experiments by n genes expression matrix
% - tflist: a list of the transcription factors (TF)
% - genenames: a list containing all gene names
% - tf_index: the indexes of the TFs in genenames
% - gold_edges_num: the gold_edges 

% Additionnally, DREAM5 datasets should contain:
% - pdf_aupr: the pdf used to compute the p-values for the AUPR
% - pdf_auroc: same for AUROC
