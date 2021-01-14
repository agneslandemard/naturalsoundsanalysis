%% Starter
% Source code for analysis described in:
% Distinct higher-order representations of natural sounds in human and ferret
% auditory cortePlotCorrelationx (2021)
% Landemard A, Bimbard C, Demené C, Shamma S, Norman-Haigneré S, Boubenec Y
% 2021, Jan. 11

%% Define paths 

global data_path additional_path analysis_path code_path

% TO ADAPT before starting
project_directory = 'naturalsoundsdata/'; % define data folder location
% Data can be downloaded here :
% https://lsp.dec.ens.fr/en/research/supporting-materials-848 ("natural
% sounds dataset")

code_path = 'naturalsoundsanalysis/'; % define code folder location

data_path = [project_directory 'fUSData/'];
additional_path = [project_directory 'AdditionalData/'];
analysis_path = [project_directory 'Analysis/'];

addpath(genpath(code_path));

%% Source code for paper main figures

open PlotPaperFigures

%% Denoising procedure
% to implement full denoising procedure: 
% this will lead to processed data used for figures (as in
% natural_myVersion folder)
open FullDenoising

% Code implementing CCA correction :
open CCACorrection

% Simulation for CCA
open SimulationCCA

%% Notice on how to get voxels in their original spatial organisation 
% from raw or denoised data
open GetOriginalVoxels 
