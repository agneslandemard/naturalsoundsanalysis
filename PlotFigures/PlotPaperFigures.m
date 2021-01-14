% Plots main figures (and some sup) from: 
% Distinct higher-order representations of natural sounds in human and ferret
% auditory cortex (2021)
% Landemard A, Bimbard C, Demené C, Shamma S, Norman-Haigneré S, Boubenec Y
% 2021, Jan. 14

% Paths need to be defined using Starter.m first

% denoised version that will be used to plot figures
version_name = 'myVersion';
NbICs = 8;

%% Plot example single voxels (Fig 1 & 2)
% Raw and denoised time-courses for example single voxels
% Nat vs Synth and test-retest time-averaged responses
% Viewing of voxels' spatial position

PlotSingleVoxels

%% Plot maps of NSE natural vs synthetic for all models, for exp I (Fig 2 E-F)
% quantification of response difference using NSE
% Top view NSE maps 
% NSE as function of distance to center of PAC

PlotNSE_NatSynth

% Human figures from Fig 2 are drawn from 
% Norman-Haignere et al., 2018

%% Plot example components for exp I (Fig 3)
% 3 components
PlotFerretComponents

%% Plot vocalization results , exp I & II (Fig 4)
% Movement (panel C)
PlotMvt

% for exp I (panels E,F)
% Maps of difference nat - synth for ferret and humans
% NSE as function to center of PAC, by category of sounds, ferrets and
% humans
% Summary boxplot
PlotHumanData

% for exp II (panels D,G)
% NSE maps
% Maps of difference nat - synth
% NSE as function to center of PAC, by category of sounds
PlotVocalizationsResults

%% Effect of denoising
% fig S2
MeasureDenoisingEffect


