% PROGRAM FOR: 2D peak delay and coda-wave attenuation tomography with
% kernels
%
% Author:  L. De Siena, December 2019
%
% Installation
% ------------
% SYSTEM: The program works on 1) a Macbook pro with High Sierra and 2)
% Ubuntu 16.04.6 LTS, both with Matlab R2017a.
%
% Necessary Toolboxes: Signal Processing, Curve Fitting,
% Image Processing (compulsory) and Mapping (optional, for geolocalisation).
%
% Two sample datasets (Mount St. Helens and Romania) can be downloaded
% at https://doi.pangaea.de/10.1594/PANGAEA.893893 to test the code.

% INSTRUCTIONS:
% The current version works following these steps:
%
% 1. Download the package at https://github.com/LucaDeSiena/MuRAT.
%
% 2. Download the two sample datasets at
% https://doi.pangaea.de/10.1594/PANGAEA.893893.
% Unzip the MSH (Mount St. Helens) and Romania datasets and put
% the folders in the Murat-master folder.
%
% 3. Build your own input file (.m) - each field is described in the
% attached README file and in the INPUT sections of the code.
% When building your example, use one of the Input
% files as template (MSH for 3D and Romania for 2D). Always start with an
% analytic 2D analysis (“pa=2”).
%
% 4. Run MuRAT2_2.m.
%
% 5. Select the name of the input file desired
% (Input_MSH.m, Input_Romania.m, Input_Pollino and Input_OlDoinyo).
%
% 6. After the L curves are produced, write the smoothing parameter.
%
% Theory:
% ----
% Peak-delays and 2D coda attenuation: Takahashi et al. 2007 (JGR);
% Calvet et al. 2013b (GJI); De Siena et al. 2016 (EPSL);
% Borleanu et al. 2017 (Tectonophysics)
%
% Kernel-based 2D coda attenuation: Del Pezzo et al. 2016 (GJI);
% De Siena et al. 2017 (GRL); Napolitano et al. 2019 (Geos. Fron.)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUTS

addpath('./Utilities_Matlab')

clear; close all; clc

[file,path]                         =   uigetfile('*.m');

if isequal(file,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(path,file)]);
end

run(fullfile(path, file))
%% CHECKS

disp('Checks and Loops')

Murat                               =   Murat_checks(Murat);

%% Seismic attributes for peak delay and Qc imaging

disp('Data Section')

Murat                               =   Murat_data(Murat);

%%  2D peak-delay and Qc TOMOGRAPHIC INVERSIONS

disp('Inversion Section')

Murat                               =   Murat_inversion(Murat);

%% Creating maps

disp('Plot Section')

Murat                               =   Murat_plot(Murat);

save('Murat.mat','Murat');
