%% Hard coded data in .m spreadsheet - 2D Analysis
% EVERYTHING MARKED BY 'E' IS TO/CAN BE EDITED

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% GENERAL - CHOICES
% Which analysis do you want to perform?
% Pick delay and Qc without kernels: pa=1 - De Siena et al. 2016 EPSL
% Pick delay and Qc with analytic kernels: pa=2 - De Siena et al. 2017 GRL
% Pick delay and Qc with Pacheco-Snieder kernels: pa=3 - Sketsiou et al.
% 2020

Murat.analysis                              =  2;%'E'

%Treshold to reduce computational time for Pacheco-Snieder kernels.
%It divides the inversion grid by the treshold.
if Murat.analysis == 3
    Murat.geometry.kernelTreshold           =  2;%'E'
end

% INPUT DATA
% Choose between P- (2) and S- (3) waves for peak delay
% RULE: First column in time files always contains origin of events
% RULE: Second column in time files always contains P-wave phase
% RULE: Third column in time files always contains S-wave phase
Murat.data.PorS                             =  2; %'E'

% Central frequency (Hz) - set it according to your spectrograms
Murat.data.centralFrequency                 =  12;%'E'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PATHS AND FIGURES
% Working directory
Murat.paths.workingdir                      =  './';%'E'

% Folder containing data
Murat.paths.datadir                         =  './sac_OlDoinyo/*Z.sac';%'E'

% Name of folder to store results and figures
Murat.paths.label                           =  'OlDoinyo';%'E'

% Name of the SAC variables where zero time and P/S pickings are saved
Murat.paths.originTime                      =  'SAChdr.times.o';%'E'
Murat.paths.PTime                           =  'SAChdr.times.t0';%'E'
Murat.paths.STime                           =  'SAChdr.times.t1';%'E'

%Figure format - 'jpeg' (fast) or 'tiff' (for publication)
Murat.figures.format                        =  'jpeg';%'E'

% Figure visibility - pick 'on' or 'off'
Murat.figures.visibility                    =  'on';%'E'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DATA DRIVEN CHOICES - WAVEFORMS
% Work with 1 vertical (1) or 2 horizontal (2) recordings, or with the
% three components(3)
% Order MUST BE: WE, SN, Z.
Murat.data.components                       =  1;%'E'

% Parameter for smoothing - must be > 2
Murat.data.smoothing                        =  8;%'E'

% Maximum window to pick pick-delays in seconds
Murat.data.maximumPD                        =  15;%'E'

% Minimum peak delay considering scattering in the area and frequency
Murat.data.minimumPD                        =  0.5;%'E'

% Lapse time for the start of the window used to measure and calculate the
% normalization energy
Murat.data.startLT                          =  [];%'E'

% Total coda window length for Qc and normalization
Murat.data.codaWindow                       =  20;%'E'

% The sped coefficient sets the spectral energy decay of the coda
% wavefield
Murat.data.spectralDecay                    =  0.5;%'E'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DATA DRIVEN CHOICES - VELOCITY
%None (0), 1D (1) or 3D (2) velocity model
Murat.geometry.availableVelocity            =  0;%'E'

% If the origin time is unknown, you can set a theoretichal velocity for
% the whole area and evaluate it from picking. It must be the velocity of
% the phase you are mapping. Velocity is in km/s.
Murat.data.averageVelocityP                 =  6;%'E'
Murat.data.averageVelocityS                 =  3;%'E'

%name of the velocity model - if options (1) or (2) 
Murat.data.namev                            =  'modv.txt';%'E'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GEOMETRY
% Import event origin time and coords of event and station from files
% (Import=1).
% Ideal is to have both in the SAC header (Import=2) and do it
% in lat/long.
Murat.geometry.import                       =  2;%'E'

%set createrays to 1 if you want to create files for rays
%WARNING: This will create A BIG .mat file
Murat.geometry.createrays                   =  0;%'E'

% UTM coordinates of the origin of the model - this must be put
% (1) at least one cell before the westernmost and southernmost
% earthquakes/station for positive latitudes and longitudes
%
% (2) at least one cell after the easternmost and northermost
% earthquakes/station for positive latitudes and longitudes
Murat.geometry.origin                       =  [35.7 -3];%'E'
Murat.geometry.end                          =  [36.3 -2.4];%'E'

%Step of the grid and number of nodes for pd and Qc imaging
Murat.geometry.gridX                        =  24;%'E'
Murat.geometry.gridY                        =  24;%'E'

%Set if in meters (1) or degrees (111)
Murat.geometry.degreesorutm                 =  111;%'E' 

%Rays measured in meters or degrees
if Murat.geometry.degreesorutm == 1
    %Set it to either 1 (km) or 1000 (meters)
    Murat.geometry.unity                    =  1000;%E
elseif Murat.geometry.degreesorutm == 111
    Murat.geometry.unity                    =  1;%E
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INVERSION
% Size anomaly for testing: twice(2) or four times node spacing.
% The input of the checkerboard must be  checked visually at the end of
% the process    
Murat.inversion.sizeCheck                   =  2;%'E'

% Values of attenuation for testing
Murat.inversion.highCheck                   =  0.02;%'E'
Murat.inversion.lowCheck                    =  0.001;%'E'

% Seconds of time-windows for non-linear inversion and corresponding
% number, set nonlinear=1 to activate;
Murat.inversion.nonlinear                   =  0;%'E'

% Uncertainty on Qc estimation
if Murat.inversion.nonlinear == 0
    % Minimum R-squared for Qc fitting
    Murat.inversion.fitT                    =  0.1;%'E'
elseif Murat.inversion.nonlinear == 1
    
    % Length of smaller time windows to compute compute coda intensity
    Murat.inversion.fitL                    =  3;%'E'
    Murat.inversion.fitT                    =  5;%'E'
    
    % Number of time windows to compute coda intensity
    Murat.inversion.fitN    =   Murat.data.codaWindow/Murat.inversion.fitL;
    %Grid search - set the minimum, maximum and spacing to search in the
    %parameter space
    m1min                                   =  0;%'E'
    m1max                                   =  0.01;%'E'
    
    Murat.inversion.minimum =   m1min; % minimum inverse Qc allowed
    Murat.inversion.maximum =   m1max;% minimum inverse Qc allowed
    total                                   =  1001;%'E'
    Murat.inversion.total   = total;% total number of Qc in the interval
    Murat.inversion.fit     = (m1min+(m1max-m1min)/(total-1)*(0:total-1))';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%