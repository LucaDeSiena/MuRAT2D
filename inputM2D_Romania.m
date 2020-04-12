%% Hard coded data in .m spreadsheet - 2D Analysis
% EVERYTHING MARKED BY 'E' IS TO/CAN BE EDITED

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% GENERAL - CHOICES
% Which analysis do you want to perform?
% Pick delay and Qc without kernels: pa=1
% Pick delay and Qc with analytic kernels: pa=2
% Pick delay and Qc with Pacheco-Snieder kernels: pa=3
Murat.analysis = 2;%'E'

%Treshold to reduce computational time for Pacheco-Snieder kernels.
%It divides the inversion grid by the treshold.
if Murat.analysis == 3
    Murat.geometry.kernelTreshold=10;%'E'
end

% INPUT DATA
% Choose between P- (2) and S- (3) waves for peak delay
% RULE: First column in time files always contains origin of events
% RULE: Second column in time files always contains P-wave phase
% RULE: Third column in time files always contains S-wave phase
Murat.data.PorS = 2; %E

% Central frequency (Hz) - set it according to your spectrograms
Murat.data.centralFrequency = 6;%'E'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PATHS AND FIGURES
% Working directory
Murat.paths.workingdir = './';%'E'

% Folder containing data
Murat.paths.datadir ='./sac_Romania/*.sac';%'E'

% Name of folder to store results and figures
Murat.paths.label = 'Romania';%'E'

% Name of the SAC variables where zero time and P/S pickings are saved
Murat.paths.originTime = 'SAChdr.times.o';%'E'
Murat.paths.PTime = 'SAChdr.times.t0';%'E'
Murat.paths.STime = [];%'E'

%Figure format - 'jpeg' (fast) or 'tiff' (for publication)
Murat.figures.format = 'jpeg';%'E'

% Figure visibility - pick 'on' or 'off'
Murat.figures.visibility = 'on';%'E'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DATA DRIVEN CHOICES - WAVEFORMS
% Work with 1 vertical (1) or 2 horizontal (2) recordings, or with the
% three components(3)
% Order MUST BE: WE, SN, Z.
Murat.data.components =  1;%'E'

% Parameter for smoothing - must be > 2
Murat.data.smoothing = 8;%'E'

% Maximum window to pick pick-delays in seconds
Murat.data.maximumPD = 10;%'E'

% Minimum peak delay considering scattering in the area and frequency
Murat.data.minimumPD = 0.5;%'E'

% Lapse time for the start of the window used to measure and calculate the
% normalization energy
Murat.data.startLT = 30;%'E'

% Total coda window length for Qc and normalization
Murat.data.codaWindow = 30;%'E'

% The sped coefficient sets the spectral energy decay of the coda
% wavefield
Murat.data.spectralDecay = 0.5;%'E'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DATA DRIVEN CHOICES - VELOCITY
%None (0), 1D (1) or 3D (2) velocity model
Murat.geometry.availableVelocity=0;%'E'

% If the origin time is unknown, you can set a theoretichal velocity for
% the whole area and evaluate it from picking. It must be the velocity of
% the phase you are mapping. Velocity is in km/s.
Murat.data.averageVelocityP = 6;%'E'
Murat.data.averageVelocityS = 3;%'E'

%name of the velocity model - if options (1) or (2) 
namev=[];%'E'

% set createrays to 1 if you want to create files for rays,
% if options (1) or (2)
% WARNING: This will create A BIG .mat file
Murat.geometry.createrays = 0;%'E'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% GEOMETRY
% Import event origin time and coords of event and station from files
% (Import=1).
% Ideal is to have both in the SAC header (Import=2) and do it
% in lat/long.
Murat.geometry.import = 2;%'E'

% UTM coordinates of the origin of the model - this must be put
% (1) at least one cell before the westernmost and southernmost
% earthquakes/station for positive latitudes and longitudes
%
% (2) at least one cell after the easternmost and northermost
% earthquakes/station for positive latitudes and longitudes
Murat.geometry.origin=[20 43];%'E'
Murat.geometry.end=[29 48];%'E'

%Step of the grid and number of nodes for pd and Qc imaging
Murat.geometry.gridX = 10;%'E'
Murat.geometry.gridY = 6;%'E'

%Set if in meters (1) or degrees (111)
Murat.geometry.degreesorutm=111;%'E'

%Rays measured in meters or degrees
if Murat.geometry.degreesorutm==1
    %Set it to either 1 (km) or 1000 (meters)
    Murat.geometry.unity = 1000;%E
elseif Murat.geometry.degreesorutm==111
    Murat.geometry.unity = 1;%E
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INVERSION
% Size anomaly for testing: twice(2) or four times node spacing.
% The input of the checkerboard must be  checked visually at the end of
% the process    
Murat.inversion.sizeCheck = 2;%'E'

% Values of attenuation for testing
Murat.inversion.highCheck = 0.02;%'E'
Murat.inversion.lowCheck = 0.001;%'E'

% Seconds of time-windows for non-linear inversion and corresponding
% number, set nonlinear=1 to activate;
Murat.inversion.nonlinear = 0;%'E'

% Uncertainty on Qc estimation
if Murat.inversion.nonlinear==0
    % Minimum R-squared for Qc fitting
    Murat.inversion.fitT = 0.1;%'E'
elseif Murat.inversion.nonlinear==1
    % Length of smaller time windows to compute compute coda intensity
    Murat.inversion.fitL = 3;%'E'
    Murat.inversion.fitT = 5;%'E'
    % Number of time windows to compute coda intensity
    Murat.inversion.fitN=Murat.data.codaWindow/Murat.inversion.fitL;
    %Grid search - set the minimum, maximum and spacing to search in the
    %parameter space
    m1min = 0;%'E'
    Murat.inversion.minimum=m1min; % minimum inverse Qc allowed
    m1max = 0.01;%'E'
    Murat.inversion.maximum=m1max;% minimum inverse Qc allowed
    total = 1001;%'E'
    Murat.inversion.total = total;% total number of Qc in the interval
    Murat.inversion.fit = (m1min + (m1max-m1min)/(total-1)*(0:total-1))';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CHECKS AND LOOPS - DO NOT EDIT

if Murat.geometry.import == 1
    
    %Name of the event file if importing from file
    namee='even.txt';%'E'
    
    %Name of the station file if importing from file
    names='staz.txt';%'E'
    
end

% How big are markers for stations and events
Murat.figures.sizeMarker=60;

% Check that user has the Mapping Toolbox installed.
hasMT = license('test', 'map_toolbox');
Murat.figures.hasMT=hasMT;
if ~hasMT
  % User does not have the toolbox installed.
  message =...
      sprintf('No Mapping Toolbox - Maps will not be geo-localised');
end

if Murat.geometry.import==1
    
    % Opening the event file. Format is:
    % column (1) = twelve numbers for the origin time of the event
    % (date+time in seconds)
    % column (2) = UTM (WE) or latitude
    % column (3) = UTM (SN) or longitude
    % column (4) = Altitude above sea level in meters
    event = fopen(namee);
    namee= textscan(event,'%s %f %f %f');
    nameeven=namee{1};
    even=[namee{2} namee{3} -namee{4}];
    fclose(event);
    
    % Opening the station file. Format is:
    % column (1) = Name of station (3 characters)
    % column (2) = UTM (WE) or latitude
    % column (3) = UTM (SN) or longitude
    % column (4) = Altitude above sea level in meters
    station=fopen(names);
    names= textscan(station,'%s %f %f %f');
    namestation=names{1};
    stat=[names{2} names{3} names{4}];
    fclose(station);
    
end

origin = Murat.geometry.origin;
nxc = Murat.geometry.gridX;
nyc = Murat.geometry.gridY;

x=linspace(Murat.geometry.origin(1),Murat.geometry.end(1),nxc);
y=linspace(Murat.geometry.origin(2),Murat.geometry.end(2),nyc);

Murat.geometry.gridStepX=x(2)-x(1);
Murat.geometry.gridStepY=y(2)-y(1);
stepgx = Murat.geometry.gridStepX;
stepgy = Murat.geometry.gridStepY;


%pre-define 2D matrix in space
XY = zeros(floor(nxc)*floor(nyc),2);


% Some figures require moving the point of plotting to half the resolution
Murat.geometry.x=x+stepgx/2;
Murat.geometry.y=y+stepgy/2;

index=0;
for i=1:nxc
    for j=1:nyc
        index=index+1;
        XY(index,1:2)=[x(i) y(j)];
    end
end
Murat.geometry.map=XY;
lxy=length(XY(:,1));

% Name of the file containing the seismograms and rays
% creates a list of the waveforms and calculates their number.
% RULE: This is necessary for the option Murat.geometry.import = 1, where
% event and station info are in two separate files. In this case,
% the SAC file-name must have a specific format:
% the first 12 characters are the origin time of the event, as per
% event-location file). Characters 18:20 are the name of station as per
% station file.
%
% For Murat.geometry.import = 2, there is no requirement on SAC filenames.
% In the case of pa=3 the names of rays created will be identical to files.

% Get name of files in traces' folder
DFolder=Murat.paths.datadir;
list = dir(DFolder);
filenames = {list.name}'; %create cell array of file names
filedir = {list.folder}'; %create cell array of file folder
listasac = filenames;
lls = length(listasac);

evestaz=zeros(lls,6);
indexray=0;
for i = 1:lls
    listasac{i,1} = cat(2,filedir{i},'/',filenames{i});
    li = filenames{i,1};
    li1 = cat(2,li(1:12),li(18:20));
        
    %Here it takes the event/station name. It is necessary to adapt the
    %numbers to where even name and station name are
    if Murat.geometry.import==1
        for ii=1:length(nameeven)
            namee1=nameeven{ii};
            namee=namee1(1:12);
            
            % Loop over stations in staz.txt file
            for ir=1:length(namestation)
                namest1=namestation{ir};
                namest=namest1(1:3);
                levst=cat(2,namee,namest);
                
                if find(strncmp(li1,levst,length(levst)))>0
                    indexray=indexray+1;
                    sst = [even(ii,:) stat(ir,:)];
                    evestaz(indexray,1:6)=sst;% Creates file of sel. coords
                    break
                end
            end
        end
    end
end
Murat.paths.listasac=listasac;
Murat.geometry.evestaz=evestaz;

if exist(cat(2,'./',Murat.paths.label),'dir')~=7
    mkdir(Murat.paths.label)
end

if Murat.geometry.availableVelocity==1
    
    modv1=load(namev); %1D velocity model iasp91
    dend=-34000; %select maximum depth range
    
    li=length(modv1(:,1));
    modv=zeros(lxy*li,5);
    index=0;
    for i=1:nxc
        for j=1:nyc
            index1=index;
            index=index+1;
            modv(index1*li+1:index1*li+li,1)=XY(index,1)*1000;
            modv(index1*li+1:index1*li+li,2)=XY(index,2)*1000;
            modv(index1*li+1:index1*li+li,3)=...
                -modv1(1:li,1)*1000;
            modv(index1*li+1:index1*li+li,4)=...
                modv1(1:li,PorS+1);
        end
    end
    
    modv(:,1)=(modv(:,1)-modv(1,1));
    modv(:,2)=(modv(:,2)-modv(1,2));
    modv(modv(:,3)<dend,:)=[];
    
elseif Murat.geometry.availableVelocity==2
    
    % This works in [x,y,z], created for UTM WGS84 and origins to zero.
    modv=load(namev); %3D velocity model from text file
    modv(:,5)=0;
    modv(:,1)=modv(:,1)+origin(1);
    modv(:,2)=modv(:,2)+origin(2);
    
    
    Murat.geometry.modv=modv;
    
    chx=find(modv(:,1)~=modv(1,1),1);
    chy=find(modv(:,2)~=modv(1,2),1);
    resol2x = abs(modv(chx,1)-modv(1,1))/2;
    resol2y = abs(modv(chy,2)-modv(1,2))/2;
    resol2z = abs(modv(2,3)-modv(1,3))/2;
    
    %Steps of the velocity models
    passox=max(modv(:,1))-min(modv(:,1));
    passoy=max(modv(:,2))-min(modv(:,2));
    passoz=max(modv(:,3))-min(modv(:,3));
    
    passo=[passox passoy passoz];
    resol=[resol2x resol2y resol2z];
    resol2=min(resol);
    Murat.geometry.resolutionMin=resol2;
    
    %Creates grid for direct waves and check for zeroes
    %Regular step of the gridD for interpolation - half of step of modv
    
    ixD=floor(passox/resol2x);%numer of x,given the step of gridD
    iyD=floor(passoy/resol2y);%numer of y,given the step of gridD
    izD=floor(passoz/resol2z);%numer of depths,given the step of gridD
    
    gridD=zeros(3,max(passo./resol));
    gridD(1,1:ixD)=origin(1):resol2x:origin(1)+passox-resol2x;
    gridD(2,1:iyD)=origin(2):resol2y:origin(2)+passoy-resol2y;
    gridD(3,1:izD)=-modv(1,3):resol2z:-modv(1,3)+passoz-resol2z;
    Murat.geometry.gridD=gridD;
    % gridD dimensions
    pvel = zeros(ixD,iyD,izD);
    
    %NUMBER OF X, Y, AND Z LAYERS
    for k=1:izD
        index=0;
        for i=1:ixD
            for j=1:iyD
                index=index+1;
                pvel(i,j,k) = modv(index,4);
            end
        end
    end
    Murat.geometry.pvel=pvel;
end

%clearvars -except Murat
save('Murat.mat','Murat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%