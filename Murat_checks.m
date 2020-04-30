%% Checks and loops necessary to set up the Murat structure
%The function takes the input and defines a seies of parameters o be used
%in the following functions. All will be saved inside the Murat variable.

function Murat=Murat_checks(Murat)
% CHECKS AND LOOPS

if Murat.geometry.import == 1
    
    %Name of the event file if importing event locations from file
    namee                                   =  'even.txt';
    
    %Name of the station file if importing event locations from file
    names                                   =  'staz.txt';
    
end

% How big are markers for stations and events
Murat.figures.sizeMarker    = 60;

% Check that user has the Mapping Toolbox installed.
hasMT = license('test', 'map_toolbox');
Murat.figures.hasMT         = hasMT;
if ~hasMT
  % User does not have the toolbox installed.
  sprintf('No Mapping Toolbox - Maps will not be geo-localised');
end

if Murat.geometry.import == 1
    
    % Opening the event file. Format is:
    % column (1) = twelve numbers for the origin time of the event
    % (date+time in seconds)
    % column (2) = UTM (WE) or latitude
    % column (3) = UTM (SN) or longitude
    % column (4) = Altitude above sea level in meters
    event                   = fopen(namee);
    namee                   = textscan(event,'%s %f %f %f');
    nameeven                = namee{1};
    even                    = [namee{2} namee{3} -namee{4}];
    fclose(event);
    
    % Opening the station file. Format is:
    % column (1) = Name of station (3 characters)
    % column (2) = UTM (WE) or latitude
    % column (3) = UTM (SN) or longitude
    % column (4) = Altitude above sea level in meters
    station                 = fopen(names);
    names                   = textscan(station,'%s %f %f %f');
    namestation             = names{1};
    stat                    = [names{2} names{3} names{4}];
    fclose(station);
    
end

%Loading the origin and steps of the grid to create plotting variables and
%define the steps of the grid.
origin                      = Murat.geometry.origin;
nxc                         = Murat.geometry.gridX;
nyc                         = Murat.geometry.gridY;

x           = linspace(Murat.geometry.origin(1),Murat.geometry.end(1),nxc);
y           = linspace(Murat.geometry.origin(2),Murat.geometry.end(2),nyc);

Murat.geometry.gridStepX    = x(2)-x(1);
Murat.geometry.gridStepY    = y(2)-y(1);

stepgx                      = Murat.geometry.gridStepX;
stepgy                      = Murat.geometry.gridStepY;

%pre-define 2D matrix in space
XY                          = zeros(floor(nxc)*floor(nyc),2);


% Some figures require moving the point of plotting to half the resolution
Murat.geometry.x            = x+stepgx/2;
Murat.geometry.y            = y+stepgy/2;

index=0;
for i=1:nxc
    for j=1:nyc
        index               = index+1;
        XY(index,1:2)       = [x(i) y(j)];
    end
end
Murat.geometry.map          = XY;
lxy                         = length(XY(:,1));

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
DFolder                     = Murat.paths.datadir;
list                        = dir(DFolder);
filenames                   = {list.name}'; %create cell array of file names
filedir                     = {list.folder}'; %create cell array of file folder
listasac                    = filenames;
lls                         = length(listasac);

%variable containing coordinates of events and stations
evestaz                     = zeros(lls,6);
indexray                    = 0;

%creating he list
for i = 1:lls
    listasac{i,1}           = cat(2,filedir{i},'/',filenames{i});
    li                      = filenames{i,1};
    li1                     = cat(2,li(1:12),li(18:20));
        
    %Here it takes the event/station name. It is necessary to adapt the
    %numbers to where even name and station name are.
    if Murat.geometry.import == 1
        for ii = 1:length(nameeven)
            namee1          = nameeven{ii};
            namee           = namee1(1:12);
            
            % Loop over stations in staz.txt file
            for ir = 1:length(namestation)
                namest1     = namestation{ir};
                namest      = namest1(1:3);
                levst       = cat(2,namee,namest);
                
                if find(strncmp(li1,levst,length(levst)))>0
                    indexray= indexray+1;
                    sst     = [even(ii,:) stat(ir,:)];
                    evestaz(indexray,1:6)=sst;% Creates file of sel. coords
                    break
                end
            end
        end
    end
end

%Saving liststo Murat
Murat.paths.listasac        = listasac;
Murat.geometry.evestaz      = evestaz;

if exist(cat(2,Murat.paths.workingdir,Murat.paths.label),'dir')~=7
    mkdir(cat(2,Murat.paths.workingdir,Murat.paths.label))
end

%STILL UNTESTED
if Murat.geometry.availableVelocity > 0
    
    if Murat.geometry.availableVelocity ==  1
        
        modv1                   = load(Murat.data.namev); %1D velocity  model iasp91
        dend                    = -34000; %select maximum depth range
        
        li                      = length(modv1(:,1));
        modv                    = zeros(lxy*li,5);
        index                   = 0;
        for i = 1:nxc
            for j = 1:nyc
                index1          = index;
                index           = index+1;
                modv(index1*li+1:index1*li+li,1) = XY(index,1)*1000;
                modv(index1*li+1:index1*li+li,2) = XY(index,2)*1000;
                modv(index1*li+1:index1*li+li,3) =...
                    -modv1(1:li,1)*1000;
                modv(index1*li+1:index1*li+li,4) =...
                    modv1(1:li,PorS+1);
            end
        end
        
        modv(:,1)               = (modv(:,1)-modv(1,1));
        modv(:,2)               = (modv(:,2)-modv(1,2));
        modv(modv(:,3)<dend,:)  = [];
        Murat.geometry.modv     = modv;
        
    elseif Murat.geometry.availableVelocity==2
        
        % This works in [x,y,z], created for UTM WGS84 and origins to zero.
        modv                    = load(Murat.data.namev); %3D velocity model from text file
        modv(:,5)               = 0;
        modv(:,1)               = modv(:,1)+origin(1);
        modv(:,2)               = modv(:,2)+origin(2);
        
        
        Murat.geometry.modv     = modv;
    end
    
    chx=find(modv(:,1)     ~= modv(1,1),1);
    chy=find(modv(:,2)     ~= modv(1,2),1);
    resol2x                 = abs(modv(chx,1)-modv(1,1))/2;
    resol2y                 = abs(modv(chy,2)-modv(1,2))/2;
    resol2z                 = abs(modv(2,3)-modv(1,3))/2;
    
    %Steps of the velocity models
    passox                  = max(modv(:,1))-min(modv(:,1));
    passoy                  = max(modv(:,2))-min(modv(:,2));
    passoz                  = max(modv(:,3))-min(modv(:,3));
    
    passo                   = [passox passoy passoz];
    resol                   = [resol2x resol2y resol2z];
    resol2                  = min(resol);
    Murat.geometry.resolutionMin = resol2;
    
    %Creates grid for direct waves and check for zeroes
    %Regular step of the gridD for interpolation - half of step of modv
    
    ixD                     = floor(passox/resol2x);%numer of x,given the step of gridD
    iyD                     = floor(passoy/resol2y);%numer of y,given the step of gridD
    izD                     = floor(passoz/resol2z);%numer of depths,given the step of gridD
    
    gridD                   = zeros(3,max(passo./resol));
    gridD(1,1:ixD)          = origin(1):resol2x:origin(1)+passox-resol2x;
    gridD(2,1:iyD)          = origin(2):resol2y:origin(2)+passoy-resol2y;
    gridD(3,1:izD)          = -modv(1,3):resol2z:-modv(1,3)+passoz-resol2z;
    Murat.geometry.gridD    = gridD;
    % gridD dimensions
    pvel                    = zeros(ixD,iyD,izD);
    
    %NUMBER OF X, Y, AND Z LAYERS
    for k = 1:izD
        index               = 0;
        for i = 1:ixD
            for j = 1:iyD
                index       = index+1;
                pvel(i,j,k) = modv(index,4);
            end
        end
    end
    Murat.geometry.pvel     = pvel;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%