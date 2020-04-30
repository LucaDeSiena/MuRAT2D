%% Seismic attributes for peak delay and Qc imaging
function Murat=Murat_data(Murat)
%PATHS
pa              =   Murat.analysis;
f0              =   Murat.paths.originTime;
fP              =   Murat.paths.PTime;
fS              =   Murat.paths.STime;

%DATA
listasac        =   Murat.paths.listasac;
lls             =   length(listasac);
tCm             =   Murat.data.startLT;
tWm             =   Murat.data.codaWindow;
cf              =   Murat.data.centralFrequency;
mintpde         =   Murat.data.minimumPD;
maxtpde         =   Murat.data.maximumPD;
nf              =   Murat.data.smoothing;
sped            =   Murat.data.spectralDecay;
vP              =   Murat.data.averageVelocityP;
vS              =   Murat.data.averageVelocityS;
PorS            =   Murat.data.PorS;
compon          =   Murat.data.components;

%GEOMETRY
evst            =   Murat.geometry.import;
XY              =   Murat.geometry.map;
degorutm        =   Murat.geometry.degreesorutm;
stepgx          =   Murat.geometry.gridStepX;
stepgy          =   Murat.geometry.gridStepY;

if pa==3
    kT          =   Murat.geometry.kernelTreshold;
end

xx              =   XY(:,1);
yy              =   XY(:,2);
lxy             =   length(xx);
x1              =   Murat.geometry.x;

% INVERSION
nonlinear       =   Murat.inversion.nonlinear;
if nonlinear == 1
    nW          =   Murat.inversion.fitL;
    ntW         =   Murat.inversion.fitN;
    m1a         =   Murat.inversion.fit;
    L1          =   Murat.inversion.total;
end
RZZ2            =   Murat.inversion.fitT;

if evst == 1
    evestaz     =   Murat.geometry.evestaz;
elseif evst == 2
    evestaz     =   zeros(lls,6);
end

%Set up variables to save
Qm              =   zeros(lls,1); %Qc with linearised or non-linear theory
RZZ             =   zeros(lls,1); %Correlation coefficient with respect to linear
peakd           =   zeros(lls,1); %Peak delay
ttheory         =   zeros(lls,1); %Theorethical time in case of missing origin
D               =   zeros(lls,1); %Hypocentral distance
tPS             =   zeros(lls,1); %travel time of chosen phase

%Predefine the inversion matrix for peak-delay and coda-Q imaging
Ac              =   zeros(lls,lxy);
Apd             =   zeros(lls,lxy);

indexlonger     =   0; %number of files longer than chosen coda window

%=========================================================================
%% Loop through SAC files
% A loop to measure Qc and peak-delay for each seismic trace located in
% a folder - it only accepts SAC files.
% The seismograms may be the output of both single and three-component
% seismic stations. In the case of more than one component the rays must be
% ordered in the folder starting with the WE component, then SN and finally
% the Z component for each ray.
%
% The program accepts SAC files and uses the information in the header.
% The header must include peaking of the direct phase of interest
% (preferred is marker "a" for P and "t0" for S).
% If evst=2 you can get event and station info
% directly from the header - this is the preferred option.

for i = 1:lls %loop through source-station pairs
    
    %Display every 200 waveforms
    if isequal(mod(i,200),0)
        disp(['Waveform number is ', num2str(i)])
    end
    
    %OPERATIONS ON WAVEFORM
    %Read seismogram and get peaking/event/station info
    [tempis,sisma,SAChdr]   =   fget_sac(listasac{i}); % Imports SAC files
    srate                   =   1/SAChdr.times.delta; %sampling frequency
    sisma                   =   detrend(sisma,1);%remove trend
    lsis                    =   length(sisma);
    tu                      =   tukeywin(lsis,0.05);%avoid windowing
    tsisma                  =   tu.*sisma;
    
    % Filter creation - in loop as sampling might change
    Wn                      =   ([cf-cf/3 cf+cf/3]/srate*2); %frequency band
    [z,p,k]                 =   butter(4,Wn,'bandpass'); %butter filter
    [sos,g]                 =   zp2sos(z,p,k); % Convert to SOS form
    %Hd                      =   dfilt.df2tsos(sos,g); % Create a dfilt object
    
    %Now using filtfilt to avoid phase shift
    fsisma                  =   filtfilt(sos,g,tsisma);% filtered waveform
    hsp                     =   hilbert(fsisma); %hilbert of the waveform
    sp                      =   smooth(abs(hsp),nf/cf*srate);%ms of the waveform
    
    t00                     =   tempis(1); %starting time of the waveform
    
    %Check if you are doing P or S peak delay
    if PorS == 2 %P-wave peak delay
        pktime              =   eval(fP);%should find a way to improve here
    elseif PorS == 3 %S-wave peak delay
        pktime              =   eval(fS);%should find a way to improve here
    end
    
    %Skip waveform if picking is not on the waveform
    if pktime < tempis(1) || pktime > tempis(end)
        continue
    end
    
    %CREATE KERNELS
    if evst == 2 %creating source-station coordinates if from SAC files
        even                =   [SAChdr.event.evlo SAChdr.event.evla...
            SAChdr.event.evdp]; %event location
        
        stati               =   [SAChdr.station.stlo SAChdr.station.stla...
            SAChdr.station.stel]; %station location
        
        %Conditions to avoid missing coordinates and calculate segments
        if stati(end) == -12345 || even(end) == -12345
            stati(end)      =   0; even(end)=0;
        end
        if even(end) > 0
            even(end)       =   -even(end);
        end
        if stati(end) < 0
            stati(end)      =   -stati(end);
        end
        if stati(end) > 2
            stati(end)      =   stati(end)/1000;
        end
        evestaz(i,:)        =   [even stati];
    end
    
    %Take coordinates and put it in variables to save
    sst                     =   evestaz(i,:);

    % Hypocentral distance
    if degorutm == 1 %UTM case
        metordeg            =   1000;
        Distance            =   sqrt((sst(4)-sst(1))^2+(sst(5)-sst(2))^2+...
            (sst(6)-sst(3))^2)/metordeg;
        D(i,1)              =   Distance/degorutm;
        
    elseif degorutm == 111 %Lat/long case - to be improved
        Distance            =   sqrt((sst(4)-sst(1))^2*degorutm^2+...
            (sst(5)-sst(2))^2*degorutm^2+(sst(6)-sst(3))^2);
        D(i,1)              =   Distance;
    end
    
    %Conditions in case the zero time is missing in the header
    if isequal(f0,[]) || isequal(eval(f0),-12345)
        ttheory(i,1)        =   Distance/vP;
        t0                  =   pktime-ttheory(i,1);
    else
        t0                  =   eval(f0);
    end
    
    % Picking time
    tPS(i,1)                =   pktime-t0;
    
    %starting-sample of the direct-wave window
    cursor1                 =   floor((pktime-t00)*srate);
    cursor2                 =   floor(cursor1+maxtpde*srate);
    
    %Measurse for peak delay - compute shorter window if waveform is cut
    if cursor2 > lsis
        pdcursor            =   lsis; %set to length of seismogram
    else
        pdcursor            =   cursor2; %end sample for peak-delay
    end
    
    % Envelope for the peak-delay measurement
    pdsp=sp(cursor1:pdcursor);
    [~,tspm]=max(pdsp);
    peakd(i,1)=tspm/srate;
    
    %Compute Qc from late lapse times - in case you miss it you are
    %measuring envelopes from the peak of the direct wave, assuming that
    %the entire waveform is diffusive (e.g., Deception Island).
    if isequal(tCm,[])
        tC                  =   (pktime-t00-t0)+tspm/srate; %in case you start from max
    else
        tC                  =   tCm;
    end
    
    %starting-sample of the coda-wave window
    cursorc3                =   floor((t0-t00+tC-1)*srate);%coda start sample
    cursorc4                =   floor(cursorc3 + tWm*srate-1);%coda end sample
    
    %Measurse for Qc - compute shorter window if waveform is cut
    if cursorc4 > lsis
        indexlonger         =   indexlonger+1;
        spcm                =   sp(cursorc3:lsis);
    else
        spcm                =   sp(cursorc3:cursorc4);
    end
    lspm                    =   length(spcm);
    
    %Compute time vector
    tm                      =   (tC+1/srate:1/srate:tC+lspm/srate)';
    
    %Only evaluate central time series
    edgeno                  =   floor(0.05*length(tm));
    tm1                     =   tm(edgeno:end-edgeno);
    spcm1                   =   spcm(edgeno:end-edgeno);
    
    %Using the linear approximation
    if nonlinear == 0
        
        %THIS IS THE DATA VECTOR WITH LINEARISED THEORY
        
        EWz                 =   spcm1.*tm1.^sped; %Calculating energy
        lspmz               =   log(EWz)/2/pi/cf; % source-station data
        Rz                  =   corrcoef([tm1,lspmz]); %sets uncertainty
        polyz               =   polyfit(tm1,lspmz,1); %Qc^-1
        
        if polyz(1)<0%only where you have a positive Qc
            Qm(i,1)         =   -polyz(1);
            RZZ(i,1)        =   abs(Rz(1,2));
        end
    
    %Solving with the grid-search
    elseif nonlinear == 1
        
        %THIS IS THE SYSTEM OF EQUATIONS FOR THE NON-LINEAR SOLUTION
        %USING THE LINEAR INVERSE QC AS STARTING MODEL
        
        tlapse              =   (tC+nW/2:nW:tC+tWm-nW/2)'; %lapse time
        d_obs               =   zeros(ntW,1);%set up the data vector
        
        %For all the windows in play...
        for k = 1:ntW
            ntm             =   (k-1)*nW*srate + 1:k*nW*srate; %take the samples in the window...
            d_obs(k,1)      =   mean(spcm(floor(ntm))); %and compute the average energy.
        end
        d_obs1              =   d_obs(1:end-1)/d_obs(end);%normalize by the last window
        
        %PREDICTED BY THEORY
        E                   =   zeros(L1,1);
        for n = 1:L1
            d_pre           =   tlapse.^(-sped).*exp(-2*pi*cf.*tlapse*m1a(n));%data predicted
            d_pre1          =   d_pre(1:end-1)/d_pre(end);%normalized
            E(n,1)          =   sum(abs(d_obs1-d_pre1));%L1 residual
        end
        [Emin, indexE]      =   min(E);%find the minimum to set the right Qc
        Qm(i,1)             =   m1a(indexE);%Qc set
        % Uncertainty from the minimum residual, the bigger it is the more
        % uncertain the model becomes
        RZZ(i,1)            =   1/Emin;
    end
    
    %Geometry required to set the length of the segments used to weight
    %peak-delay mapping
    miix                    =   min(sst(1),sst(4));
    maax                    =   max(sst(1),sst(4));
    miiy                    =   min(sst(2),sst(5));
    maay                    =   max(sst(2),sst(5));
    if abs(maax-miix) > abs(x1(2)-x1(1)) && abs(maay-miiy) > abs(yy(2)-yy(1))
        fxy                 =   find(xx>miix & xx<maax & yy>miiy & yy<maay);
    elseif abs(maax-miix) > abs(x1(2)-x1(1)) && abs(maay-miiy) < abs(yy(2)-yy(1))
        fxy                 =   find(xx>miix & xx<maax & yy>miiy & maay<yy);
    elseif abs(maax-miix) < abs(x1(2)-x1(1)) && abs(maay-miiy) > abs(yy(2)-yy(1))
        fxy                 =   find(xx>miix & maax<xx & yy>miiy & yy<maay);
    else
        fxy                 =   find(xx>miix & maax<xx & yy>miiy & maay<yy);
    end
        
    lfxy                    =   length(fxy);
    
    xd                      =   linspace(miix,maax,lfxy);
    yd                      =   linspace(miiy,maay,lfxy);
    dist                    =   sqrt((maax-miix)^2+(maay-miiy)^2)/(lfxy+1);
    
    for l = 1:length(xd)
        pointer             =   [xd(l), yd(l)];
        targ                =   XY(:,1:2);
        %compute Euclidean distances:
        distances           =   sqrt(sum(bsxfun(@minus, targ, pointer).^2,2));
        %find the smallest distance and use that as an index into B:
        mind                =   distances==min(distances);
        Apd(i,mind)         =   Apd(i,mind)+dist;
        Ac(i,mind)          =   Ac(i,mind)+dist;      
    end
    
    
    %Inversion matrix for Qc
    if pa == 2 % Approximate kernels
        
        F                   =   kernel_analytic(sst,XY); %see function
        %Inversion matrix for Qc
        Ac(i,:)             =   F;
        
    elseif pa == 3 %Pacheco-Snieder kernels
        
        [K_grid,r_grid]     =...
            kernels_diffusive(tC+tWm/2,sst(1:3),sst(4:6),XY,degorutm,vS,kT); %see function
        K1                  =   K_grid;
        
        %removing the eventual poles - see function
        K1(isnan(K1))           =   [];
        K_grid(isnan(K_grid))   =  max(K1);
        
        %Assigning kernels to the grid
        rK                      =   [r_grid K_grid];
        rK1                     =   XY;
        for ixy = 1:lxy
            x                   =   XY(ixy,1);
            y                   =   XY(ixy,2);
            iixy                =   find(rK(:,1)>x-stepgx/2 &...
                rK(:,1)<x+stepgx/2 & rK(:,2)>y-stepgy/2 ...
                & rK(:,2)<y+stepgy/2);
            if iixy
                rK1(ixy,3)      =   mean(rK(iixy,4));
            end
        end
        
        %Normalization
        if ~isequal(sum(rK1(:,3)),0)
            Ac(i,1:lxy)         =   rK1(:,3)/sum(rK1(:,3));
        end
        
    end
end

%% Setting up the final data vectors with checks on values

%Find the amount of Qc unavailable in the data
noQm                            =   sum(Qm==0)/length(Qm)*100;

%Conditions for the linear and non-linear inversions are opposite
if nonlinear == 0
    noRZZ                       =   sum(RZZ<=RZZ2)/length(RZZ)*100;
elseif nonlinear == 1
    noRZZ                       =   sum(RZZ>=RZZ2)/length(RZZ)*100;
end

if noQm > 20
    Qless                       =   ...
        ['Attention: ',num2str(noQm),'% of your Q are =0'];
    warning(Qless)
end

if noRZZ>20
    Rless                       =   ['Attention: ',num2str(noRZZ),...
        '% of your coerrelation coefficients are lower than treshold'];
    warning(Rless)
end

%% Setting up the data vector in case of 2- and 3-components data
icomp                           =   0;
ll                              =   lls;
lA                              =   length(Ac(:,1));

%Only one waveform component
if compon > 1 
    Ac                          =   Ac(1:compon:lA,:);
    Apd                         =   Apd(1:compon:lA,:);
    Murat.inversion.APeakDelay  =   Apd;
    Murat.inversion.AQCoda      =   Ac;
end
lu                              =   D(1:compon:lA);
time0                           =   tPS(1:compon:lA);

% 2 components (typically [WE, SN])
if compon ==  2
    lsig                        =   ll/2;
    Qm1                         =   zeros(lsig,1);
    RZZ1                        =   zeros(lsig,1); 
    peakd1                      =   zeros(lsig,1);
    
    for i = 1:compon:(ll-1)
        icomp                   =   icomp+1;
        Qm(icomp,1)            =   (Qm(i,1)+Qm(i+1))/2;
        RZZ(icomp,1)           =   (RZZ(i,1)+RZZ(i+1))/2;
        peakd(icomp,1)         =   (peakd(i,1)+peakd(i+1))/2;
    end
    
    Qm                          =   Qm1;
    RZZ                         =   RZZ1;
    peakd                       =   peakd1;

% 3 components (WE, SN, Z)
elseif compon == 3
    lsig                        =   ll/3;
    Qm1                         =   zeros(lsig,1); %
    RZZ1                        =   zeros(lsig,1); %
    peakd1                      =   zeros(lsig,1); %
    
    for i = 1:3:(ll-2)
        icomp                   =   icomp+1;
        
        %Averaging the horizontals before the vertical
        Qm(icomp,1)            =   ((Qm1(i)+Qm1(i+1))/2 + Qm1(i+2))/2;
        RZZ(icomp,1)           =   ((RZZ1(i)+RZZ1(i+1))/2 + RZZ1(i+2))/2;
        peakd(icomp,1)         =   ((peakd(i)+peakd(i+1))/2 +peakd(i+2))/2;
    end
    
    Qm                          =   Qm1;
    RZZ                         =   RZZ1;
    peakd                       =   peakd1;
end

%% Peak-delay, for reference see e.g. Calvet et al. 2013, Tectonophysics
timepd                          =   time0;
% setting up the trehold for acceptance to less than the period
nop                             =   find(peakd<1/cf);

if nop
    peakd(nop)                  =   [];
    lu(nop)                     =   [];
    Apd(nop,:)                  =   [];
    timepd(nop,:)               =   [];
end

%De Siena et al. 2016 / using Vp/Vs to map max of S waves in case of P
vpvs                            =   vP/vS;
l10pd                           =  log10(peakd);
%Case for P waves
if PorS == 2
    l10l                        =   log10(timepd*vpvs);
    time0                       =   time0*vpvs;
    
% Case for S waves
elseif PorS == 3
    l10l                        =   log10(timepd);
end
fitrobust                       =   fit(l10l,l10pd,'poly1','Robust','on');
pab                             =   [fitrobust.p1 fitrobust.p2];
l10pdt                          =   polyval(pab,l10l);
lpdelta                         =   l10pd-l10pdt;

% To remove outliers
I                               = abs(lpdelta) > 2*std(lpdelta) | ...
    peakd >= maxtpde | peakd < mintpde;
outlierspd                      =   excludedata(l10l,lpdelta,'indices',I);

%Same for Qc
%Remove anomalous Q remembering that there are two different types of
%inversion
if nonlinear == 0
    retainQm                    =   find(Qm>0 & RZZ>RZZ2);
    mQm                         =   mean(Qm(retainQm));
    outliersc                   =   Qm > mQm+2*std(Qm(retainQm));
    retainQm                    =   find(Qm>0 & RZZ>RZZ2 & outliersc==0);
    discardQm=find(Qm<=0 | RZZ<RZZ2 | outliersc==1);
    
elseif nonlinear == 1
    retainQm                    =   find(Qm>0 & RZZ<RZZ2);
    mQm                         =   mean(Qm(retainQm));
    outliersc                   =   Qm > mQm+2*std(Qm(retainQm));
    retainQm                    =   find(Qm>0 & RZZ<RZZ2 & outliersc==0);
    discardQm                   =   find(Qm<=0 | RZZ>RZZ2 | outliersc==1);
end

%save quantities for the data vector
Murat.data.hypocentralDistance  =   lu;
Murat.data.measuredQc           =   Qm;
Murat.data.uncertainty          =   RZZ;
Murat.data.logTravelPD          =   l10l;
Murat.data.logPeakDelay         =   lpdelta;
Murat.data.outliersPeakDelay    =   outlierspd;
Murat.data.theoreticalTravelTime=   time0;
Murat.data.averageQc            =   mQm;
Murat.data.retainQm             =   retainQm;
Murat.data.discardQm            =   discardQm;
Murat.data.fitrobust            =   fitrobust;
Murat.data.peakd                =   peakd;
Murat.inversion.APeakDelay      =   Apd;
Murat.inversion.AQCoda          =   Ac;
Murat.geometry.evestaz          =   evestaz;