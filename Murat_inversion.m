%%  2D peak-delay and Qc TOMOGRAPHIC INVERSIONS
% Using Tikhonov inversions and L curves

function Murat=Murat_inversion(Murat)
%PATHS and FIGURES
FPath               =   Murat.paths.workingdir;
FLabel              =   Murat.paths.label;
fformat             =   Murat.figures.format;

%DATA
Qm                  =   Murat.data.measuredQc;
retainQm            =   Murat.data.retainQm;
RZZ                 =   Murat.data.uncertainty;
lpdelta             =   Murat.data.logPeakDelay;
outlierspd          =   Murat.data.outliersPeakDelay;

%GEOMETRY
XY                  =   Murat.geometry.map;
nxc                 =   Murat.geometry.gridX;
nyc                 =   Murat.geometry.gridY;

% INVERSION
sizea               =   Murat.inversion.sizeCheck;
latt                =   Murat.inversion.lowCheck;
hatt                =   Murat.inversion.highCheck;
Apd                 =   Murat.inversion.APeakDelay;
Ac                  =   Murat.inversion.AQCoda;
nonlinear           =   Murat.inversion.nonlinear;

pd                  =   XY(:,1:2);
pd(:,3)             =   -1000;
Qc                  =   pd;

%Peak delay mapping - weighted average
lpdelta_o           =   lpdelta(outlierspd==0); %only without outliers
Apd                 =   Apd(outlierspd==0,:); %only without outliers
lApd                =   size(Apd);
Apd1                =   Apd;
mpd                 =   zeros(lApd(2),1); %final peak delay map

% For loop to multiply each Q measured by the weight (segment length)
for j = 1:lApd(2)
        Apd1(:,j)   =   Apd(:,j).*lpdelta_o;
end

% For loop to sum and normalize
for j = 1:lApd(2)
    if sum(Apd1(:,j)) ~= 0
        mpd(j,1)    =   sum(Apd1(:,j))/sum(Apd(:,j));
    end
end

%Assign to 4th column of what you save
pd(:,4)             =   mpd;

%Qc mapping
%Remove anomalous Qm
stat                =   Qm(retainQm);
Ac1                 =   Ac(retainQm,:);
RZZ1                =   RZZ(retainQm,1);

% Weighting depends on RZZ
if nonlinear == 0
    
    W1              =   RZZ1<0.3;
    W2              =   RZZ1<0.5;
    W3              =   RZZ1<0.7;
    W4              =   W3+W2+W1;
    Wc              =   diag(1./(W4+1));% weights
% For the nonlinear case there is no weighting
elseif nonlinear == 1
    
    Wc              =   1;
    
end

%Weighted tikhonov
dcW                 =   Wc*stat;
Gc                  =   Wc*Ac1;
[Uc,Sc,Vc]          =   svd(Gc);

% smooting parameter is user defined
LcQc                =   figure('Name','L-curve Qc','NumberTitle','off');
l_curve(Uc,diag(Sc),dcW,'Tikh')
% define from L-curve
tik0_regC           =   input('Your personal smoothing parameter for coda ');
FName               =   'Lc_Qc';
saveas(LcQc,fullfile(FPath, FLabel, FName), fformat);
close(LcQc)

% picard plot
PpQc                =   figure('Name','Picard-plot Qc',...
    'NumberTitle','off','visible','off');
picard(Uc,diag(Sc),dcW);
FName               =   'Picard_Qc';
saveas(PpQc,fullfile(FPath, FLabel, FName), fformat);

% invert
mtik0C              =   tikhonov(Uc,diag(Sc),Vc,dcW,tik0_regC);
Qc(:,4)             =   mtik0C;

%Testing - Creating 2D checkerboard matrix
nxc1                =   nxc/sizea;
nyc1                =   nyc/sizea;
I                   =   imresize(xor(mod(1:nyc1, 2).', mod(1:nxc1, 2)),...
    [sizea*nyc1 sizea*nxc1], 'nearest');
Qc(:,5)             =   I(1:end);
Qc(Qc(:,5)==1,5)    =   latt*100;
Qc(Qc(:,5)==0,5)    =   hatt*100;

% invert checkerboard
Qc5                 =   Qc(:,5);
re                  =   Gc*Qc5;
mcheckc             =   tikhonov(Uc,diag(Sc),Vc,re,tik0_regC);
Qc(:,6)             =   mcheckc;
Qc(:,7)             =   Gc(1,:);

% save in Murat
Murat.inversion.Qc          =   Qc;
Murat.inversion.peakDelay   =   pd;

% save peak-delay
FName                       =   'peakdelay.txt';
save(fullfile(FPath, FLabel, FName), 'pd','-ascii');
% save Qc
FName                       =   'Qc.txt';
save(fullfile(FPath, FLabel, FName), 'Qc','-ascii');