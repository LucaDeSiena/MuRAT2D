%%  2D peak-delay and Qc TOMOGRAPHIC INVERSIONS
function Murat=Murat_inversion(Murat)
%PATHS and FIGURES
FPath=Murat.paths.workingdir;
FLabel=Murat.paths.label;
fformat=Murat.figures.format;

%DATA
Qm=Murat.data.measuredQc;
retainQm=Murat.data.retainQm;
RZZ=Murat.data.uncertainty;
lpdelta=Murat.data.logPeakDelay;
outlierspd=Murat.data.outliersPeakDelay;

%GEOMETRY
XY=Murat.geometry.map;
nxc=Murat.geometry.gridX;
nyc=Murat.geometry.gridY;

% INVERSION
sizea=Murat.inversion.sizeCheck;
latt=Murat.inversion.lowCheck;
hatt=Murat.inversion.highCheck;
Apd=Murat.inversion.APeakDelay;
Ac=Murat.inversion.AQCoda;
nonlinear=Murat.inversion.nonlinear;

pd=XY(:,1:2);
pd(:,3)=-1000;
Qc=pd;

%Peak delay mapping
lpdelta_o=lpdelta(outlierspd==0);
Apd=Apd(outlierspd==0,:);
%Apd(Apd~=0)=1;
lApd=size(Apd);
Apd1=Apd;
for i=1:lApd(1)
    Apd1(i,:)=Apd1(i,:)/sum(Apd(i,:))*lpdelta_o(i);
end


dp=lpdelta_o;
[Up,Sp,Vp]=svd(Apd);

% smooting parameter is user defined
Lcpd=figure('Name','L-curve Peak Delay','NumberTitle','off');
l_curve(Up,diag(Sp),dp,'Tikh')
tik0_regPD=input('Your personal smoothing parameter for peak delay ');
FName = 'Lc_PeakDelay';
saveas(Lcpd,fullfile(FPath, FLabel, FName), fformat);
close(Lcpd)

% picard plot
Pppd=figure('Name','Picard-plot Peak Delay',...
    'NumberTitle','off','visible','off');
picard(Up,diag(Sp),dp);
FName = 'Picard_PeakDelay';
saveas(Pppd,fullfile(FPath, FLabel, FName), fformat);

mtik0PD=tikhonov(Up,diag(Sp),Vp,dp,tik0_regPD);
pd(:,4)=mtik0PD;

%Testing - Creating 2D checkerboard matrix
nxc1=nxc/sizea;
nyc1=nyc/sizea;
I = imresize(xor(mod(1:nyc1, 2).', mod(1:nxc1, 2)),...
    [sizea*nyc1 sizea*nxc1], 'nearest');
pd(:,5)=I(1:end);
pd(pd(:,5)==1,5)=latt*100;
pd(pd(:,5)==0,5)=hatt*100;

pd5= pd(:,5);
rePD = Apd*pd5;
mcheckpd=tikhonov(Up,diag(Sp),Vp,rePD,tik0_regPD);
pd(:,6)=mcheckpd;
pd(:,7)=Apd(1,:);

%Qc mapping
%Remove anomalous Qm
stat=Qm(retainQm);
Ac1=Ac(retainQm,:);
RZZ1=RZZ(retainQm,1);
if nonlinear==0
    
    W1=RZZ1<0.3;
    W2=RZZ1<0.5;
    W3=RZZ1<0.7;
    W4=W3+W2+W1;
    Wc=diag(1./(W4+1));% weights
    
elseif nonlinear==1
    
    Wc=1;
    
end

dcW=Wc*stat;
Gc=Wc*Ac1;
[Uc,Sc,Vc]=svd(Gc);

% smooting parameter is user defined
LcQc=figure('Name','L-curve Qc','NumberTitle','off');
l_curve(Uc,diag(Sc),dcW,'Tikh')
tik0_regC=input('Your personal smoothing parameter for coda ');
FName = 'Lc_Qc';
saveas(LcQc,fullfile(FPath, FLabel, FName), fformat);
close(LcQc)

% picard plot
PpQc=figure('Name','Picard-plot Qc',...
    'NumberTitle','off','visible','off');
picard(Uc,diag(Sc),dcW);
FName = 'Picard_Qc';
saveas(PpQc,fullfile(FPath, FLabel, FName), fformat);

mtik0C=tikhonov(Uc,diag(Sc),Vc,dcW,tik0_regC);
Qc(:,4)=mtik0C;

%Testing - Creating 2D checkerboard matrix
Qc(:,5)=pd(:,5);

Qc5= Qc(:,5);
re = Gc*Qc5;
mcheckc=tikhonov(Uc,diag(Sc),Vc,re,tik0_regC);
Qc(:,6)=mcheckc;
Qc(:,7)=Gc(1,:);

Murat.inversion.Qc=Qc;
Murat.inversion.peakDelay=pd;

% save peak-delay
FName = 'peakdelay.txt';
save(fullfile(FPath, FLabel, FName), 'pd','-ascii');
% save Qc
FName = 'Qc.txt';
save(fullfile(FPath, FLabel, FName), 'Qc','-ascii');