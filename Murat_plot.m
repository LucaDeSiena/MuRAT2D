%%  Creates maps
function Murat=Murat_plot(Murat)
hasMT=Murat.figures.hasMT;

pa=Murat.analysis;

%PATHS and FIGURES
FPath=Murat.paths.workingdir;
FLabel=Murat.paths.label;
visib=Murat.figures.visibility;
sz=Murat.figures.sizeMarker;
fformat=Murat.figures.format;

%DATA
lls=length(Murat.paths.listasac);
Qm=Murat.data.measuredQc;
outlierspd=Murat.data.outliersPeakDelay;
time0=Murat.data.theoreticalTravelTime;
mQm=Murat.data.averageQc;
retainQm=Murat.data.retainQm;
l10l=Murat.data.logTravelPD;
fitrobust=Murat.data.fitrobust;
peakd=Murat.data.peakd;

%GEOMETRY
origin=Murat.geometry.origin;
stepgx=Murat.geometry.gridStepX;
stepgy=Murat.geometry.gridStepY;
nxc=Murat.geometry.gridX;
nyc=Murat.geometry.gridY;
x=Murat.geometry.x;
y=Murat.geometry.y;
degorutm=Murat.geometry.degreesorutm;
evestaz=Murat.geometry.evestaz;
Ac=Murat.inversion.AQCoda;
pd=Murat.inversion.peakDelay;
Qc=Murat.inversion.Qc;

pdel=zeros(length(x),length(y));
QQc=zeros(length(x),length(y));
pdchi=zeros(length(x),length(y));
pdcho=zeros(length(x),length(y));
QQchi=zeros(length(x),length(y));
QQcho=zeros(length(x),length(y));
[X,Y]=meshgrid(x,y);
index=0;
for i=1:length(x)
    for j=1:length(y)
        
        index=index+1;
        pdel(i,j)=pd(index,4);
        QQc(i,j)=Qc(index,4);
        pdchi(i,j)=pd(index,5);
        pdcho(i,j)=pd(index,6);
        QQchi(i,j)=Qc(index,5);
        QQcho(i,j)=Qc(index,6);
        
    end
end

if degorutm==111 && hasMT
    load coastlines coastlat coastlon
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS

rays=figure('Name','Rays','NumberTitle','off','visible',visib,...
    'Position',[300,200,800,600]);
if degorutm==111 && hasMT==1
    hold on
    geoshow(coastlat,coastlon);
    xlim([origin(1) origin(1)+nxc*stepgx]);
    ylim([origin(2) origin(2)+nyc*stepgy]);
end

for nn=1:lls
    hold on
    plot([evestaz(nn,1) evestaz(nn,4)],...
        [evestaz(nn,2) evestaz(nn,5)],'k-')
end
hold on
scatter(evestaz(:,1),evestaz(:,2),sz,'c','MarkerEdgeColor',...
    [1 1 1], 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)
hold on
scatter(evestaz(:,4),evestaz(:,5),sz,'^','MarkerEdgeColor',...
    [1 1 1], 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)
hold off
grid on
ax = gca;
ax.GridLineStyle = '-';
ax.GridColor = 'k';
ax.GridAlpha = 1;
ax.LineWidth = 1;
FName = 'Rays';
saveas(rays,fullfile(FPath, FLabel, FName), fformat);

if pa>1
    
    Qcsen=figure('Name','Qc sensitivity, first source-station pair',...
        'NumberTitle','off','visible',visib,'Position',[300,200,800,600]);
    Qcss=Ac(1,:);
    Qcs=zeros(nxc,nyc);
    
    index=0;
    for i=1:length(x)
        for j=1:length(y)
            index=index+1;
            Qcs(i,j)=Qcss(index);
        end
    end
    
    [X,Y]=meshgrid(x,y);
    contourf(X,Y,Qcs')
    
    colorbar
    grid on
    ax = gca;
    ax.GridLineStyle = '-';
    ax.GridColor = 'k';
    ax.GridAlpha = 1;
    ax.LineWidth = 1;
    if degorutm==111 && hasMT==1
        hold on
        geoshow(coastlat,coastlon);
        xlim([origin(1) origin(1)+nxc*stepgx]);
        ylim([origin(2) origin(2)+nyc*stepgy]);
    end
    
    FName = 'Qc_sensitivity';
    saveas(Qcsen,fullfile(FPath, FLabel, FName), fformat);
end

%plot to check that Qc is constant with travel time and peak delays
%increase with travel time

Qcpd=figure('Name','Qc and peak-delays','NumberTitle','off',...
    'visible',visib,'Position',[300,200,800,600]);
subplot(2,1,1)
plot(time0(retainQm),Qm(retainQm),'o',...
    'MarkerSize',6,'MarkerEdgeColor',[0 0 0])
hold on
plot(time0(retainQm),mQm*ones(length(time0(retainQm)),1),'r-',...
    'LineWidth',2)
title('Dependence of Qc^{-1} on travel time');
xlabel('Travel time (s)','FontSize',12,'FontWeight','bold','Color','k')
ylabel('Qc^{-1}','FontSize',12,...
    'FontWeight','bold','Color','k')
legend({'Qc^{-1}',cat(2,'<Qc> = ',num2str(1/mQm))},...
    'Location','northeast')

subplot(2,1,2)
plot(fitrobust,'k--',l10l,log10(peakd),'ko',outlierspd,'r*')
xti=xticks;
xtl=length(xti);
xt=cell(xtl,1);
for i=1:length(xti)
    xt(i,1)={10^xti(i)};
end
nameText=cat(2,'A(f) = ',num2str(fitrobust.p2),', B(f) = ',...
    num2str(fitrobust.p1));
xticklabels(xt)
title('Dependence of log. peak delays on log. travel time');
xlabel('Travel time (s)','FontSize',12,'FontWeight','bold',...
    'Color','k')
ylabel('Log. Peak delay','FontSize',12,...
    'FontWeight','bold','Color','k')
legend('Location','southeast')
text(mean(l10l),max(log10(peakd)),nameText,'HorizontalAlignment','center',...
    'FontSize',14)
FName = 'Qc_Peak_Delay';
saveas(Qcpd, fullfile(FPath, FLabel, FName), fformat);


pdmap=figure('Name','Peak-delay map','NumberTitle','off',...
    'visible',visib,'Position',[300,200,800,600]);
contourf(X,Y,pdel');
axis equal
view(2)
colormap(autumn);
hcb=colorbar;
title(hcb,'Log. Peak Delay','FontSize',14,'FontWeight','bold','Color','k');

xlabel('WE','FontSize',12,'FontWeight','bold','Color','k')
ylabel('SN','FontSize',12,'FontWeight','bold','Color','k')
title('Peak-delay variations',...
    'FontSize',12,'FontWeight','bold','Color','k');
hold on
scatter(evestaz(:,1),evestaz(:,2),sz,'c','MarkerEdgeColor',...
    [1 1 1], 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)
hold on
scatter(evestaz(:,4),evestaz(:,5),sz,'^','MarkerEdgeColor',...
    [1 1 1], 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)
if degorutm==111 && hasMT==1
    hold on
    geoshow(coastlat,coastlon);
    xlim([origin(1) origin(1)+nxc*stepgx]);
    ylim([origin(2) origin(2)+nyc*stepgy]);
end
hold off
FName = 'Peak_delay_map';
saveas(pdmap,fullfile(FPath, FLabel, FName), fformat);

Qcmap=figure('Name','Qc map','NumberTitle','off',...
    'visible',visib,'Position',[300,200,800,600]);
contourf(X,Y,QQc');
axis equal
view(2)
colormap(copper)

hcb=colorbar;
title(hcb,'Qc^{-1}','FontSize',14,'FontWeight','bold','Color','k');

xlabel('WE','FontSize',12,'FontWeight','bold','Color','k')
ylabel('SN','FontSize',12,'FontWeight','bold','Color','k')
title('Coda attenuation variations','FontSize',12,'FontWeight','bold',...
    'Color','k');
hold on
scatter(evestaz(:,1),evestaz(:,2),sz,'c','MarkerEdgeColor',...
    [1 1 1], 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)
hold on
scatter(evestaz(:,4),evestaz(:,5),sz,'^','MarkerEdgeColor',...
    [1 1 1], 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)
if degorutm==111 && hasMT==1
    hold on
    geoshow(coastlat,coastlon);
    xlim([origin(1) origin(1)+nxc*stepgx]);
    ylim([origin(2) origin(2)+nyc*stepgy]);
end
hold off
FName = 'Qc_map';
saveas(Qcmap,fullfile(FPath, FLabel, FName), fformat);

%Checkerboard Peak Delay
PDcheck=figure('Name','Peak Delay checkerboard test','NumberTitle',...
    'off','visible',visib,'Position',[300,200,1000,600]);
subplot(1,2,1)
contourf(X,Y,pdchi');
axis equal
view(2)
colormap(gray)

hcb=colorbar;
title(hcb,'log10 peak delay','FontSize',14,'FontWeight','bold','Color','k');

xlabel('WE','FontSize',12,'FontWeight','bold','Color','k')
ylabel('SN','FontSize',12,'FontWeight','bold','Color','k')
title('Qc checherboard input',...
    'FontSize',12,'FontWeight','bold','Color','k');

hold on
scatter(evestaz(:,1),evestaz(:,2),sz,'c','MarkerEdgeColor',...
    [1 1 1], 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)
hold on
scatter(evestaz(:,4),evestaz(:,5),sz,'^','MarkerEdgeColor',...
    [1 1 1], 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)
if degorutm==111 && hasMT==1
    hold on
    geoshow(coastlat,coastlon);
    xlim([origin(1) origin(1)+nxc*stepgx]);
    ylim([origin(2) origin(2)+nyc*stepgy]);
end
hold off

subplot(1,2,2)
contourf(X,Y,pdcho');
axis equal
view(2)
colormap(gray)

hcb=colorbar;
title(hcb,'log10 peak delay','FontSize',14,'FontWeight','bold','Color','k');

xlabel('WE','FontSize',12,'FontWeight','bold','Color','k')
ylabel('SN','FontSize',12,'FontWeight','bold','Color','k')
title('Peak Delay checherboard output',...
    'FontSize',12,'FontWeight','bold','Color','k');

hold on
scatter(evestaz(:,1),evestaz(:,2),sz,'c','MarkerEdgeColor',...
    [1 1 1], 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)
hold on
scatter(evestaz(:,4),evestaz(:,5),sz,'^','MarkerEdgeColor',...
    [1 1 1], 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)
hold off
if degorutm==111 && hasMT==1
    hold on
    geoshow(coastlat,coastlon);
    xlim([origin(1) origin(1)+nxc*stepgx]);
    ylim([origin(2) origin(2)+nyc*stepgy]);
end
hold off

FName = 'PD_checkerboard';
saveas(PDcheck,fullfile(FPath, FLabel, FName), fformat);

%Checkerboard Qc
Qccheck=figure('Name','Qc checkerboard','NumberTitle',...
    'off','visible',visib,'Position',[300,200,1000,600]);
subplot(1,2,1)
contourf(X,Y,QQchi');
axis equal
view(2)
colormap(gray)

hcb=colorbar;
title(hcb,'Qc^{-1}','FontSize',14,'FontWeight','bold','Color','k');

xlabel('WE','FontSize',12,'FontWeight','bold','Color','k')
ylabel('SN','FontSize',12,'FontWeight','bold','Color','k')
title('Qc checherboard input',...
    'FontSize',12,'FontWeight','bold','Color','k');

hold on
scatter(evestaz(:,1),evestaz(:,2),sz,'c','MarkerEdgeColor',...
    [1 1 1], 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)
hold on
scatter(evestaz(:,4),evestaz(:,5),sz,'^','MarkerEdgeColor',...
    [1 1 1], 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)
if degorutm==111 && hasMT==1
    hold on
    geoshow(coastlat,coastlon);
    xlim([origin(1) origin(1)+nxc*stepgx]);
    ylim([origin(2) origin(2)+nyc*stepgy]);
end
hold off

subplot(1,2,2)
contourf(X,Y,QQcho');
axis equal
view(2)
colormap(gray)

hcb=colorbar;
title(hcb,'Qc^{-1}','FontSize',14,'FontWeight','bold','Color','k');

xlabel('WE','FontSize',12,'FontWeight','bold','Color','k')
ylabel('SN','FontSize',12,'FontWeight','bold','Color','k')
title('Qc checherboard output',...
    'FontSize',12,'FontWeight','bold','Color','k');

hold on
scatter(evestaz(:,1),evestaz(:,2),sz,'c','MarkerEdgeColor',...
    [1 1 1], 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)
hold on
scatter(evestaz(:,4),evestaz(:,5),sz,'^','MarkerEdgeColor',...
    [1 1 1], 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)
hold off
if degorutm==111 && hasMT==1
    hold on
    geoshow(coastlat,coastlon);
    xlim([origin(1) origin(1)+nxc*stepgx]);
    ylim([origin(2) origin(2)+nyc*stepgy]);
end
hold off

FName = 'Qc_checkerboard';
saveas(Qccheck,fullfile(FPath, FLabel, FName), fformat);

%Parameter analysis
pdd =  abs(pd(:,4))>10^(-10);
pdef=pd(pdd,:);
Qcef=Qc(pdd,:);



pdef(:,4)=pdef(:,4)-mean(pdef(:,4));
Qcef(:,4)=Qcef(:,4)-mean(Qcef(:,4));
mipdm=min(pdef(:,4));
mapdm=max(pdef(:,4));
miQcm=min(Qcef(:,4));
maQcm=max(Qcef(:,4));
trepd=0.01*std(pdef(:,4));
treQc=0.01*std(Qcef(:,4));

Qps=Qcef(:,4);
pdps=pdef(:,4);

param_plot=figure('Name','Parameter space separation',...
    'NumberTitle','off','visible',visib,'Position',[300,200,800,600]);

par=pdef(:,1:2);
par(:,3)=pdef(:,3);

c=Qps<-treQc & pdps<-trepd;
par(c,4)=1;
scatter(Qps(c),pdps(c),65,'filled','MarkerFaceColor',[0 0.8 0])
hold on
line([0 0],[mipdm-trepd mapdm+trepd],'Color',[0 0 0],...
    'LineWidth',3)
hold on
line([miQcm-treQc maQcm+treQc],[0 0],'Color',[0 0 0],...
    'LineWidth',3)
hold on
c=Qps<-treQc & pdps>trepd;
par(c,4)=2;
scatter(Qps(c),pdps(c),65,'filled','MarkerFaceColor',[0 0.6 1])
hold on
c=Qps>treQc & pdps<-trepd;
par(c,4)=3;
scatter(Qps(c),pdps(c),65,'filled','MarkerFaceColor',[1 0.6 0])
hold on
c=Qps>treQc & pdps>trepd;
par(c,4)=4;
scatter(Qps(c),pdps(c),65,'filled','MarkerFaceColor',[1 0 0])
hold on
c=(Qps>-treQc & Qps<treQc) | (pdps>-trepd & pdps<trepd);
par(c,4)=0;
scatter(Qps(c),pdps(c),85,'filled','MarkerFaceColor',[0.7 0.7 0.7],...
    'MarkerEdgeColor',[1 1 1],'LineWidth',2)
hold off
xlim([miQcm-treQc maQcm+treQc])
ylim([mipdm-trepd mapdm+trepd])
xlabel('Qc','FontSize',12,'FontWeight','bold','Color','k')
ylabel('Log. peak delay','FontSize',12,'FontWeight','bold','Color','k')
title('Parameter space plot',...
    'FontSize',12,'FontWeight','bold','Color','k');
FName = 'Parameter_space_variations';
saveas(param_plot,fullfile(FPath, FLabel, FName), fformat);

para=pd(:,1:3);
for k=1:length(par(:,1))
    px=par(k,1);
    py=par(k,2);
    pp=par(k,4);
    pf= pd(:,1)==px & pd(:,2)==py;
    para(pf,4)=pp;
end

index=0;
param=zeros(size(QQc));
for i=1:length(x)
    for j=1:length(y)
        index=index+1;
        param(i,j)=para(index,4)-5;
    end
end

mparam=figure('Name','Parameter separation map',...
    'NumberTitle','off','visible',visib,'Position',[300,200,800,600]);
surf(X,Y,param');
view(2)
un_X = unique(param);

fu=find(un_X==-1,1);
if isempty(fu)
    cmap = [0.7 0.7 0.7;  0 0.8 0; 0 0.6 1; 1 0.6 0];
    HTick={'Average','Ls La','Hs La','Ls Ha'};
else
    cmap = [0.7 0.7 0.7;  0 0.8 0; 0 0.6 1; 1 0.6 0; 1 0 0];
    HTick={'Average','Ls La','Hs La','Ls Ha','Hs Ha'};
end
colormap(cmap)

colorbar('Ticks',un_X,'TickLabels',HTick);

hold on
scatter(evestaz(:,1),evestaz(:,2),sz,'c','MarkerEdgeColor',...
    [1 1 1], 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)
hold on
scatter(evestaz(:,4),evestaz(:,5),sz,'^','MarkerEdgeColor',...
    [1 1 1], 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)
if degorutm==111 && hasMT==1
    hold on
    geoshow(coastlat,coastlon);
    xlim([origin(1) origin(1)+nxc*stepgx]);
    ylim([origin(2) origin(2)+nyc*stepgy]);
end
hold off
xlabel('WE','FontSize',12,'FontWeight','bold','Color','k')
ylabel('SN','FontSize',12,'FontWeight','bold','Color','k')
axis square
title('Parameter separation','FontSize',12,'FontWeight','bold',...
    'Color','k');
FName = 'Parameter_map';
saveas(mparam,fullfile(FPath, FLabel, FName), fformat);