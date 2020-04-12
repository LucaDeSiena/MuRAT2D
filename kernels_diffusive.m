% script_grid_tomography_coda_v1.m
%
% Script for computing the K integral (eq 5) for a 3D grid
% to be used for coda tomography
% This script will write two result files (plain text):
% An ascii file describing the coordinates of the grid
% An ascii file with the K function estimated at each point of the grid
%
% Execution time is also reported in this script
%
% Estructure of the file describing the coordinates:
%     Nx    -> number of x values in the grid
%     Ny    -> number of y values in the grid
%     Nz    -> number of z values in the grid
%     x(1)  -> first x coordinate
%     x(2)  ...
%     x(3)
%     ...
%     x(Nx) -> last x coordinate
%     y(1)  -> first y coordinate
%     y(2)
%     ...
%     y(Ny) -> last y coordinate
%     z(1)  -> first z coordinate
%     z(2)
%     ...
%     z(Nz) -> last z coordinate
%
% Estructure of the file with the K function results
%     K(1)      1st point for x(1),y(1),z(1)
%     K(2)      2nd point for x(2),y(1),z(1)
%     K(3)      3rd point for x(3),y(1),z(1)
%     ....
%     K(Nx * Ny * Nz)  last point for x(Nx),y(Ny),z(Nz)
%
% The FLAGS DRAW_FIGURES and DRAW_FINAL_FIGURES allow enabling or disabling drawing figures 
% The FLAG WRITE_FILE_RESULTS enables or disables writing the results in the files 
%
% Values to be configured:
% FLAGS
% FILENAME_BASE to define the name of the output filenames
% Geophysics constants: v, B0, Le_1
% Resolution in space, resulution in time, Tmax: DR,DT,T
% Limits of the spatial grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K_grid,r_grid1]=...
    kernels_diffusive(T,even,staz,XY,degorutm,v,kT)
%% FLAGS
exec_time(1)=now;

DRAW_FIGURES = 0;    % Flag for drawing or not figures
DRAW_FINAL_FIGURES = 0; % Flag for drawing only final figures (out of estimation of execution time)
WRITE_FILE_RESULTS = 0; % Flag for writing or not results in files
% Note that drawing plots and writing results consumes sinificant
% computation time!

%% CONFIGURATION OF THE SCRIPT
% Output filenames
% Geophysical constants
B0 = 0.5;%0.72; % constant Paasschen function- Del Pezzo et al. 1995, Campi Flegrei
Le_1 = 0.02;%0.35; % constant Paasschen function- Del Pezzo et al. 1995, Campi Flegrei
DT=0.5;
x1=unique(XY(:,1));
stepgx=x1(2)-x1(1);
stepgy=XY(2,2)-XY(1,2);

if even(3)<0
    even(3)=abs(even(3));
end

if degorutm==1
    metordeg=1/1000;
    metordegz=metordeg;
elseif degorutm==111
    metordeg=111;
    metordegz=1;
end

DR=(stepgx+stepgy)/kT*metordeg; % Spacing
%origin of the grid where we compute kernels

origin=[(even(1)+staz(1))/2 (even(2)+staz(2))/2 0]; 
S=v*T+10;

xs=(even(1)-origin(1))*metordeg;
ys=(even(2)-origin(2))*metordeg;
zs=even(3)*metordegz; 
xr=(staz(1)-origin(1))*metordeg;
yr=(staz(2)-origin(2))*metordeg;
zr=0;

%% INITIAL COMPUTATIONS
% Distance between source and receiver
dx=xr-xs; dy=yr-ys; dz=zr-zs;
D0 = sqrt(dx^2 + dy^2 + dz^2);

% GRID DEFINITION
x_grid=-S+DR:DR:S;  % vector of positions x (in km)
y_grid=-S+DR:DR:S;  % vector of positions y (in km)
z_grid=0:DR:60; % vector of positions z (in km) from DR km a.s.l.
Nx=length(x_grid); % number of positions x
Ny=length(y_grid); % number of positions y
Nz=length(z_grid); % number of positions z
Nxyz=Nx*Ny*Nz;
if Nxyz>50e6
    fprintf('Warning: a lot of points are going to be included in the grid (%d points)\n',Nx*Ny*Nz)
    fprintf('Press any key to continue, or CRTL-C to abort\n')
    pause
end    
r_grid=zeros(Nx*Ny*Nz,3);
idx=0;
for ix=1:Nx
    for iy=1:Ny
        for iz=1:Nz
            idx=idx+1;
            r_grid(idx,1)=x_grid(ix);
            r_grid(idx,2)=y_grid(iy);
            r_grid(idx,3)=z_grid(iz);
        end
    end
end

% Distances r1 and r2, from each point of the grid to source or to receptor
x=r_grid(:,1); y=r_grid(:,2); z=r_grid(:,3);
d1x=x-xs; d1y=y-ys; d1z=z-zs;
d2x=x-xr; d2y=y-yr; d2z=z-zr;
r1_grid=sqrt(d1x.^2+d1y.^2+d1z.^2);
r2_grid=sqrt(d2x.^2+d2y.^2+d2z.^2);

cond_ellipsoid=abs(r1_grid+r2_grid-T*v)<=DR;   % elements of the grid in the ellipsoid
cond_interior=r1_grid+r2_grid<T*v;%-1.6*DR;      % elements of the grid inside the ellipsoid
exec_time(2)=now;

%% CATALOG OF PAASSCHENS FUNCTIONS
R=DR:DR:T*v;    % vector of distances
NT=ceil((T+DT)/DT); % max number of points in coda
NR=length(R);   % number of distances

t0_R=zeros(NR,1); % to store t0 for each r
A_R=zeros(NR,1); % to store A for each r
N_R=zeros(NR,1);   % to store the number of points in coda for each r
coda_R=zeros(NR,NT); % to store the coda for each r

for idx_r=1:NR
    r=R(idx_r);
    [t0,A,N,coda]=paasschens_function_v1(r,v,B0,Le_1,DT,T);
    t0_R(idx_r)=t0;
    A_R(idx_r)=A;
    N_R(idx_r)=N;
    coda_R(idx_r,1:N)=coda(1:N);
    
%     if mod(idx_r,100)==0
%         fprintf('Finished E(r,t) for r=%f (%d of %d)\n',r,idx_r,NR);
%     end
end
exec_time(3)=now;

if DRAW_FIGURES
    figure(1);
    imagesc(0:DT:T,R,log10(coda_R));
    set(gca,'ydir','normal');
    xlabel('time (s)');
    ylabel('distance (km)')
    title('Paasschens catalog: normalized coda (shifted in time) for each r')
end

%% CATALOG OF r1 - r2 CONVOLUTIONS EVALUATED AT T
R1=R;
R2=R;
coda_coda=zeros(NR,NR);
coda_delta=zeros(NR,NR);
tau_vect=(0:(NT-1))*DT;
for idx_r1=1:NR
    for idx_r2=1:idx_r1
        r1=R1(idx_r1);
        r2=R2(idx_r2);
        % delta_T is the overlapping interval between E(r1,t) and  E(r2,T-r)
        delta_T=T-t0_R(idx_r1) - t0_R(idx_r2);
        % N_mues is the number of samples involved in the overlapping interval
        N_mues = round(delta_T/DT);
        if r1+r2>=(D0-2*DR) && r1+r2<(T+DT)*v && N_mues>=1
            
            % Computing the contributions (coda1 x delta2) and (delta1 x coda 2)
            %coda1=A_R(idx_r1)*A_coda_R(idx_r2)*interp1(tau_vect,coda_R(idx_r2,:),delta_T);
            %coda2=A_R(idx_r2)*A_coda_R(idx_r1)*interp1(tau_vect,coda_R(idx_r1,:),delta_T);
            coda1=A_R(idx_r1)*coda_R(idx_r2,N_mues+1); % Important: N_mues+1
            coda2=A_R(idx_r2)*coda_R(idx_r1,N_mues+1);
            coda_delta(idx_r1,idx_r2)= coda1+coda2;
            
            % Computing the contribution of (coda1 x coda2)
            % integral coda1(u) coda(T-u) du            
            x1=coda_R(idx_r1,1:N_mues+1);
            x2=coda_R(idx_r2,N_mues+1:-1:1);            
            coda_coda(idx_r1,idx_r2)=sum(x1.*x2)*DT;
            
            % simmetrical calculations
            coda_delta(idx_r2,idx_r1)=coda_delta(idx_r1,idx_r2);
            coda_coda(idx_r2,idx_r1)=coda_coda(idx_r1,idx_r2);
        end
    end
%     if mod(idx_r1,100)==0
%         fprintf('Finished Ks(r1,r2,T) for r1=%f (%d of %d)\n',r1,idx_r1,NR);
%     end
    
end
K_integral=coda_coda+coda_delta;

%% Plots and verification of results
if DRAW_FIGURES
    figure(2)
    imagesc(R1,R2,log10(coda_coda));
    set(gca,'ydir','normal');
    xlabel('distance r1 (km)');
    ylabel('distance r2 (km)')
    title('K integral at T: contribution coda1-coda2')
    
    figure(3)
    imagesc(R1,R2,log10(coda_delta));
    set(gca,'ydir','normal');
    xlabel('distance r1 (km)');
    ylabel('distance r2 (km)')
    title('K integral at T: contribution coda1-delta2 + coda2-delta1')
    
    figure(4)
    imagesc(R1,R2,log10(K_integral));
    set(gca,'ydir','normal');
    xlabel('distance r1 (km)');
    ylabel('distance r2 (km)')
    title('K integral at T: both contributions')
    
    figure(5)
    semilogy(R2,coda_coda(8,:),R2,coda_delta(8,:),R2,K_integral(8,:));
    legend('coda-coda','coda-delta','total')
    title(sprintf('K integral for r1=%.3f km as a function of r2',R1(8)))
    xlabel('distance r2 (km)');
    ylabel('K integral')
end

% Results for r1=R(idx_r1) r2=R(idx_r2)
% idx_r1=27; idx_r2=12;
% fprintf('r1=%.2f km   r2=%.2f km   T=%.3f s   K_integ=%g   dt=%f\n',...
%     R(idx_r1),R2(idx_r2),T,K_integral(idx_r1,idx_r2),DT)

% computing K(x,y,z) by linear interpolation for those points inside the ellipsoid
c1=cond_ellipsoid; c2=cond_interior;
K_grid=zeros(Nx*Ny*Nz,1);
exec_time(4)=now;
K_grid(c2)=interp2(R1,R2,K_integral,r1_grid(c2),r2_grid(c2));
exec_time(5)=now;
K_grid(isnan(K_grid))=max(K_grid);
r_grid1=zeros(Nx*Ny*Nz,3);
r_grid1(:,1)=origin(1)+(r_grid(:,1)-DR)/metordeg;
r_grid1(:,2)=origin(2)+(r_grid(:,2)-DR)/metordeg;
r_grid1(:,3)=-r_grid(:,3);
x_grid1=r_grid1(:,1);  % vector of positions x (in km)
y_grid1=r_grid1(:,2);  % vector of positions x (in km)

% PENDING TASK: computing deltas at the ellipsoid (but this is analytical: easy)
% Finally, as discussed, this contribution should not be included in the integral.
% Discussed by Edoardo del Pezzo and Angel de la Torre, Granada, Oct 26th and 27th, 2016

%% WRITING RESULTS IN OUTPUT FILES
if WRITE_FILE_RESULTS
    file_K=sprintf('%s_K_values.txt',FILENAME_BASE);
    file_coordinates=sprintf('%s_coordinates.txt',FILENAME_BASE);
    f_out=fopen(file_K,'wt');
    fprintf(f_out,'%g\n',K_grid);
    fclose(f_out);
    f_out=fopen(file_coordinates,'wt');
    fprintf(f_out,'%d\n',[Nx Ny Nz]);
    fprintf(f_out,'%f\n',x_grid);
    fprintf(f_out,'%f\n',y_grid);
    fprintf(f_out,'%f\n',z_grid);
    fclose(f_out);
end

if DRAW_FIGURES
    %drawing the ellipsoid, the source, the receptor and the grid points inside the ellipsoid
    figure(6)
    plot3(r_grid(c2,1),r_grid(c2,2),r_grid(c2,3),'.b',r_grid(c1,1),r_grid(c1,2),r_grid(c1,3),'.r',...
        xs,ys,zs,'*k',xr,yr,zr,'ok') 
end

exec_time(6)=now;   % estimation of execution time finishes here
Exec_time=diff(exec_time)/24/60/60;
% fprintf('Execution time = %.3f seconds (preparing grid)\n',Exec_time(1))
% fprintf('Execution time = %.3f seconds (catalog E(r,t) - Paasschens)\n',Exec_time(2))
% fprintf('Execution time = %.3f seconds (catalog K(r1,r2) - convolutions)\n',Exec_time(3))
% fprintf('Execution time = %.3f seconds (interpolation)\n',Exec_time(4))
% fprintf('Execution time = %.3f seconds (writing results)\n',Exec_time(5))
% fprintf('         Total = %.3f seconds\n',sum(Exec_time))

%% DRAWING FINAL RESULTS
if DRAW_FINAL_FIGURES
    c_z0=r_grid1(:,3)==0;  % condition for belonging to plane z=0
%     c_z0=r_grid1(:,3)==0;  % condition for belonging to plane z=0
    K_plane=K_grid(c_z0); % extracting data corresponding to the plane
    x_plane=r_grid1(c_z0,1); 
    y_plane=r_grid1(c_z0,2);
    K_plane=reshape(K_plane,Ny,Nx); % rearranging data
    x_plane=reshape(x_plane,Ny,Nx);
    y_plane=reshape(y_plane,Ny,Nx);
    
    figure(7); clf
    imagesc(x_grid1,y_grid1,log(K_plane))
    set(gca,'ydir','normal'); axis equal
    hcb=colorbar;
    title(hcb,'Ka (s/km^3)')
    xlabel('WE')
    ylabel('SN')
    figure(8); clf
    surf(x_plane,y_plane,log(K_plane))
    colorbar
    
end

return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%