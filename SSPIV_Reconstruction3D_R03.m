%%% Scanning Stereo PIV 3Dimensional Reconstruction %%%%%%%%%%%%%%%%%%%%%%%
%%% Originally made by Fernando Zigunov %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Editted and Added onto by Rhylan A. Huss %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% February 28th 2023 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reads in the Davis Vector files (*.vc7) and the Traverse Location file
% (*.mat) before performing a spatio-temporal average on the data set to
% reconstruct the 3-dimensional flow field. The code can also calculate the
% vorticity, divergence, Q-Criterion, and Lambda2 values for the flow
% field. The code can output to TecPlot and/or Paraview. 

% The paper written on the process: 
% https://doi.org/10.1007/s00348-023-03596-w

% Enjoy!

clear; clc; close all;

% load('04_Traverse Locations\TraversePosition_Run5406.mat','laserZPlanes');
load('04_Traverse Locations\TraversePosition_Run5405.mat','laserZPlanes');

[f,p]=uigetfile('*.vc7','Please select the files to reconstruct','Multiselect','on');

%% Define Constants for non-dimensionalizing variables
Winf=0.4*343; %[mps] must match velocity units
D=60; %[mm] must match positional units

%% Opens all files and puts them in a 3D variable
for i=1:length(f)                                                                   
    
    VYavg=loadvec(fullfile(p,f{i}));
    vx=VYavg.vx'; vy=VYavg.vy'; vz=VYavg.vz';
    
    VX(:,:,i)=vx; VY(:,:,i)=vy; VZ(:,:,i)=vz; 
    
%     imagesc(v.x, v.y,v.vz'); colormap(flipud(jet)); title(laserZPlanes(i))
%     daspect([1 1 1]); set(gca,'ydir','normal')
%     caxis([-150 25]);
%     drawnow; pause(0.4)    
    disp(i)
end

%Checks which one is the longer dataset:
if length(laserZPlanes)>length(f)
    disp('Cropping the Z planes to the length of the data...');
    laserZPlanes=laserZPlanes(1:length(f));
else
    disp('Cropping the data to the length of the Z planes...');
    VX=VX(:,:,1:length(laserZPlanes));
    VY=VY(:,:,1:length(laserZPlanes));
    VZ=VZ(:,:,1:length(laserZPlanes));
end

%Builds the coordinate variables
x=VYavg.x + 0;
y=VYavg.y - 24.4; %Only offset from calibration is in the Y direction
z=laserZPlanes + 0;

%% Crops the data to a valid viewing region
%===1-3 (yz) direction===
imagesc(squeeze(VZ(:,165,:))); caxis([-150 25]);colormap jet;
title('Select a CROPPING RECTANGLE for the dataset in the 1-3 direction');
Rcrop=getrect();
Rcrop=round(Rcrop);
%Only uses the Y direction of the crop (Assumes all vector fields included have some useful data
VX=VX(Rcrop(2):(Rcrop(2)+Rcrop(4)),:,:);
VY=VY(Rcrop(2):(Rcrop(2)+Rcrop(4)),:,:);
VZ=VZ(Rcrop(2):(Rcrop(2)+Rcrop(4)),:,:);
y=y(Rcrop(2):(Rcrop(2)+Rcrop(4)));

%===1-2 (xy) direction===
imagesc(VZ(:,:,1)); caxis([-150 25]);colormap jet;
title('Select a CROPPING RECTANGLE for the dataset in the 1-2 direction');
Rcrop=getrect();
Rcrop=round(Rcrop);

VX=VX(Rcrop(2):(Rcrop(2)+Rcrop(4)),Rcrop(1):(Rcrop(1)+Rcrop(3)),:);
VY=VY(Rcrop(2):(Rcrop(2)+Rcrop(4)),Rcrop(1):(Rcrop(1)+Rcrop(3)),:);
VZ=VZ(Rcrop(2):(Rcrop(2)+Rcrop(4)),Rcrop(1):(Rcrop(1)+Rcrop(3)),:);
y=y(Rcrop(2):(Rcrop(2)+Rcrop(4)));
x=x(Rcrop(1):(Rcrop(1)+Rcrop(3)));

%% Now we're ready to build the averaged "blurred" version of the data
%1. Defines the user parameters
zStep=0.6; %Spacing between Z planes in the final data
zBoxSize=3; %Size of the averaging box, mm (+/- zBoxSize/2)
xyOutlierSize=3; %Size of the outlier detection box to improve averages

%2. Turns the zero-vectors into NaN's so we can locate/ignore these data
%points more easily
VX(VX==0)=nan; VY(VY==0)=nan; VZ(VZ==0)=nan;

%3. Starts averaging operation
zNew=min(z):zStep:max(z);
VXavg=nan(size(VX,1),size(VX,2),length(zNew)); VYavg=nan(size(VX,1),size(VX,2),length(zNew)); VZavg=nan(size(VX,1),size(VX,2),length(zNew)); %Original variables created full of NaNs in case some data doesn't get filled

for k=1:length(zNew)
    kVals=find(z>(zNew(k)-zBoxSize/2) & z<(zNew(k)+zBoxSize/2)); %Planes to use for average according to traverse measurements
    %Outlier detection and averaging operation
    for i=1:size(VX,1)
        iVals=(i-floor(xyOutlierSize/2)):(i+floor(xyOutlierSize/2));
        iVals(iVals<1)=[]; iVals(iVals>size(VX,1))=[];
        for j=1:size(VX,2)
            jVals=(j-floor(xyOutlierSize/2)):(j+floor(xyOutlierSize/2));
            jVals(jVals<1)=[]; jVals(jVals>size(VX,2))=[];
            
            %Gets all vector values for this pixel
            vxVals=VX(iVals,jVals,kVals); 
            vyVals=VY(iVals,jVals,kVals); 
            vzVals=VZ(iVals,jVals,kVals); 
            
            notNan=~isnan(vxVals) & ~isnan(vyVals) & ~isnan(vzVals);
            vxVals=vxVals(notNan); vyVals=vyVals(notNan); vzVals=vzVals(notNan);
            
            %Performs outlier rejection on the vectors
            outliers=isoutlier([vxVals vyVals vzVals],'quartiles',1);
            outliers=sum(outliers,2)>0;
            
            %Averages the vectors in the set and records result
            VXavg(i,j,k)=mean(vxVals(~outliers));
            VYavg(i,j,k)=mean(vyVals(~outliers));
            VZavg(i,j,k)=mean(vzVals(~outliers));
            
        end
        %Updates progress tracker
        progPercent=100*i/size(VX,1);
        disp(['z Plane ' num2str(k) '/' num2str(length(zNew)) '; Progress:' num2str(progPercent, '%0.2f') '%']);
        
    end
end

%Creates variables of the coordinates
[X, Y, Z]=meshgrid(x,y,zNew);

%Scans the planes created to get a sense of how it looks like
for k=1:length(zNew)
    [~,Wz]=curl(X(:,:,k),Y(:,:,k),VXavg(:,:,k),VYavg(:,:,k));
    
    imagesc(Wz); caxis([-1 1]*25);colormap redblue;
    drawnow; pause(0.05);
end
stop;

%% Creates a mask before exporting to Tecplot to remove bad vectors
%Recommend to run manually and add/exclude vectors as needed

Mask=isnan(VXavg); %1=Exclude vectors

%===1-3 (yz) direction===
imagesc(squeeze(VXavg(:,80,:))); caxis([-150 25]);colormap jet;title('Select region');
M=roipoly();
MM=reshape(M,[size(M,1) 1 size(M,2)]);
MM=repmat(MM,1,size(Mask,2),1);

Mask=Mask | MM; %MM is an exclusive region
%Mask=Mask & ~MM; %MM is an inclusive region

%===1-2 (xy) direction===
imagesc(squeeze(VZavg(:,:,75))); caxis([-150 25]);colormap jet;title('Select region');
M=roipoly();
MM=repmat(M,1,1,size(Mask,3));

Mask=Mask | MM; %MM is an exclusive region
%Mask=Mask & ~MM; %MM is an inclusive region

%===2-3 (xz) direction===
imagesc(squeeze(nanmean(VZavg,1))); caxis([-150 25]);colormap jet;title('Select region');
M=roipoly();
MM=reshape(M,[1 size(M,1) size(M,2)]);
MM=repmat(MM,size(Mask,1),1,1);

Mask=Mask | MM; %MM is an exclusive region
%Mask=Mask & ~MM; %MM is an inclusive region
   
%% Computes Vorticity

[OmegaX, OmegaY, OmegaZ]=curl(X/D,Y/D,Z/D,VXavg/Winf,VYavg/Winf,VZavg/Winf);
OmegaMag=sqrt(OmegaX.^2+OmegaY.^2+OmegaZ.^2);
Div=divergence(X/D,Y/D,Z/D,VXavg/Winf,VYavg/Winf,-VZavg/Winf);

%% Computes Q-Criterion
% Note the derivatives in this section are for normal cartesian coordinate
% systems. This will not work on Fernando's weird "SkEwEd" system.

disp('Calculating Q-Criterion...')

%Calculates second order accurate derivatives for all components
dudx = zeros(size(VXavg)); dvdx = dudx; dwdx = dudx; dudy = dudx; 
dvdy = dudx; dwdy = dudx; dudz = dudx; dvdz = dudx; dwdz = dudx;

dudx(:,2:end-1,:) = (VXavg(:,3:end,:)-VXavg(:,1:end-2,:))./(X(:,3:end,:)-X(:,1:end-2,:)) * D/Winf;
dvdx(:,2:end-1,:) = (VYavg(:,3:end,:)-VYavg(:,1:end-2,:))./(X(:,3:end,:)-X(:,1:end-2,:)) * D/Winf;
dwdx(:,2:end-1,:) = (VZavg(:,3:end,:)-VZavg(:,1:end-2,:))./(X(:,3:end,:)-X(:,1:end-2,:)) * D/Winf;

dudy(2:end-1,:,:) = (VXavg(3:end,:,:)-VXavg(1:end-2,:,:))./(Y(3:end,:,:)-Y(1:end-2,:,:)) * D/Winf;
dvdy(2:end-1,:,:) = (VYavg(3:end,:,:)-VYavg(1:end-2,:,:))./(Y(3:end,:,:)-Y(1:end-2,:,:)) * D/Winf;
dwdy(2:end-1,:,:) = (VZavg(3:end,:,:)-VZavg(1:end-2,:,:))./(Y(3:end,:,:)-Y(1:end-2,:,:)) * D/Winf;

dudz(:,:,2:end-1) = (VXavg(:,:,3:end)-VXavg(:,:,1:end-2))./(Z(:,:,3:end)-Z(:,:,1:end-2)) * D/Winf;
dvdz(:,:,2:end-1) = (VYavg(:,:,3:end)-VYavg(:,:,1:end-2))./(Z(:,:,3:end)-Z(:,:,1:end-2)) * D/Winf;
dwdz(:,:,2:end-1) = (VZavg(:,:,3:end)-VZavg(:,:,1:end-2))./(Z(:,:,3:end)-Z(:,:,1:end-2)) * D/Winf;

Q = -0.5*(dudx.^2 + dvdy.^2 + dwdz.^2 + 2*dudy.*dvdx + 2*dudz.*dwdx + 2*dvdz.*dwdy);

disp('Q-Criterion Calculated.')

%% Calculate Lambda2 Values
disp('Calculating Lambda2...')

% Preallocate because we are good coders :)
lambda2 = NaN(size(VXavg)); 
A = NaN(3);

for i=1:size(VXavg,1)
    fprintf('i index %i out of %i\n',i,size(VXavg,1));
    for j=1:size(VXavg,2)
        for k=1:size(VXavg,3)
                    A = [...
                        dudx(i,j,k) dudy(i,j,k) dudz(i,j,k); 
                        dvdx(i,j,k) dvdy(i,j,k) dvdz(i,j,k); 
                        dwdx(i,j,k) dwdy(i,j,k) dwdz(i,j,k)];
                    
                    SR = (A +A')/2; %Strain rate tensor
                    OR = (A -A')/2; %Vorticity tensor

                    try
                        ls = sort(eig(SR^2 +OR^2));
                        lambda2(i,j,k) = ls(2);  %\lambda_2-criterion
                    catch
                        lambda2(i,j,k) = NaN;
                    end
        end
    end
end

disp('Lambda2 Calculated')

%% Saves to MatLab
%Saves the memory into a mat file for future review
save('SSPIV_Run5405_M0p4_02.mat', '-v7.3')

%% Saves to Tecplot
disp('Writing to *.plt file');
tdata.Nvar = 13; %Sets number of variables
tdata.varnames = {'x','y','z','u','v','w','omegax','omegay','omegaz','div','Q','lambda2','b'}; %Sets variable names
tdata.vformat(1:13) = 1; %Sets variables to float

tdata.cubes(1).zonename= 'SSPIV' ;
%Data input (See mat2tecplot.m)
tdata.cubes(1).x=X;
tdata.cubes(1).y=Y;
tdata.cubes(1).z=Z;
tdata.cubes(1).v(1,:,:,:)= VXavg;
tdata.cubes(1).v(2,:,:,:)= VYavg;
tdata.cubes(1).v(3,:,:,:)= -VZavg; %Negative sign due to calibration direction
tdata.cubes(1).v(4,:,:,:)= OmegaX;
tdata.cubes(1).v(5,:,:,:)= OmegaY;
tdata.cubes(1).v(6,:,:,:)= OmegaZ;
tdata.cubes(1).v(7,:,:,:)= Div;
tdata.cubes(1).v(8,:,:,:)= Q;
tdata.cubes(1).v(9,:,:,:)= lambda2;
tdata.cubes(1).v(10,:,:,:)= Mask;

BINARYfilename = 'SSPIV_Run5406_M0p3_02.plt';
mat2tecplot(tdata,BINARYfilename);
    
disp('*.plt File written.');  

%% Saves to Paraview 
disp('Writing to *.vtk file');

VTKfilename = 'SSPIV_Run5405_M0p4_02.vtk';

% VXavg(isnan(VXavg))=0; VYavg(isnan(VYavg))=0; VZavg(isnan(VZavg))=0;
% OmegaX(isnan(OmegaX))=0; OmegaY(isnan(OmegaY))=0; OmegaZ(isnan(OmegaZ))=0; OmegaMag(isnan(OmegaMag))=0;
% Div(isnan(Div))=0; 
% Q(isnan(Q))=0; lambda2(isnan(lambda2))=0;

se=strel('sphere',2);
Mask2=imdilate(Mask,se);

vtkwrite(VTKfilename, 'structured_grid', X, Y, Z, ...
'vectors', 'V', VXavg, VYavg, -VZavg,...
'scalars', 'OmegaX', OmegaX,...
'scalars', 'OmegaY', OmegaY,...
'scalars', 'OmegaZ', OmegaZ,...
'scalars', 'Div', Div,...
'scalars', 'qcriterion',Q,...
'scalars', 'lambda2',lambda2,...
'scalars', 'Mask', Mask2,'binary')

disp('*.vtk File written.');  

