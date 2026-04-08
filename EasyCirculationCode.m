%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Quick Circulation Grab %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Rhylan A Huss %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% October 7th 2025 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Follow the instructions in this shit.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% Declaration of the Random

%Variables
D = 146.05; % mm 
phi = 45;
nu = 1.5339e-05; %m^2/s

% Colormaps
mustang = [1 1 1; jet];
temp = redblue(540);
blue = [temp(1:270,:)];
red = [temp(270:end,:)];

% Default Figure Formats
set(groot,'defaultAxesFontName','Times New Roman')
set(groot,'defaultAxesFontSize',14)
set(groot,'defaultAxesFontWeight','Bold')
set(groot,'defaultFigurePosition',[100 100 700 600])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Start Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% CHANGE THE FILE PATH HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path = {...
    'A:\02 - Cavity&Door\SharpEdgeModel_shoaib1\Baseline_70K\SubOverTimeAvg_sL=all_01\StereoPIV_MP(3x32x32_75%ov)_GPU_09corr'};

% Freestream Velocity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Winf = [...
    02.63;
    03.68;
    04.19;
    04.73;
    07.35;
    12.60;
    15.75;
    16.80;
    17.85;
    19.95;
    21.01]; % m/s

% Reynolds Number %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ReD = Winf*(D/1000)/nu;

%% Load PIV Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load and Average the data set %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear vxi vyi uvz uvy
files = dir(fullfile(path{1}, '*.vc7'));

% Preallocates Arrays %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic; Vec = loadvec(fullfile(path{1},files(1).name)); 
TotalFiles = min(1000,length(files));
vxi = NaN([size(Vec.vx'),TotalFiles]); 
vyi = NaN([size(Vec.vx'),TotalFiles]); 
vzi = NaN([size(Vec.vx'),TotalFiles]); 
Counter = zeros(size(Vec.vx'));

% Loop Through .vc7 files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Reading in ReD: %iK\nPath: %s\n\n',round(ReD(1),-4)/1000,path{1});
for I=1:TotalFiles %Loopity Loop that dataset
  
    % Load in Velocity Components and Uncertainties %%%%%%%%%%%%%%%%%%%
    IMX = readimx(fullfile(path{1},files(I).name)); Frames = IMX.Frames{1,1}; 
    Components = Frames.ComponentNames;
    if I == 1 
        U0i = findComponent(Components,'U0');
        V0i = findComponent(Components,'V0');
        W0i = findComponent(Components,'W0');
        CVi = findComponent(Components,'TS:Correlation value');
        Eni  = findComponent(Components,'ENABLED');
        try UVi = findComponent(Components,'TS:Uncertainty V'); catch; UVi = NaN; end
    end
    % ADDITIONAL POST-PROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Enabled (Geometric Mask)
    En = ((Frames.Components{Eni,1}.Planes{1,1} * Frames.Components{Eni,1}.Scale.Slope) + Frames.Components{Eni,1}.Scale.Offset)';
    % Correlation Value
    CV = ((Frames.Components{CVi,1}.Planes{1,1} * Frames.Components{CVi,1}.Scale.Slope) + Frames.Components{CVi,1}.Scale.Offset)';
    % Vector Uncertainty
    if ~isnan(UVi); UV = ((Frames.Components{UVi,1}.Planes{1,1} * Frames.Components{UVi,1}.Scale.Slope) + Frames.Components{UVi,1}.Scale.Offset)';
    else; UV = ones(size(CV)); end

    % Load in Component Vectors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vxi(:,:,I)= -((Frames.Components{U0i,1}.Planes{1,1} * Frames.Components{U0i,1}.Scale.Slope) + Frames.Components{U0i,1}.Scale.Offset)'; 
    vyi(:,:,I)= -((Frames.Components{V0i,1}.Planes{1,1} * Frames.Components{V0i,1}.Scale.Slope) + Frames.Components{V0i,1}.Scale.Offset)'; 
    vzi(:,:,I)= -((Frames.Components{W0i,1}.Planes{1,1} * Frames.Components{W0i,1}.Scale.Slope) + Frames.Components{W0i,1}.Scale.Offset)'; 

    if I == 1
        % Rebuild Coordinate System %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        x = -(Vec.x); y = -(Vec.y);
        x=(x - 0); y=y - (99.225 + 51.5); % (76.2 + 51.5)
    end

    % Define NaN Mask
    M = (double(En)); M(M==0) = NaN; %double(CV > 0.2).*
    vxi(:,:,I) =  M .* vxi(:,:,I); 
    vyi(:,:,I) =  M .* vyi(:,:,I);
    vzi(:,:,I) =  M .* vzi(:,:,I);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Gaussian Smoothing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    filtWidth = 5;
    filtSigma = 1.5;
    imageFilter=fspecial('gaussian',filtWidth,filtSigma);
    vxWorking = nanconv(vxi(:,:,I),imageFilter, 'nanout');
    vyWorking = nanconv(vyi(:,:,I),imageFilter, 'nanout');
    vzWorking = nanconv(vzi(:,:,I),imageFilter, 'nanout');
    vxWorking(vxWorking==0) = NaN; vyWorking(vyWorking==0) = NaN; vzWorking(vzWorking==0) = NaN;

    % Finalize Output Vectors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vxi(:,:,I) = vxWorking; vyi(:,:,I) = vyWorking; vzi(:,:,I) = vzWorking;

    % Display Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clc; 
    fprintf('Reading in ReD: %iK\nPath: %s\n\n',round(ReD(1),-4)/1000,path{1});
    fprintf('Read File: %s\n',files(I).name);
    fprintf('Vector Plot %i/%i\n',I,TotalFiles);
    fprintf('Time Elapsed: %.0f of %.0f minutes\n',toc/60,(TotalFiles/I) * (toc/60));
end

% Compute Average Excluding NaN's (Note Matlabs 'omitnan' is bogus)
vxWorking = vxi; vyWorking = vyi; vzWorking = vzi;
vxWorking(isnan(vxWorking)) = 0; vyWorking(isnan(vyWorking)) = 0; vzWorking(isnan(vzWorking)) = 0;
ActiveVx = sum(~isnan(vxi),3); ActiveVy = sum(~isnan(vyi),3); ActiveVz = sum(~isnan(vzi),3); 
vxWorking = sum(vxWorking,3)./ActiveVx; vyWorking = sum(vyWorking,3)./ActiveVy; vzWorking = sum(vzWorking,3)./ActiveVz;

% Fill Pot Holes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vxAvg = inpaint_nans(vxWorking,2);
vyAvg = inpaint_nans(vyWorking,2);
vzAvg = inpaint_nans(vzWorking,2);
vMag = sqrt(vxAvg.^2 + vyAvg.^2 + vzAvg.^2);

% Find Correct Velocity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,K] = min(abs(Winf - mean(vMag,'all')));
Winf = Winf(K);

%% Calculate Vorticity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X,Y] = meshgrid(x,y);
clipMask = (X/D<0.6).*(X/D>-0.6).*(Y/D>-1.35).*(Y/D<-0.66);
[vortz] = VortZ(X/D,Y/D,vxAvg/Winf(1),vyAvg/Winf(1));
vortz = clipMask.*vortz;

%% Find Vortices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The Left Vortex Now
figure(69); clf; hold on;
contourf(vortz,[-30:2:30],'edgecolor','none')
clim([-30 30]); colormap redblue
daspect([1 1 1]);

title('Select POSITIVE Vortex');
clc; fprintf('Select POSITIVE Vortex...\n');
roiLeft = drawrectangle;
roiLeft = round(roiLeft.Position);

jIndex = roiLeft(2):roiLeft(2)+roiLeft(4);
iIndex = roiLeft(1):roiLeft(1)+roiLeft(3);

% Find Left Vortex Center
[vortexLeft,G1] = IDVortex_Gamma1(-x(iIndex),-y(jIndex),...
    vxAvg(jIndex,iIndex),vyAvg(jIndex,iIndex));
G1Left = max(G1,[],'all'); vortexLeft = -vortexLeft;

% The Right Vortex Now
figure(69); clf; hold on;
contourf(vortz,[-30:2:30],'edgecolor','none')
clim([-30 30]); colormap redblue
daspect([1 1 1]);

title('Select NEGATIVE Vortex');
clc; fprintf('Select NEGATIVE Vortex...\n');
roiRight = drawrectangle;
roiRight = round(roiRight.Position);

jIndex = roiRight(2):roiRight(2)+roiRight(4);
iIndex = roiRight(1):roiRight(1)+roiRight(3);

% Find Right Vortex Center
[vortexRight,G1] = IDVortex_Gamma1(-x(iIndex),-y(jIndex),...
    vxAvg(jIndex,iIndex),vyAvg(jIndex,iIndex));
G1Right = max(G1,[],'all'); vortexRight = -vortexRight;

figure(69); clf; hold on;
contourf(x/D,y/D,vortz,[-30:2:30],'edgecolor','none')
n = 5; quiver(x(1:n:end)/D,y(1:n:end)/D,vxAvg(1:n:end,1:n:end),vyAvg(1:n:end,1:n:end))
clim([-30 30]); colormap redblue
xlim([-0.6 0.6]); ylim([-1 -0.35]); daspect([1 1 1]);
xlabel('x/D'); ylabel('y/D');

plot(vortexLeft(1)/D,vortexLeft(2)/D,'ro')
plot(vortexRight(1)/D,vortexRight(2)/D,'ro')

%% Left Calculate Circulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate Circulation
[~,i] = min(abs(x-vortexLeft(1)),[],'all');
[~,j] = min(abs(y-vortexLeft(2)),[],'all');

leftHalf = vortz>=0;

rim = 0;

dx = x(2) - x(1); %mm
dy = y(2) - y(1); %mm

counter = 1;
Circulation = NaN(50,2);
while counter <= 2 || (Circulation(counter,2) - Circulation(counter-1,2))/Circulation(counter,2) > 0.03
    counter = counter+1;
    rim = rim + abs(dx);
    tempMaskLeft = zeros(size(vortz));
    for R = 1:length(reshape(vortz,[],1))
        if sqrt((X(R) - vortexLeft(1))^2 + (Y(R) - vortexLeft(2))^2) <= rim^2
            tempMaskLeft(R) = 1;
        end
    end
    Circulation(counter,:) = [rim/D,sum(leftHalf.*tempMaskLeft.*vortz,'all','omitnan')*dx*dy/D^2];
    figure(6); clf; hold on;
    contourf(x,y,leftHalf.*tempMaskLeft.*vortz);
    pause(0.1);
end

R0Left = max(abs(Circulation(:,1)));
CirculationLeft(1) = max(abs(Circulation(:,2)));

clc;
fprintf('Left Vortex:\nR0: %f\nCirculation: %f\n\n',R0Left,CirculationLeft(1));

%% Right Calculate Circulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate Circulation
[~,i] = min(abs(x-vortexRight(1)),[],'all');
[~,j] = min(abs(y-vortexRight(2)),[],'all');

rightHalf =vortz<=0;

rim = 1;

dx = x(2) - x(1);
dy = y(2) - y(1);

counter = 1;
Circulation = NaN(50,2);
while counter <= 2 || (Circulation(counter,2) - Circulation(counter-1,2))/Circulation(counter,2) > 0.03
    counter = counter+1;
    rim = rim + abs(dx);
    tempMaskRight = zeros(size(vortz));
    for R = 1:length(reshape(vortz,[],1))
        if sqrt((X(R) - vortexRight(1))^2 + (Y(R) - vortexRight(2))^2) <= rim^2
            tempMaskRight(R) = 1;
        end
    end
    Circulation(counter,:) = [rim/D,sum(rightHalf.*tempMaskRight.*vortz,'all','omitnan')*dx*dy/D^2];
    figure(6); clf; hold on;
    contourf(x,y,rightHalf.*tempMaskRight.*vortz);
    pause(0.1);
end

R0Right = max(abs(Circulation(:,1)));
CirculationRight(1) = max(abs(Circulation(:,2)));

fprintf('Right Vortex:\nR0: %f\nCirculation: %f\n\n',R0Right,CirculationRight(1));

%% Consolidate Circulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gamma(1) = (abs(CirculationLeft(1)) + abs(CirculationRight(1)))/2;
fprintf('Total Circulation: %f\n\n',Gamma(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate Vorticity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vortz] = VortZ(X,Y,Vx,Vy)
    % Calculate Vorticity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dx = X(1,2)-X(1,1); dy = Y(2,1)-Y(1,1);
    % Calculate Gradients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dudy = zeros(size(Vx)); dvdx = zeros(size(Vx));
    dudy(2:end-1,:) = (Vx(3:end,:)-Vx(1:end-2,:))./(2*dy);
    dvdx(:,2:end-1) = (Vy(:,3:end)-Vy(:,1:end-2))./(2*dx);
    % Output Vorticity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vortz = dvdx - dudy;
end

% Find Component in .im7 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function index = findComponent(Components,String)
    for i = 1:length(Components)
        if contains(Components{i,1},String)
            index = i;
            break
        end
    end
end

% NaN Convolution to Fill Holes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c = nanconv(a, k, varargin)
    % NANCONV Convolution in 1D or 2D ignoring NaNs.
    %   C = NANCONV(A, K) convolves A and K, correcting for any NaN values
    %   in the input vector A. The result is the same size as A (as though you
    %   called 'conv' or 'conv2' with the 'same' shape).
    %
    %   C = NANCONV(A, K, 'param1', 'param2', ...) specifies one or more of the following:
    %     'edge'     - Apply edge correction to the output.
    %     'noedge'   - Do not apply edge correction to the output (default).
    %     'nanout'   - The result C should have NaNs in the same places as A.
    %     'nonanout' - The result C should have ignored NaNs removed (default).
    %                  Even with this option, C will have NaN values where the
    %                  number of consecutive NaNs is too large to ignore.
    %     '2d'       - Treat the input vectors as 2D matrices (default).
    %     '1d'       - Treat the input vectors as 1D vectors.
    %                  This option only matters if 'a' or 'k' is a row vector,
    %                  and the other is a column vector. Otherwise, this
    %                  option has no effect.
    %
    %   NANCONV works by running 'conv2' either two or three times. The first
    %   time is run on the original input signals A and K, except all the
    %   NaN values in A are replaced with zeros. The 'same' input argument is
    %   used so the output is the same size as A. The second convolution is
    %   done between a matrix the same size as A, except with zeros wherever
    %   there is a NaN value in A, and ones everywhere else. The output from
    %   the first convolution is normalized by the output from the second 
    %   convolution. This corrects for missing (NaN) values in A, but it has
    %   the side effect of correcting for edge effects due to the assumption of
    %   zero padding during convolution. When the optional 'noedge' parameter
    %   is included, the convolution is run a third time, this time on a matrix
    %   of all ones the same size as A. The output from this third convolution
    %   is used to restore the edge effects. The 'noedge' parameter is enabled
    %   by default so that the output from 'nanconv' is identical to the output
    %   from 'conv2' when the input argument A has no NaN values.
    %
    % See also conv, conv2
    %
    % AUTHOR: Benjamin Kraus (bkraus@bu.edu, ben@benkraus.com)
    % Copyright (c) 2013, Benjamin Kraus
    % $Id: nanconv.m 4861 2013-05-27 03:16:22Z bkraus $
    
    % Process input arguments
    for arg = 1:nargin-2
        switch lower(varargin{arg})
            case 'edge'; edge = true; % Apply edge correction
            case 'noedge'; edge = false; % Do not apply edge correction
            case {'same','full','valid'}; shape = varargin{arg}; % Specify shape
            case 'nanout'; nanout = true; % Include original NaNs in the output.
            case 'nonanout'; nanout = false; % Do not include NaNs in the output.
            case {'2d','is2d'}; is1D = false; % Treat the input as 2D
            case {'1d','is1d'}; is1D = true; % Treat the input as 1D
        end
    end
    
    % Apply default options when necessary.
    if(exist('edge','var')~=1); edge = false; end
    if(exist('nanout','var')~=1); nanout = false; end
    if(exist('is1D','var')~=1); is1D = false; end
    if(exist('shape','var')~=1); shape = 'same';
    elseif(~strcmp(shape,'same'))
        error([mfilename ':NotImplemented'],'Shape ''%s'' not implemented',shape);
    end
    
    % Get the size of 'a' for use later.
    sza = size(a);
    
    % If 1D, then convert them both to columns.
    % This modification only matters if 'a' or 'k' is a row vector, and the
    % other is a column vector. Otherwise, this argument has no effect.
    if(is1D);
        if(~isvector(a) || ~isvector(k))
            error('MATLAB:conv:AorBNotVector','A and B must be vectors.');
        end
        a = a(:); k = k(:);
    end
    
    % Flat function for comparison.
    o = ones(size(a));
    
    % Flat function with NaNs for comparison.
    on = ones(size(a));
    
    % Find all the NaNs in the input.
    n = isnan(a);
    
    % Replace NaNs with zero, both in 'a' and 'on'.
    a(n) = 0;
    on(n) = 0;
    
    % Check that the filter does not have NaNs.
    if(any(isnan(k)));
        error([mfilename ':NaNinFilter'],'Filter (k) contains NaN values.');
    end
    
    % Calculate what a 'flat' function looks like after convolution.
    if(any(n(:)) || edge)
        flat = conv2(on,k,shape);
    else flat = o;
    end
    
    % The line above will automatically include a correction for edge effects,
    % so remove that correction if the user does not want it.
    if(any(n(:)) && ~edge); flat = flat./conv2(o,k,shape); end
    
    % Do the actual convolution
    c = conv2(a,k,shape)./flat;
    
    % If requested, replace output values with NaNs corresponding to input.
    if(nanout); c(n) = NaN; end
    
    % If 1D, convert back to the original shape.
    if(is1D && sza(1) == 1); c = c.'; end

end
