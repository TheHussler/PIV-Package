%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Process Planar Fuck %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Rhylan A Huss %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% August 18th 2023 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This code serves as scratch paper for a delerious mind. There are
% multiple plots representing several aspects of the flow geometry. All of
% my Planar Plots for these runs are made from these plots. The code
% forever and always will be sloppy, but it will work. Thank you for
% listening. - The aforementioned delerious mind

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

% Colormaps
mustang = [1 1 1; jet];
temp = redblue(540);
blue = [temp(1:270,:)];
red = [temp(270:end,:)];

%Variables
D = 146.05; % mm 
phi = 32;
nu = 1.5339e-05; %m^2/s

% Default Figure Formats
set(groot,'defaultAxesFontName','Times New Roman')
set(groot,'defaultAxesFontSize',14)
set(groot,'defaultAxesFontWeight','Bold')
set(groot,'defaultFigurePosition',[100 100 800 600])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% File Path Declaration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 32 Degree Model: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FilePath = A:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Planar\32DegreeRounded_01
% ReD = 40K  => 23 => 04.19
% ReD = 95K  => 03 => 09.98
% ReD = 110K => 07 => 11.55
% ReD = 120K => 08 => 12.60
% ReD = 130K => 09 => 13.65
% ReD = 140K => 10 => 14.70
% ReD = 150K => 11 => 15.75
% ReD = 160K => 12 => 16.80
% ReD = 170K => 13 => 17.85
% ReD = 180K => 14 => 18.90
% ReD = 190K => 15 => 19.95
% ReD = 200K => 16 => 21.00
% ReD = Trns => 25 => ~~~~~

% File Path %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% path = {...
%     'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Planar\32DegreeRounded_02\32DegRound_23\SubOverTimeAvg_sL=all\PIV_MP(3x16x16_75%ov)_GPU_01\PostProc';                  % ReD = 40K  => 23 => 04.19
%     % 'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Planar\32DegreeRounded_01\32Deg_TEST_03\SubAvg\TR_PIV_MP(2x16x16_75%ov_ImgCorr)';                          % ReD = 95K  => 03 => 09.98
%     'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Planar\32DegreeRounded_01\32Deg_TEST_07\SubAvg\PIV_MP(3x16x16_75%ov)_GPU_02\PostProc';  % ReD = 110K => 07 => 11.55
%     'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Planar\32DegreeRounded_01\32Deg_TEST_08\SubOverTimeAvg_sL=all\PIV_MP(3x16x16_75%ov)_GPU\PostProc_02';  % ReD = 120K => 08 => 12.60
%     'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Planar\32DegreeRounded_01\32Deg_TEST_09\SubOverTimeAvg_sL=all\PIV_MP(3x16x16_75%ov)_GPU\PostProc';     % ReD = 130K => 09 => 13.65
%     'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Planar\32DegreeRounded_01\32Deg_TEST_10\SubOverTimeAvg_sL=all\PIV_MP(3x16x16_75%ov)_GPU\PostProc';     % ReD = 140K => 10 => 14.70
%     'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Planar\32DegreeRounded_01\32Deg_TEST_11\SubOverTimeAvg_sL=all\PIV_MP(3x16x16_75%ov)_GPU\PostProc';     % ReD = 150K => 11 => 15.75
%     'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Planar\32DegreeRounded_01\32Deg_TEST_12\SubOverTimeAvg_sL=all\PIV_MP(3x16x16_75%ov)_GPU\PostProc';     % ReD = 160K => 12 => 16.80
%     'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Planar\32DegreeRounded_01\32Deg_TEST_13\SubOverTimeAvg_sL=all\PIV_MP(3x16x16_75%ov)_GPU\PostProc';     % ReD = 170K => 13 => 17.85
%     'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Planar\32DegreeRounded_01\32Deg_TEST_14\SubOverTimeAvg_sL=all\PIV_MP(3x16x16_75%ov)_GPU\PostProc_12';     % ReD = 180K => 14 => 18.90
%     'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Planar\32DegreeRounded_01\32Deg_TEST_15\SubOverTimeAvg_sL=all\PIV_MP(3x16x16_75%ov)_GPU\PostProc_01';     % ReD = 190K => 15 => 19.95
%     'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Planar\32DegreeRounded_01\32Deg_TEST_16\SubOverTimeAvg_sL=all\PIV_MP(3x16x16_75%ov)_GPU\PostProc_01'};    % ReD = 200K => 16 => 21.00

path = {...
    'R:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Planar\32DegreeRounded_02\32DegRound_23\SubOverTimeAvg_sL=all\PIV_MP(3x16x16_75%ov)_GPU_01\PostProc';                  % ReD = 40K  => 23 => 04.19
    % 'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Planar\32DegreeRounded_01\32Deg_TEST_03\SubAvg\TR_PIV_MP(2x16x16_75%ov_ImgCorr)';                          % ReD = 95K  => 03 => 09.98
    'R:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Planar\32DegreeRounded_01\32Deg_TEST_07\SubAvg\PIV_MP(3x16x16_75%ov)_GPU_02\PostProc';  % ReD = 110K => 07 => 11.55
    'R:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Planar\32DegreeRounded_01\32Deg_TEST_08\SubOverTimeAvg_sL=all\PIV_MP(3x16x16_75%ov)_GPU\PostProc_02';  % ReD = 120K => 08 => 12.60
    'R:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Planar\32DegreeRounded_01\32Deg_TEST_09\SubOverTimeAvg_sL=all\PIV_MP(3x16x16_75%ov)_GPU\PostProc';     % ReD = 130K => 09 => 13.65
    'R:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Planar\32DegreeRounded_01\32Deg_TEST_10\SubOverTimeAvg_sL=all\PIV_MP(3x16x16_75%ov)_GPU\PostProc';     % ReD = 140K => 10 => 14.70
    'R:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Planar\32DegreeRounded_01\32Deg_TEST_11\SubOverTimeAvg_sL=all\PIV_MP(3x16x16_75%ov)_GPU\PostProc';     % ReD = 150K => 11 => 15.75
    'R:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Planar\32DegreeRounded_01\32Deg_TEST_12\SubOverTimeAvg_sL=all\PIV_MP(3x16x16_75%ov)_GPU\PostProc';     % ReD = 160K => 12 => 16.80
    'R:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Planar\32DegreeRounded_01\32Deg_TEST_13\SubOverTimeAvg_sL=all\PIV_MP(3x16x16_75%ov)_GPU\PostProc';     % ReD = 170K => 13 => 17.85
    'R:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Planar\32DegreeRounded_01\32Deg_TEST_14\SubOverTimeAvg_sL=all\PIV_MP(3x16x16_75%ov)_GPU\PostProc_12';     % ReD = 180K => 14 => 18.90
    'R:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Planar\32DegreeRounded_01\32Deg_TEST_15\SubOverTimeAvg_sL=all\PIV_MP(3x16x16_75%ov)_GPU\PostProc_01';     % ReD = 190K => 15 => 19.95
    'R:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Planar\32DegreeRounded_01\32Deg_TEST_16\SubOverTimeAvg_sL=all\PIV_MP(3x16x16_75%ov)_GPU\PostProc_01'};    % ReD = 200K => 16 => 21.00


% path = {...
%     'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Planar\45DegRound_05\45Deg_Round11\SubOverTimeAvg_sL=all_01\PIV_MP(3x16x16_75%ov)_GPU\PostProc';   % ReD = 120K => 08 => 12.60
%     'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Planar\45DegRound_05\45Deg_Round13\SubOverTimeAvg_sL=all\PIV_MP(3x16x16_75%ov)_GPU\PostProc';      % ReD = 150K => 11 => 15.75
%     'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Planar\45DegRound_05\45Deg_Round14\SubOverTimeAvg_sL=all\PIV_MP(3x16x16_75%ov)_GPU\PostProc';      % ReD = 190K => 15 => 19.95
%     'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Planar\45DegRound_05\45Deg_Round15\SubOverTimeAvg_sL=all\PIV_MP(3x16x16_75%ov)_GPU\PostProc'};     % ReD = 40K  => 23 => 04.19

% Freestream Velocity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Winf = [...
    04.19;
    % 09.98;
    11.55;
    12.60;
    13.65;
    14.70;
    15.75;
    16.80;
    17.85;
    18.90;
    19.95;
    21.00]; % m/s

% Winf = [...
% 12.60;
% 15.75;
% 19.95;
% 04.19]; % m/s

% Reynolds Number %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ReD = Winf*(D/1000)/nu;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Start Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load PIV Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load and Average the data set %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SaveBool = input(sprintf('Save Data? [1/0]\n\n'));

clear Data
for redI = 1%:length(path)
    clear vzi vyi uvz uvy
    files = dir(fullfile(path{redI}, '*.vc7'));

    % Preallocates Arrays %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic; Vec = loadvec(fullfile(path{redI},files(1).name)); 
    TotalFiles = min(1000,length(files));
    vzi = NaN([size(Vec.vx'),TotalFiles]); vyi = NaN([size(Vec.vx'),TotalFiles]); 
    Counter = zeros(size(Vec.vx'));

    % Loop Through .vc7 files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Reading in ReD: %iK\nPath: %s\n\n',round(ReD(redI),-4)/1000,path{redI});
    for I=1:TotalFiles %Loopity Loop that dataset
      
        % Load in Velocity Components and Uncertainties %%%%%%%%%%%%%%%%%%%
        IMX = readimx(fullfile(path{redI},files(I).name)); Frames = IMX.Frames{1,1}; 
        Components = Frames.ComponentNames;
        if I == 1 
            U0i = findComponent(Components,'U0');
            V0i = findComponent(Components,'V0');
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
        vzi(:,:,I)=  ((Frames.Components{U0i,1}.Planes{1,1} * Frames.Components{U0i,1}.Scale.Slope) + Frames.Components{U0i,1}.Scale.Offset)'; 
        vyi(:,:,I)= -((Frames.Components{V0i,1}.Planes{1,1} * Frames.Components{V0i,1}.Scale.Slope) + Frames.Components{V0i,1}.Scale.Offset)'; 

        if I == 1
            % Rebuild Coordinate System %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            z = (Vec.x-Vec.x(1));               y = -(Vec.y-Vec.y(1));
            % z = z - max(abs(z))*(275/1280);     y = y + max(abs(y))*(95/800);
            z = z - max(abs(z))*(390/1280);      y = y + max(abs(y))*(278/800);
            % z = z - max(abs(z))*(412/1280);     y = y + max(abs(y))*((800-549)/800);
            vTest = sqrt(vzi(:,:,I).^2 + vyi(:,:,I).^2);
    
            % Generate Masks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Rounded-Edge Model Vertices
            R = (18.487/60)*D;     %mm
            phj = (0:0.1:phi)';  %deg
            modelVert = [R*cosd(90-phj)-R*tand(phi/2),R*sind(90-phj)-R];
            whiteWash = [min(z)/D,max(y)/D; 0.07,max(y)/D; -0.27,0; modelVert/D; 1/tand(phi),-1; 1/tand(phi),min(y)/D; min(z)/D, min(y)/D];
            modelVert = [min(z)/D,0; modelVert/D; 1/tand(phi),-1; min(z)/D, -1];

            % Generates a geometric mask for plotting
            figure(2); clf; imagesc(z/D,y/D,vTest); set(gca,'ydir','normal');
            ROI = images.roi.Polygon(gca,'Position',modelVert);
            mask = double(~createMask(ROI)); mask(mask==0) = NaN;
            
            % Generates a geometric mask for BL Slices
            figure(2); clf; imagesc(z/D,y/D,vTest); set(gca,'ydir','normal');
            ROI = images.roi.Polygon(gca,'Position',whiteWash);
            Mask = double(~createMask(ROI)); Mask(Mask==0) = NaN; close 2;

            % Local Crop
            Locs = [128, 12; 128, 61; 160, 68; 160, 12; 128, 12];
            [X,Y] = meshgrid(1:320,1:200);
            [in,on] = inpolygon(X(:),Y(:),Locs(:,1),Locs(:,2));
            in = reshape(in,size(X));

        end

        % Define NaN Mask
        M = (double(CV > 0.175).*double(UV<2.75).*double(En)); M(M==0) = NaN;
        vzi(:,:,I) =  M .* vzi(:,:,I); vyi(:,:,I) =  M .* vyi(:,:,I);
    
        % Allowable Vector Range %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % M = ((abs(vzi(:,:,I)) < 1.2*Winf(redI)).*(abs(vyi(:,:,I)) < 1.2*Winf(redI))); M(M==0) = NaN;
        % % vzi(:,:,I) = M.*vzi(:,:,I); vyi(:,:,I) = M.*vyi(:,:,I);

        % Allowable Vector Range %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        M = ((abs(vzi(:,:,I)) < 1.5*Winf(redI)).*(abs(vyi(:,:,I)) < 1.2*Winf(redI))); M(M==0) = NaN;
        vzi(:,:,I) = M.*vzi(:,:,I); vyi(:,:,I) = M.*vyi(:,:,I);

        % % Local Allowable Vector Range %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % M = ~in + in.*(sqrt(vzi(:,:,I).^2 + vyi(:,:,I).^2) > 1.0*Winf(redI)); M(M==0) = NaN;
        % vzi(:,:,I) = M.*vzi(:,:,I); vyi(:,:,I) = M.*vyi(:,:,I);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Gaussian Smoothing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        filtWidth = 5;
        filtSigma = 3;
        imageFilter=fspecial('gaussian',filtWidth,filtSigma);
        vzWorking = nanconv(vzi(:,:,I),imageFilter, 'nanout');
        vyWorking = nanconv(vyi(:,:,I),imageFilter, 'nanout');
        vzWorking(vzWorking==0) = NaN; vyWorking(vyWorking==0) = NaN;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Median Filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % vzWorking = vzi(:,:,I); vyWorking = vyi(:,:,I);
        % vzWorking = fast_mediannan(vzWorking, 3);
        % vyWorking = fast_mediannan(vyWorking, 3);
        % vzWorking(vzWorking==0) = NaN; vyWorking(vyWorking==0) = NaN;
        
        % % Local Allowable Vector Range %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % M = ~in + in.*(sqrt(vzWorking.^2 + vyWorking.^2) > 0.9*Winf(redI)); M(M==0) = NaN;
        % vzWorking = M.*vzWorking; vyWorking = M.*vyWorking;
        % Finalize Output Vectors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        vzi(:,:,I) = vzWorking.*Mask; vyi(:,:,I) = vyWorking.*Mask;

        % figure(69);
        % contourf(sqrt(vzi(:,:,I).^2 + vyi(:,:,I).^2))

        % Display Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        clc; 
        fprintf('Reading in ReD: %iK\nPath: %s\n\n',round(ReD(redI),-4)/1000,path{redI});
        fprintf('Read File: %s\n',files(I).name);
        fprintf('Vector Plot %i/%i\n',I,TotalFiles);
        fprintf('Time Elapsed: %.0f of %.0f minutes\n',toc/60,(TotalFiles/I) * (toc/60));

        % pause;
    end 

    % Compile Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clc; disp('Compiling Data...')

    vzWorking = vzi; vyWorking = vyi;
    vzWorking(isnan(vzWorking)) = 0; vyWorking(isnan(vyWorking)) = 0;
    ActiveVz = sum(~isnan(vzi),3); ActiveVy = sum(~isnan(vyi),3); 
    vzWorking = sum(vzWorking,3)./ActiveVz; vyWorking = sum(vyWorking,3)./ActiveVy;

    % Fill Holes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vzWorking = inpaint_nans(vzWorking,2);
    vyWorking = inpaint_nans(vyWorking,2);

    % Store and Clear Data Accordingly %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Store Residuals
    Data.ReD = ReD(redI);
    Data.Winf = Winf(redI);
    Data.Coordinates.z = z; Data.Coordinates.y = y;
    Data.Mask = Mask;

    % Store Vector Quantities
    Data.Vec.Inst.Vz = vzi;
    Data.Vec.Inst.Vy = vyi;
    Data.Vec.Avg.Vz = Mask.*vzWorking;
    Data.Vec.Avg.Vy = Mask.*vyWorking;
    Data.Vec.Avg.Vmag = sqrt(Data.Vec.Avg.Vz.^2+Data.Vec.Avg.Vy.^2);

    % Store Scalar Quantities
    uvz = vzi - Data.Vec.Avg.Vz; uvy = vyi - Data.Vec.Avg.Vy;
    Data.Vec.StdX = Mask.*sqrt((1/TotalFiles) * sum(uvz - vzWorking,3,'omitnan'));
    Data.Vec.StdX = Mask.*sqrt((1/TotalFiles) * sum(uvy - vyWorking,3,'omitnan'));
    Data.Vec.RMSX = Mask.*sqrt((1/TotalFiles) * sum(vzi,3,'omitnan'));
    Data.Vec.RMSY = Mask.*sqrt((1/TotalFiles) * sum(vyi,3,'omitnan'));
    Data.Vec.I = sqrt(mean(uvz.^2,3,'omitnan') + mean(uvy.^2,3,'omitnan'))/Winf(redI);

    % figure(10)
    % contourf(Data(redI).Coordinates.z/D,Data(redI).Coordinates.y/D,Data(redI).Vec.Avg.Vmag/Data(redI).Winf,[-0.125:0.0625:1.5], 'LineStyle','none');
    % daspect([1,1,1]); fprintf('\nSanity Check...\n');
    % colormap jet

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Plot and Play %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    sigma = 1.2;
    Vz = Data.Vec.Avg.Vz; Vy = Data.Vec.Avg.Vy;
    Z = Data.Coordinates.z; Y = Data.Coordinates.y;
    [Z,Y] = meshgrid(Z,Y);

    % Calculate Vorticity and Divergence %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dx = Z(1,2)-Z(1,1); dy = Y(2,1)-Y(1,1);
    
    dwdy = zeros(size(Vz)); dvdx = zeros(size(Vz));
    dwdy(2:end-1,:) = (Vz(3:end,:)-Vz(1:end-2,:))./(2*dy) * D/Winf(redI);
    dvdx(:,2:end-1) = (Vy(:,3:end)-Vy(:,1:end-2))./(2*dx) * D/Winf(redI);
    
    vortx = dwdy - dvdx;

    skip = 5;
    plotVorticity(Z/D,Y/D,vortx,0.5,Mask,modelVert,whiteWash,blue)
    hold on; quiver(Z(1:skip:end,1:skip:end)/D,Y(1:skip:end,1:skip:end)/D,...
        Vz(1:skip:end,1:skip:end),Vy(1:skip:end,1:skip:end));
    plotVelocity(Z/D,Y/D,Data.Vec.Avg.Vmag/Winf(redI),1.2,Mask,modelVert,whiteWash,jet)
    hold on; quiver(Z(1:skip:end,1:skip:end)/D,Y(1:skip:end,1:skip:end)/D,...
        Vz(1:skip:end,1:skip:end),Vy(1:skip:end,1:skip:end));

    % Sanity Check %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Sanity Check...\n\n'); % pause;
    SavePath = sprintf('PlanarPIVData_%i_ReD%iK_nImg_%i_R03.mat',phi,round(Data.ReD/1000,-1),TotalFiles);
    if SaveBool; save(SavePath,"Data"); fprintf('Saved...\n\n');
    else; fprintf('No Save!\n\n'); 
    end
    
end
disp('done...')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot Vorticity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = plotVorticity(z,y,vortx,sigma,mask,modelVert,whiteWash,coloringBook)
    
    vortx(isnan(vortx)) = 0;

    figure(101); clf; hold on; 
    set(gca,'Visible','off'); set(gca,'XTick',[]); set(gca,'YTick',[]);
    ax1 = axes;
    contourf(z,y, mask.*imgaussfilt(vortx,sigma),[-10,0:2.5:30], 'LineStyle','none');
    % xlim([-0.2 1.75]); ylim([-1.10 0.25]); daspect([1 1 1]);
    xlim([-0.2 1.25]); ylim([-1.10 0.25]); daspect([1 1 1]);
    xlabel('z/D'); ylabel('y/D');
    clim([0 30])
    
    hold on; fill(whiteWash(:,1),whiteWash(:,2),[1 1 1],'EdgeColor','none');
    
    view(2); ax2 = axes;
    patch(modelVert(:,1),modelVert(:,2),modelVert(:,2),'LineWidth',2,'EdgeColor',[0.5,0.5,0.5])
    % xlim([-0.2 1.75]); ylim([-1.10 0.25]); daspect([1 1 1]);
    xlim([-0.2 1.25]); ylim([-1.10 0.25]); daspect([1 1 1]);
    ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; ax2.ZTick = [];
    ax2.XTickLabel = {}; ax2.YTickLabel = {};
    
    linkaxes([ax1,ax2])
    ax1.InnerPosition = [0.06,0.11,0.775,0.815];
    ax2.InnerPosition = [0.06,0.11,0.775,0.815];
    
    %%Give each one its own colormap
    colormap(ax2,'gray')
    colormap(ax1,flip(coloringBook))
    set(gca,'ydir','normal')
    c = colorbar(ax1,'Position',[.82 .11 .0675 .815]);
    c.Label.String = '\omega_x(D/W_\infty)';
    
    set(gcf, 'Position', [2.6679e+03 275.8571 987.4286 600]);

end

% Plot Vorticity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = plotVelocity(z,y,vmag,sigma,mask,modelVert,whiteWash,coloringBook)
    
    vmag(isnan(vmag)) = 0;

    figure(102); clf; hold on; 
    set(gca,'Visible','off'); set(gca,'XTick',[]); set(gca,'YTick',[]);
    ax1 = axes;
    contourf(z,y, mask.*vmag,[-10,0:0.05:1.25,5], 'LineStyle','none');
    % xlim([-0.2 1.75]); ylim([-1.10 0.25]); daspect([1 1 1]);
    xlim([-0.2 1.25]); ylim([-1.10 0.25]); daspect([1 1 1]);
    xlabel('z/D'); ylabel('y/D');
    clim([0 1.25])
    
    hold on; fill(whiteWash(:,1),whiteWash(:,2),[1 1 1],'EdgeColor','none');
    
    view(2); ax2 = axes;
    patch(modelVert(:,1),modelVert(:,2),modelVert(:,2),'LineWidth',2,'EdgeColor',[0.5,0.5,0.5])
    % xlim([-0.2 1.75]); ylim([-1.10 0.25]); daspect([1 1 1]);
    xlim([-0.2 1.25]); ylim([-1.10 0.25]); daspect([1 1 1]);
    ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; ax2.ZTick = [];
    ax2.XTickLabel = {}; ax2.YTickLabel = {};
    
    linkaxes([ax1,ax2])
    ax1.InnerPosition = [0.06,0.11,0.775,0.815];
    ax2.InnerPosition = [0.06,0.11,0.775,0.815];
    
    %%Give each one its own colormap
    colormap(ax2,'gray')
    colormap(ax1,flip(coloringBook))
    set(gca,'ydir','normal')
    c = colorbar(ax1,'Position',[.82 .11 .0675 .815]);
    c.Label.String = '|W|/W_\infty';
    
    set(gcf, 'Position', [2.6679e+03 275.8571 987.4286 600]);

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

% First Median Filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function M = mediannan(A, sz)

    if nargin<2
        sz = 5;
    end
    if length(sz)==1
        sz = [sz sz];
    end
    if any(mod(sz,2)==0)
        error('kernel size SZ must be odd)')
    end
    margin=(sz-1)/2;
    AA = nan(size(A)+2*margin);
    AA(1+margin(1):end-margin(1),1+margin(2):end-margin(2))=A;
    [iB jB]=ndgrid(1:sz(1),1:sz(2));
    is=sub2ind(size(AA),iB,jB);
    [iA jA]=ndgrid(1:size(A,1),1:size(A,2));
    iA=sub2ind(size(AA),iA,jA);
    idx=bsxfun(@plus,iA(:).',is(:)-1);
    
    B = sort(AA(idx),1);
    j=any(isnan(B),1);
    last = zeros(1,size(B,2))+size(B,1);
    [trash last(j)]=max(isnan(B(:,j)),[],1);
    last(j)=last(j)-1;
    
    M = nan(1,size(B,2));
    valid = find(last>0);
    mid = (1 + last)/2;
    i1 = floor(mid(valid));
    i2 = ceil(mid(valid));
    i1 = sub2ind(size(B),i1,valid);
    i2 = sub2ind(size(B),i2,valid);
    M(valid) = 0.5*(B(i1) + B(i2));
    M = reshape(M,size(A));

end

% Fast Median Filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = fast_mediannan(A, windowSize)
    % A=inpaint_nans(A,2);
    nanMask = isnan(A);
    A(nanMask) = 0;
    count = medfilt2(~nanMask, [windowSize windowSize], 'symmetric');
    sumFiltered = medfilt2(A, [windowSize windowSize], 'symmetric');
    out = sumFiltered ./ max(count, 1); % Avoid division by 0
    out(nanMask & count == 0) = NaN; % Preserve NaNs where appropriate
    % out=inpaint_nans(out,2);
    % out = regionfill(out, isnan(out));
end

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
    if(is1D)
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

% ROTZ  Create a 3x3 rotation matrix about the Z-axis.
function R = rotz(theta)
% ROTZ  Create a 3x3 rotation matrix about the Z-axis.
%   R = ROTZ(theta) returns the 3x3 rotation matrix for a rotation
%   by theta degrees about the Z-axis.

    % Convert degrees to radians
    thetaRad = deg2rad(theta);

    % Construct rotation matrix
    R = [cos(thetaRad), -sin(thetaRad), 0;
         sin(thetaRad),  cos(thetaRad), 0;
         0,              0,             1];
end




