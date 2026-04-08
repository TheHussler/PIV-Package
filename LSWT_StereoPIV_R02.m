%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Process Stereo Yoyo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Rhylan A Huss %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% September 30th 2023 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This code serves as scratch paper for a delerious mind. There are
% multiple plots representing several aspects of the flow geometry. All of
% my STEREO Plots for these runs are made from these plots. The code
% forever and always will be sloppy, but it will work. Thank you for
% listening. - The aforementioned delerious mind

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

path = {...
    'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Stereo\32DegRoundStereo_03\32Deg_Round_19\SubOverTimeAvg_sL=all_01\StereoPIV_MP(3x16x16_75%ov)_GPU\PostProc';
    'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Stereo\32DegRoundStereo_03\32Deg_Round_18\SubOverTimeAvg_sL=all\StereoPIV_MP(3x16x16_75%ov)_GPU\PostProc';
    'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Stereo\32DegRoundStereo_03\32Deg_Round_17\SubOverTimeAvg_sL=all\StereoPIV_MP(3x16x16_75%ov)_GPU\PostProc';
    'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Stereo\32DegRoundStereo_03\32Deg_Round_16\SubOverTimeAvg_sL=all\StereoPIV_MP(3x16x16_75%ov)_GPU\PostProc';
    'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Stereo\32DegRoundStereo_03\32Deg_Round_15\SubOverTimeAvg_sL=all\StereoPIV_MP(3x16x16_75%ov)_GPU\PostProc';
    'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Stereo\32DegRoundStereo_03\32Deg_Round_14\SubOverTimeAvg_sL=all\StereoPIV_MP(3x16x16_75%ov)_GPU\PostProc'};

% path = {...
    % 'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Stereo\45DegRoundStereo_01\45DegRound_05\SubOverTimeAvg_sL=all\StereoPIV_MP(3x16x16_75%ov)_GPU';
    % 'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Stereo\45DegRoundStereo_01\45DegRound_03\SubOverTimeAvg_sL=all\StereoPIV_MP(3x16x16_75%ov)_GPU';
    % 'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Stereo\45DegRoundStereo_01\45DegRound_02\SubOverTimeMin_sL=all_01\StereoPIV_MP(3x16x16_75%ov)_GPU_01';
    % 'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Stereo\45DegRoundStereo_02\45Deg_Round_10\SubOverTimeMin_sL=all\StereoPIV_MP(2x32x32_75%ov)'};
    % 'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Stereo\45DegRoundStereo_01\45DegRound_04\SubOverTimeAvg_sL=all_01\StereoPIV_MP(3x16x16_75%ov)_GPU'};

% Freestream Velocity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Winf = [...
    04.19;
    12.60;
    15.75;
    16.80;
    17.85;
    19.95]; % m/s

% Winf = [...
%     04.19;
%     12.60;
%     15.75;
%     19.95]; % m/s

% Reynolds Number %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ReD = Winf*(D/1000)/nu;

%% Load PIV Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SaveBool = 0;%input('Save Data? [1/0]\n\n');
% Load and Average the data set %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for redI = 1%4%1:length(path)
    clear vxi vyi uvz uvy
    files = dir(fullfile(path{redI}, '*.vc7'));

    % Preallocates Arrays %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic; Vec = loadvec(fullfile(path{redI},files(1).name)); 
    TotalFiles = min(10000,length(files));
    vxi = NaN([size(Vec.vx'),TotalFiles]); 
    vyi = NaN([size(Vec.vx'),TotalFiles]); 
    vzi = NaN([size(Vec.vx'),TotalFiles]); 
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

        % figure(69);
        % contourf(x,y,sqrt(vxi(:,:,I).^2 + vyi(:,:,I).^2 + vzi(:,:,I).^2))
        % clim([0 Winf(redI)*1.5])

        % Display Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        clc; 
        fprintf('Reading in ReD: %iK\nPath: %s\n\n',round(ReD(redI),-4)/1000,path{redI});
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

    % Calculate Standard Deviation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % vxWorking = vxi; vyWorking = vyi; vzWorking = vzi;
    % vxWorking(isnan(vxWorking)) = 0; vyWorking(isnan(vyWorking)) = 0; vzWorking(isnan(vzWorking)) = 0;
    % ActiveVx = sum(~isnan(vxi),3); ActiveVy = sum(~isnan(vyi),3); ActiveVz = sum(~isnan(vzi),3); 

    % stdX = real(sqrt(sum((vxWorking - repmat(vxAvg,[1,1,length(vxWorking(1,1,:))]) ).^2,3)./(ActiveVx-1))); 
    % stdY = real(sqrt(sum((vyWorking - repmat(vyAvg,[1,1,length(vyWorking(1,1,:))]) ).^2,3)./(ActiveVy-1))); 
    % stdZ = real(sqrt(sum((vzWorking - repmat(vzAvg,[1,1,length(vzWorking(1,1,:))]) ).^2,3)./(ActiveVz-1))); 

    % stdMag = sqrt(stdX.^2+stdY.^2+stdZ.^2);


    %% Calculate Vorticity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [X,Y] = meshgrid(x,y);
    clipMask = (X/D<0.6).*(X/D>-0.6).*(Y/D>-1.35).*(Y/D<-0.66);
    [vortz] = VortZ(X/D,Y/D,vxAvg/Winf(redI),vyAvg/Winf(redI));
    vortz = clipMask.*vortz;

    %% Sanity Check %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure(69); clf; hold on;
    contourf(x/D,y/D,vortz,[-10:0.125:10],'edgecolor','none')
    clim([-8 8]); colormap redblue
    xlim([-0.6 0.601]); ylim([-1.35 -0.65]); daspect([1 1 1]);
    xlabel('x/D'); ylabel('y/D');
    PlotCircle(0,-0.5,0.5); plot(0,0,'k+','MarkerSize',15);

    Slice = streamslice(x/D,y/D,vxAvg/Winf(redI),vyAvg/Winf(redI),3,'noarrows');
    set(Slice,'color',[0.4 0.4 0.4])
    set(Slice,'LineWidth',1)

    set(gcf,'Position',[569 339.2857 681.1429 668.5714]);
    set(gca,'Box','on');

    c = colorbar;
    c.Location = 'NorthOutside';
    c.Label.String = '\omega_z(D/W_\infty)';

    %%
    figure(70); clf; hold on;
    contourf(x/D,y/D,clipMask.*vMag/Winf(redI),[0,0.5:0.005:1.1,2],'edgecolor','none')
    clim([0.5 1.1]); colormap turbo
    xlim([-0.6 0.601]); ylim([-1.35 -0.66]); daspect([1 1 1]);
    xlabel('x/D'); ylabel('y/D');
    PlotCircle(0,-0.5,0.5); plot(0,0,'k+','MarkerSize',15);

    set(gcf,'Position',[1301 613 681.1429 396.5714]);
    set(gca,'Box','on');

    figure(71); clf; hold on;
    contourf(x/D,y/D,clipMask.*stdMag/Winf(redI),[0:0.005:0.4,5],'edgecolor','none')
    clim([0 0.4]); colormap turbo
    xlim([-0.6 0.601]); ylim([-1.35 -0.66]); daspect([1 1 1]);
    xlabel('x/D'); ylabel('y/D');
    PlotCircle(0,-0.5,0.5); plot(0,0,'k+','MarkerSize',15);

    set(gcf,'Position',[1301 613 681.1429 396.5714]);
    set(gca,'Box','on');

    % fprintf('Sanity Check...\n\n'); pause;

    %% Find Vortices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % The Left Vortex Now
    figure(69); clf; hold on;
    contourf(vortz,[-30:2:30],'edgecolor','none')
    clim([-30 30]); colormap redblue
    daspect([1 1 1]);
    
    title('Select POSITIVE Vortex');
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
    roiRight = drawrectangle;
    roiRight = round(roiRight.Position);
    
    jIndex = roiRight(2):roiRight(2)+roiRight(4);
    iIndex = roiRight(1):roiRight(1)+roiRight(3);
    
    % Find Right Vortex Center
    [vortexRight,G1] = IDVortex_Gamma1(-x(iIndex),-y(jIndex),...
        vxAvg(jIndex,iIndex),vyAvg(jIndex,iIndex));
    G1Right = max(G1,[],'all'); vortexRight = -vortexRight;
    
    % figure(69); clf; hold on;
    % contourf(x/D,y/D,vortz,[-30:2:30],'edgecolor','none')
    % n = 5; quiver(x(1:n:end)/D,y(1:n:end)/D,vxAvg(1:n:end,1:n:end),vyAvg(1:n:end,1:n:end))
    % clim([-30 30]); colormap redblue
    % xlim([-0.6 0.6]); ylim([-1 -0.35]); daspect([1 1 1]);
    % xlabel('x/D'); ylabel('y/D');

    plot(vortexLeft(1)/D,vortexLeft(2)/D,'ko')
    plot(vortexRight(1)/D,vortexRight(2)/D,'ko')
    
    % close all;
    
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
    CirculationLeft(redI) = max(abs(Circulation(:,2)));
    
    clc;
    fprintf('Left Vortex:\nR0: %f\nCirculation: %f\n\n',R0Left,CirculationLeft(redI));
    
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
    CirculationRight(redI) = max(abs(Circulation(:,2)));
    
    fprintf('Right Vortex:\nR0: %f\nCirculation: %f\n\n',R0Right,CirculationRight(redI));

    %% Consolidate Circulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Gamma(redI) = (abs(CirculationLeft(redI)) + abs(CirculationRight(redI)))/2;
    fprintf('Total Circulation: %f\n\n',Gamma(redI));

    skip = 10;
    quiver(x(1:skip:end)/D,y(1:skip:end)/D,vxAvg(1:skip:end,1:skip:end),vyAvg(1:skip:end,1:skip:end))

    figure(69); clf; hold on;
    contourf(x/D,y/D,vortz,[-20:0.5:20],'edgecolor','none')
    clim([-20 20]); colormap redblue
    % xlim([-0.6 0.601]); ylim([-1.15 0.05]); daspect([1 1 1]);
    xlim([-0.6 0.601]); ylim([-1.35 -0.66]); daspect([1 1 1]);
    xlabel('x/D'); ylabel('y/D');
    PlotCircle(0,-0.5,0.5); plot(0,0,'k+','MarkerSize',15);

    plot(vortexLeft(1)/D,vortexLeft(2)/D,'ko','MarkerFaceColor','k')
    plot(vortexRight(1)/D,vortexRight(2)/D,'ko','MarkerFaceColor','k')

    % set(gcf,'Position',[569 339.2857 681.1429 668.5714]);
    set(gcf,'Position',[1301 613 681.1429 396.5714]);
    set(gca,'Box','on');

    figure(70); clf; hold on;
    contourf(x/D,y/D,clipMask.*vMag/Winf(redI),[0,0.5:0.005:1.1,2],'edgecolor','none')
    clim([0.5 1.1]); colormap turbo
    xlim([-0.6 0.601]); ylim([-1.35 -0.66]); daspect([1 1 1]);
    xlabel('x/D'); ylabel('y/D');
    PlotCircle(0,-0.5,0.5); plot(0,0,'k+','MarkerSize',15);

    plot(vortexLeft(1)/D,vortexLeft(2)/D,'ko','MarkerFaceColor','k')
    plot(vortexRight(1)/D,vortexRight(2)/D,'ko','MarkerFaceColor','k')

    set(gcf,'Position',[1301 613 681.1429 396.5714]);
    set(gca,'Box','on');

    figure(71); clf; hold on;
    contourf(x/D,y/D,clipMask.*stdMag/Winf(redI),[0:0.005:0.4,5],'edgecolor','none')
    clim([0 0.4]); colormap turbo
    xlim([-0.6 0.601]); ylim([-1.35 -0.66]); daspect([1 1 1]);
    xlabel('x/D'); ylabel('y/D');
    PlotCircle(0,-0.5,0.5); plot(0,0,'k+','MarkerSize',15);

    plot(vortexLeft(1)/D,vortexLeft(2)/D,'ko','MarkerFaceColor','k')
    plot(vortexRight(1)/D,vortexRight(2)/D,'ko','MarkerFaceColor','k')

    set(gcf,'Position',[1301 613 681.1429 396.5714]);
    set(gca,'Box','on');

    %% Sanity Check! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pause;

    %% Save Bool %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if SaveBool
        saveStr = sprintf('Stereo_%ideg_ReD%iK_VortZ.fig',phi,round(ReD(redI)/1000,-1));
        figure(69); saveas(gcf,saveStr);
    
        saveStr = sprintf('Stereo_%ideg_ReD%iK_VMag.fig',phi,round(ReD(redI)/1000,-1));
        figure(70); saveas(gcf,saveStr);
    
        saveStr = sprintf('Stereo_%ideg_ReD%iK_Std.fig',phi,round(ReD(redI)/1000,-1));
        figure(71); saveas(gcf,saveStr);
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot Circle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = PlotCircle(xo,yo,R)
    phi = 0:0.5:365;
    x = R*cosd(phi) + xo;
    y = R*sind(phi) + yo;

    hold on; plot(x,y,'k--');
end

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

function h = viscirclesMatlabDum(varargin)

varargin = matlab.images.internal.stringToChar(varargin);

[ax, centers, radii, numCircles, options] = parseInputs(varargin{:});

if isempty(centers)
    h = [];
    return;
end

isHoldOn = ishold(ax);
hold(ax,'on');
cObj = onCleanup(@()preserveHold(isHoldOn,ax)); % Preserve original hold state

thetaResolution = 2; 
theta=(0:thetaResolution:360)'*pi/180;

x = radii' .* cos(theta);
x = x + (centers(:,1))';
x = cat(1,x,nan(1,numCircles));
x = x(:);

y = radii' .* sin(theta);
y = y + (centers(:,2))';
y = cat(1,y,nan(1,numCircles));
y = y(:);

% Create hggroup object that will contain the two circles as children
h = hggroup('Parent', ax);

if options.EnhanceVisibility
    % Draw the thicker background white circle
        
    thickEdgeColor = 'w';    
    thickLineWidth = options.LineWidth + 1;
    if (strcmpi(options.LineStyle,'none'))
        thickLineStyle = 'none';
    else
        thickLineStyle = '-';
    end

    line(x,y,'Parent',h, ...
        'Color',thickEdgeColor, ...
        'LineWidth',thickLineWidth, ...
        'LineStyle',thickLineStyle);
    
end

% Draw the thinner foreground colored circle
line(x,y,'Parent',h, ...
   'Color',options.Color, ...
   'LineWidth',options.LineWidth, ...
   'LineStyle',options.LineStyle);

end

% -------------------------------------------------------------------------

function [ax, centers, radii, numCircles, options] = parseInputs(varargin)

narginchk(2, 11);

needNewAxes = 0;

first_string = min(find(cellfun(@ischar, varargin), 1, 'first'));
if isempty(first_string)
    first_string = length(varargin) + 1;
end

if first_string == 3
    % viscircles(centers, radii)    
    needNewAxes = 1;   
    centers = varargin{1};
    radii = varargin{2};
    
elseif first_string == 4
    % viscircles(ax, centers, radii)
    ax = varargin{1};
    ax = validateAxes(ax);
    
    centers = varargin{2};
    radii = varargin{3};
    
else
    error(message('images:validate:invalidSyntax'))
end

% Handle remaining name-value pair parsing
name_value_pairs = varargin(first_string:end);

num_pairs = numel(name_value_pairs);
if (rem(num_pairs, 2) ~= 0)
    error(message('images:validate:missingParameterValue'));
end

% Do not change the order of argument names listed below
args_names = {'Color','LineWidth','LineStyle','EnhanceVisibility'};
arg_default_values = {'red', 2, '-', true};

% Set default parameter values
for i = 1: numel(args_names)
    options.(args_names{i}) = arg_default_values{i};
end

% Support for older arguments - do not change the order of argument names listed below
args_names = cat(2,args_names, {'EdgeColor', 'DrawBackgroundCircle'});

for i = 1:2:num_pairs
    arg = name_value_pairs{i};
    if ischar(arg)        
        idx = find(strncmpi(arg, args_names, numel(arg)));
        if isempty(idx)
            error(message('images:validate:unknownInputString', arg))
        elseif numel(idx) > 1
            error(message('images:validate:ambiguousInputString', arg))
        elseif numel(idx) == 1
            if(idx == 5) % If 'EdgeColor' is specified
                idx = 1; % Map to 'Color' 
            elseif(idx == 6) % If 'DrawBackgroundCircle' is specified
                idx = 4; % Map to 'EnhanceVisibility'
            end
            options.(args_names{idx}) = name_value_pairs{i+1};
        end    
    else
        error(message('images:validate:mustBeString')); 
    end
end

% Validate parameter values. Let LINE do the validation for EdgeColor,
% LineStyle and LineWidth.
[centers, radii] = validateCentersAndRadii(centers, radii, first_string); 
numCircles = size(centers,1);

options.EnhanceVisibility = validateEnhanceVisibility( ...
    options.EnhanceVisibility);

% If required, create new axes after parsing
if(needNewAxes)    
    ax = gca;
end

end

% -------------------------------------------------------------------------

function preserveHold(wasHoldOn,ax)
% Function for preserving hold behavior on exit
if ~wasHoldOn
    hold(ax,'off');
end

end

function ax = validateAxes(ax)

if ~ishghandle(ax)
    error(message('images:validate:invalidAxes','AX'))
end

objType = get(ax,'type');
if ~strcmp(objType,'axes')
    error(message('images:validate:invalidAxes','AX'))
end

end

function [centers, radii] = validateCentersAndRadii(centers, radii, ...
                                                    first_string)

if(~isempty(centers))
    validateattributes(centers,{'numeric'},{'nonsparse','real', ...
        'ncols',2}, mfilename,'centers',first_string-2);
    validateattributes(radii,{'numeric'},{'nonsparse','real','nonnegative', ...
        'vector'}, mfilename,'radii',first_string-1);
    
    if(~isscalar(radii))
        if(size(centers,1) ~= length(radii))
            error(message('images:validate:unequalNumberOfRows','CENTERS','RADII'))
        end
    end
    centers = double(centers);
    radii   = double(radii(:)); % Convert to a column vector
end

end

function doEnhanceVisibility = validateEnhanceVisibility(doEnhanceVisibility)

if ~(islogical(doEnhanceVisibility) || isnumeric(doEnhanceVisibility)) ...
        || ~isscalar(doEnhanceVisibility)
    error(message('images:validate:invalidLogicalParam', ...
        'EnhanceVisibility', 'VISCIRCLES', 'EnhanceVisibility'))
end

doEnhanceVisibility = logical(doEnhanceVisibility);

end
