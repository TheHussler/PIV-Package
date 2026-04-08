%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Process Stereo Unsteady %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
Garnet = [0.4706 0.1843 0.2510];
Gold = [0.8078 0.7216 0.5333];

% Default Figure Formats
set(groot,'defaultAxesFontName','Times New Roman')
set(groot,'defaultAxesFontSize',14)
set(groot,'defaultAxesFontWeight','Bold')
set(groot,'defaultFigurePosition',[100 100 700 600])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Start Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% path = {...
%     'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Stereo\32DegRoundStereo_03\32Deg_Round_19\SubOverTimeAvg_sL=all_01\StereoPIV_MP(3x16x16_75%ov)_GPU\PostProc';
%     'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Stereo\32DegRoundStereo_03\32Deg_Round_18\SubOverTimeAvg_sL=all_01\StereoPIV_MP(3x16x16_75%ov)_GPU\PostProc';
%     'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Stereo\32DegRoundStereo_03\32Deg_Round_17\SubOverTimeAvg_sL=all_03\StereoPIV_MP(3x16x16_75%ov)_GPU\PostProc';
%     'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Stereo\32DegRoundStereo_03\32Deg_Round_16\SubOverTimeAvg_sL=all\StereoPIV_MP(3x16x16_75%ov)_GPU\PostProc';
%     'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Stereo\32DegRoundStereo_03\32Deg_Round_14\SubOverTimeAvg_sL=all_01\StereoPIV_MP(3x16x16_75%ov)_GPU\PostProc'};

path = {...
    'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Stereo\45DegRoundStereo_01\45DegRound_05\SubOverTimeAvg_sL=all\StereoPIV_MP(3x16x16_75%ov)_GPU';
    'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Stereo\45DegRoundStereo_01\45DegRound_03\SubOverTimeAvg_sL=all\StereoPIV_MP(3x16x16_75%ov)_GPU';
    'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Stereo\45DegRoundStereo_01\45DegRound_02\SubOverTimeMin_sL=all_01\StereoPIV_MP(3x16x16_75%ov)_GPU_01';
    'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Stereo\45DegRoundStereo_01\45DegRound_04\SubOverTimeAvg_sL=all_01\StereoPIV_MP(3x16x16_75%ov)_GPU'};


% path = {...
%     'R:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Stereo\32DegRoundStereo_03\32Deg_Round_19\SubOverTimeAvg_sL=all_01\StereoPIV_MP(3x16x16_75%ov)_GPU\PostProc';
%     'R:\01 - Lofted Cylinder\03 - LSWT HSPIV\02_LSWT PIV\Stereo\32DegRoundStereo_03\32Deg_Round_14\SubOverTimeAvg_sL=all_01\StereoPIV_MP(3x16x16_75%ov)_GPU\PostProc'};

% Freestream Velocity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Winf = [...
%     04.19;
%     12.60;
%     15.75;
%     16.80;
%     19.95]; % m/s

Winf = [...
    04.19;
    12.60;
    15.75;
    19.95]; % m/s

% Reynolds Number %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ReD = Winf*(D/1000)/nu;

%% Load PIV Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SaveBool = input('Save Data? [1/0]\n\n');
VortBool = 0;%input('Track Vortex Movement? [1/0]\n\n');
CompBool = 0;%input('Compress Data for POD/SPOD? [1/0]\n\n');
% Load and Average the data set %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for redI = 1%4:length(path)
    clear vxi vyi uvz uvy
    files = dir(fullfile(path{redI}, '*.vc7'));

    % Preallocates Arrays %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Vec = loadvec(fullfile(path{redI},files(1).name)); 
    TotalFiles = length(files);
    vortexLeft = NaN(TotalFiles,2); vortexRight = NaN(TotalFiles,2);
    % vxi = NaN([size(Vec.vx'),TotalFiles]); 
    % vyi = NaN([size(Vec.vx'),TotalFiles]); 
    % vzi = NaN([size(Vec.vx'),TotalFiles]); 
    Counter = zeros(size(Vec.vx'));

    % figure(72); clf; hold on;
    % xlim([-0.6 0.601]); ylim([-1.35 0.05]); daspect([1 1 1]);
    % xlabel('x/D'); ylabel('y/D');
    % PlotCircle(0,-0.5,0.5); plot(0,0,'k+','MarkerSize',15);
    % set(gcf,'Position',[646.7143  113.5714  681.1429  668.5714]);
    % set(gca,'Box','on');

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
        % vxi(:,:,I)= -((Frames.Components{U0i,1}.Planes{1,1} * Frames.Components{U0i,1}.Scale.Slope) + Frames.Components{U0i,1}.Scale.Offset)'; 
        % vyi(:,:,I)= -((Frames.Components{V0i,1}.Planes{1,1} * Frames.Components{V0i,1}.Scale.Slope) + Frames.Components{V0i,1}.Scale.Offset)'; 
        % vzi(:,:,I)= -((Frames.Components{W0i,1}.Planes{1,1} * Frames.Components{W0i,1}.Scale.Slope) + Frames.Components{W0i,1}.Scale.Offset)'; 
        vxi= double(-((Frames.Components{U0i,1}.Planes{1,1} * Frames.Components{U0i,1}.Scale.Slope) + Frames.Components{U0i,1}.Scale.Offset)'); 
        vyi= double(-((Frames.Components{V0i,1}.Planes{1,1} * Frames.Components{V0i,1}.Scale.Slope) + Frames.Components{V0i,1}.Scale.Offset)'); 
        vzi= double(-((Frames.Components{W0i,1}.Planes{1,1} * Frames.Components{W0i,1}.Scale.Slope) + Frames.Components{W0i,1}.Scale.Offset)'); 

        if I == 1
            % Rebuild Coordinate System %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            x = -(Vec.x); y = -(Vec.y);
            x=(x - 0); y=y - (99.225 + 51.5); % (76.2 + 51.5)
            [X,Y] = meshgrid(x,y);
        end

        % Define NaN Mask
        M = (double(En)); M(M==0) = NaN; %double(CV > 0.2).*
        vxi =  M .* vxi; 
        vyi =  M .* vyi;
        vzi =  M .* vzi;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Gaussian Smoothing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % filtWidth = 5;
        % filtSigma = 1.5;

        filtWidth = 15;
        filtSigma = 1.5;
        imageFilter=fspecial('gaussian',filtWidth,filtSigma);
        vxWorking = nanconv(vxi,imageFilter, 'nanout');
        vyWorking = nanconv(vyi,imageFilter, 'nanout');
        vzWorking = nanconv(vzi,imageFilter, 'nanout');
        vxWorking(vxWorking==0) = NaN; vyWorking(vyWorking==0) = NaN; vzWorking(vzWorking==0) = NaN;

        % Finalize Output Vectors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        vxi = vxWorking; vyi = vyWorking; vzi = vzWorking;
        
        % Calculate Vorticity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        clipMask = (X/D<0.6).*(X/D>-0.6).*(Y/D>-1.35).*(Y/D<-0.66);
        [vortzi] = VortZ(X/D,Y/D,vxWorking/Winf(redI),vyWorking/Winf(redI));

        %% Sanity Check %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % figure(69); clf; hold on;
        % contourf(x/D,y/D,vortzi,[-30:2:30],'edgecolor','none')
        % clim([-20 20]); colormap redblue
        % xlim([-0.6 0.601]); ylim([-1.35 0.05]); daspect([1 1 1]);
        % xlabel('x/D'); ylabel('y/D');
        % PlotCircle(0,-0.5,0.5); plot(0,0,'k+','MarkerSize',15);
        % % l = streamslice(X/D,Y/D,vxi,vyi,5,'noarrows');
        % % set(l,'LineWidth',1.2)
        % % set(l,'Color','k');
        % set(gcf,'Position',[646.7143  113.5714  681.1429  668.5714]);
        % set(gca,'Box','on');

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Calculate Vortex Center %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if VortBool
            rim = 0.15*D;
    
            vortCut = 1.5;
            vxWorking(isnan(vxWorking)) = 0;
            vyWorking(isnan(vyWorking)) = 0;
            SecondaryMask = (X/D>-0.45).*(X/D<0.45).*(Y/D>-1.3).*(Y/D<-0.75);
            Neg = (X/D<0.1) .* ((vortzi<0).*(abs(vortzi)>vortCut)).*SecondaryMask;
            Pos = (X/D>-0.1) .* ((vortzi>0).*(abs(vortzi)>vortCut)).*SecondaryMask; 
            vxWorking = vxWorking.*SecondaryMask; vyWorking = vyWorking.*SecondaryMask;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Positive Vortex %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if I ==1
                figure(1); clf; hold on;
                contourf(vortzi,[-30:2:30],'edgecolor','none')
                l = streamslice(vxi,vyi,5,'noarrows');
                set(l,'LineWidth',1.2)
                set(l,'Color','k');
                clim([-30 30]); colormap redblue
                daspect([1 1 1]);
        
                title('Select POSITIVE Vortex');
                roiLeft=drawrectangle; roiLeft=round(roiLeft.Position); close 1;
                jIndex = roiLeft(2):roiLeft(2)+roiLeft(4);
                iIndex = roiLeft(1):roiLeft(1)+roiLeft(3);
    
                % Find Left Vortex Center
                [vortexLefti,G1] = IDVortex_Gamma1(-x(iIndex),-y(jIndex),...
                    vxWorking(jIndex,iIndex),vyWorking(jIndex,iIndex));
                G1Left = max(G1,[],'all'); vortexLeft(I,:) = -vortexLefti;
    
            else
                try
                    % Automatically Selected Mask
                    tempMaskLeft = zeros(size(vortzi));
                    for R = 1:length(reshape(vortzi,[],1))
                        if sqrt((X(R) - vortexLeft(I-1,1))^2 + (Y(R) - vortexLeft(I-1,2))^2) <= rim
                            tempMaskLeft(R) = 1;
                        end
                    end
                    tempMaskLeft = Pos.*tempMaskLeft;
        
                    % figure()
                    % contourf(tempMaskLeft.*vortzi,[-30:2:30],'edgecolor','none')
                    % clim([-30 30]); colormap redblue
                    % daspect([1 1 1]);
        
                    % Find Right Vortex Center
                    [vortexLefti,G1] = IDVortex_Gamma1(-x,-y,...
                        tempMaskLeft.*vxWorking,tempMaskLeft.*vyWorking);
                    G1Left = max(G1,[],'all'); vortexLeft(I,:) = -vortexLefti;
                catch
                    rim = rim + 0.15*D;
                    % Automatically Selected Mask
                    tempMaskLeft = zeros(size(vortzi));
                    for R = 1:length(reshape(vortzi,[],1))
                        if sqrt((X(R) - vortexLeft(I-1,1))^2 + (Y(R) - vortexLeft(I-1,2))^2) <= rim
                            tempMaskLeft(R) = 1;
                        end
                    end
                    tempMaskLeft = Pos.*tempMaskLeft;
        
                    % figure()
                    % contourf(tempMaskLeft.*vortzi,[-30:2:30],'edgecolor','none')
                    % clim([-30 30]); colormap redblue
                    % daspect([1 1 1]);
        
                    % Find Right Vortex Center
                    [vortexLefti,G1] = IDVortex_Gamma1(-x,-y,...
                        tempMaskLeft.*vxWorking,tempMaskLeft.*vyWorking);
                    G1Left = max(G1,[],'all'); vortexLeft(I,:) = -vortexLefti;
                end
    
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Negative Vortex %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if I == 1
                figure(1); clf; hold on;
                contourf(vortzi,[-30:2:30],'edgecolor','none')
                l = streamslice(vxi,vyi,5,'noarrows');
                set(l,'LineWidth',1.2)
                set(l,'Color','k');
                clim([-30 30]); colormap redblue
                daspect([1 1 1]);
        
                % Manually Select Vortex
                title('Select NEGATIVE Vortex');
                roiRight=drawrectangle; roiRight=round(roiRight.Position); close 1;
                jIndex = roiRight(2):roiRight(2)+roiRight(4);
                iIndex = roiRight(1):roiRight(1)+roiRight(3);
    
                % Find Right Vortex Center
                [vortexRighti,G1] = IDVortex_Gamma1(-x(iIndex),-y(jIndex),...
                    vxWorking(jIndex,iIndex),vyWorking(jIndex,iIndex));
                G1Right = max(G1,[],'all'); vortexRight(I,:) = -vortexRighti;
    
            else
                try
                    % Automatically Selected Mask
                    tempMaskRight = zeros(size(vortzi));
                    for R = 1:length(reshape(vortzi,[],1))
                        if sqrt((X(R) - vortexRight(I-1,1))^2 + (Y(R) - vortexRight(I-1,2))^2) <= rim
                            tempMaskRight(R) = 1;
                        end
                    end
                    tempMaskRight = Neg.*tempMaskRight;
        
                    % figure()
                    % contourf(tempMaskRight.*vortzi,[-30:2:30],'edgecolor','none')
                    % clim([-30 30]); colormap redblue
                    % daspect([1 1 1]);
        
                    % Find Right Vortex Center
                    [vortexRighti,G1] = IDVortex_Gamma1(-x,-y,...
                        tempMaskRight.*vxWorking,tempMaskRight.*vyWorking);
                    G1Right = max(G1,[],'all'); vortexRight(I,:) = -vortexRighti;
                catch
                    rim = rim + 0.15*D;
                    % Automatically Selected Mask
                    tempMaskRight = zeros(size(vortzi));
                    for R = 1:length(reshape(vortzi,[],1))
                        if sqrt((X(R) - vortexRight(I-1,1))^2 + (Y(R) - vortexRight(I-1,2))^2) <= rim
                            tempMaskRight(R) = 1;
                        end
                    end
                    tempMaskRight = Neg.*tempMaskRight;
        
                    % figure()
                    % contourf(tempMaskRight.*vortzi,[-30:2:30],'edgecolor','none')
                    % clim([-30 30]); colormap redblue
                    % daspect([1 1 1]);
        
                    % Find Right Vortex Center
                    [vortexRighti,G1] = IDVortex_Gamma1(-x,-y,...
                        tempMaskRight.*vxWorking,tempMaskRight.*vyWorking);
                    G1Right = max(G1,[],'all'); vortexRight(I,:) = -vortexRighti;
                end    
            end

            % figure(72); hold on;
            % if I == 1
            %     meanLeft = plot(mean(vortexLeft(:,1),'omitmissing')/D,mean(vortexLeft(:,2),'omitmissing')/D,'r+');
            %     meanRight = plot(mean(vortexRight(:,1),'omitmissing')/D,mean(vortexRight(:,2),'omitmissing')/D,'r+');
            %     pltLeft = plot(vortexLefti(1)/D,vortexLefti(2)/D,'r.');
            %     pltRight = plot(vortexRighti(1)/D,vortexRighti(2)/D,'r.');
            % else
            %     pltLeft.Marker = '.'; pltLeft.Color = [0 0 0]; 
            %     pltRight.Marker = '.'; pltRight.Color = [0 0 0];
            %     delete(meanLeft); delete(meanRight);
            %     meanLeft = plot(mean(vortexLeft(:,1),'omitmissing')/D,mean(vortexLeft(:,2),'omitmissing')/D,'r+');
            %     meanRight = plot(mean(vortexRight(:,1),'omitmissing')/D,mean(vortexRight(:,2),'omitmissing')/D,'r+');
            %     pltLeft = plot(vortexLeft(I,1)/D,vortexLeft(I,2)/D,'ro','MarkerFaceColor','r');
            %     pltRight = plot(vortexRight(I,1)/D,vortexRight(I,2)/D,'ro','MarkerFaceColor','r');
            % end
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end

        if CompBool
            if I == 1
                % PreAllocate
                scale = 0.5;
                vxComp = NaN([length(files),size(imresize(vxi,scale))]);
                vyComp = NaN([length(files),size(imresize(vxi,scale))]);
                vzComp = NaN([length(files),size(imresize(vxi,scale))]);
            end
            % Store Compressed Images
            vxComp(I,:,:) = imresize(vxi,scale);
            vyComp(I,:,:) = imresize(vyi,scale);
            vzComp(I,:,:) = imresize(vzi,scale);
        end

        % Display Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        clc;

        if I==1
            tic; 
        else
            fprintf('Reading in ReD: %iK\nPath: %s\n\n',round(ReD(redI),-4)/1000,path{redI});
            fprintf('Read File: %s\n',files(I).name);
            fprintf('Vector Plot %i/%i\n',I,TotalFiles);
            fprintf('Time Elapsed: %.0f of %.0f minutes\n',toc/60,(TotalFiles/I) * (toc/60));
        end
    end 

    %% Plot Spectra %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if VortBool
        %% Calculate and Plot Vort Tracking Spectra %%%%%%%%%%%%%%%%%%%%%%%
        % fs = 5000; df=2; winSize=round(fs/df);
        % fs = 5000; df=4; winSize=round(fs/df);
        fs = 5000; df=7.5; winSize=round(fs/df);
        
    
        [pxLeft,f]=pwelch((vortexLeft(:,1) - mean(vortexLeft(:,1),'omitmissing'))/D, hann(winSize), round(winSize*0.75), winSize, fs);
        [pyLeft,~]=pwelch((vortexLeft(:,2) - mean(vortexLeft(:,2),'omitmissing'))/D, hann(winSize), round(winSize*0.75), winSize, fs);
    
        [pxRight,~]=pwelch((vortexRight(:,1) - mean(vortexRight(:,1),'omitmissing'))/D, hann(winSize), round(winSize*0.75), winSize, fs);
        [pyRight,~]=pwelch((vortexRight(:,2) - mean(vortexRight(:,2),'omitmissing'))/D, hann(winSize), round(winSize*0.75), winSize, fs);
    
        St=f*(D/1000)/Winf(redI);
    
        % figure(400); clf;
        % loglog(f,pxLeft,'-','Color',Garnet,'LineWidth',2); hold on; loglog(f,pyLeft,'--','Color',Garnet,'LineWidth',2); 
        % loglog(f,pxRight,'-','Color',Gold,'LineWidth',2); loglog(f,pyRight,'--','Color',Gold,'LineWidth',2);
        % xlabel('f [Hz]'); ylabel('PSD'); grid on
    
        figure(601); clf;
        loglog(St,pxLeft,'-','Color',Garnet,'LineWidth',2); hold on; loglog(St,pyLeft,'--','Color',Garnet,'LineWidth',2);
        loglog(St,pxRight,'-','Color',Gold,'LineWidth',2); loglog(St,pyRight,'--','Color',Gold,'LineWidth',2);
        % semilogy(St,pxLeft,'-','Color',Garnet,'LineWidth',2); hold on; loglog(St,pyLeft,'--','Color',Garnet,'LineWidth',2);
        % semilogy(St,pxRight,'-','Color',Gold,'LineWidth',2); loglog(St,pyRight,'--','Color',Gold,'LineWidth',2);
        legend('X - Left Vort.', 'Y - Left Vort.', 'X - Right Vort.', 'Y - Right Vort.')
        legend('Location','southwest')
        xlabel('St_D'); ylabel('PSD'); grid on
        xlim([0.1 10])
        set(gcf,'Position',[1301 613 681.1429 396.5714]);
        set(gca,'Box','on');
    
        %% Plot all Locations on a Mean Plot
    
        VortLeftAvg = mean(vortexLeft);
        WorkingLeftX = (vortexLeft(:,1) - VortLeftAvg(1))/D; WorkingLeftY = (vortexLeft(:,2) - VortLeftAvg(2))/D;
    
        % Smooth the shape with a median filter.
        windowWidth = 50; cutoff = 0.07;
    
	    smoothedx = medfilt1(WorkingLeftX, windowWidth);
	    smoothedy = medfilt1(WorkingLeftY, windowWidth);
	    residuals = sqrt((WorkingLeftX-smoothedx) .^ 2 + (WorkingLeftY - smoothedy) .^ 2);
        In = double(residuals < cutoff); In(In==0) = NaN;
    
        WorkingLeftX = D*WorkingLeftX.*In + VortLeftAvg(1);
        WorkingLeftY = D*WorkingLeftY.*In + VortLeftAvg(2);
    
        VortRightAvg = mean(vortexRight);
        WorkingRightX = (vortexRight(:,1)-VortRightAvg(1))/D; WorkingRightY = (vortexRight(:,2)-VortRightAvg(2))/D;
    
	    smoothedx = medfilt1(WorkingRightX, windowWidth);
	    smoothedy = medfilt1(WorkingRightY, windowWidth);
	    residuals = sqrt((WorkingRightX-smoothedx) .^ 2 + (WorkingRightY - smoothedy) .^ 2);
        In = double(residuals < cutoff); In(In==0) = NaN;
    
        WorkingRightX = D*WorkingRightX.*In + VortRightAvg(1);
        WorkingRightY = D*WorkingRightY.*In + VortRightAvg(2);
    
        figure(404); clf; hold on
        PlotCircle(0,-0.5,0.5)
        xlim([-0.6 0.601]); ylim([-1.35 -0.66]); daspect([1 1 1]);
        set(gcf,'Position',[1301 613 681.1429 396.5714]);
        set(gca,'Box','on');
    
        % plot(vortexLeft(:,1),vortexLeft(:,2),'r.')
        plot(WorkingLeftX/D,WorkingLeftY/D,'.','MarkerFaceColor',Garnet,'MarkerEdgeColor',Garnet)
        plot(VortLeftAvg(1)/D,VortLeftAvg(2)/D,'o','markerfacecolor','k')
        % plot(mean(WorkingLeftX,'omitmissing'),mean(WorkingLeftY,'omitmissing'),'bo','markerfacecolor','b')
    
        % plot(vortexRight(:,1),vortexRight(:,2),'r.')
        plot(WorkingRightX/D,WorkingRightY/D,'.','MarkerFaceColor',Gold,'MarkerEdgeColor',Gold)
        plot(VortRightAvg(1)/D,VortRightAvg(2)/D,'o','markerfacecolor','k')
        % plot(mean(WorkingRightX,'omitmissing'),mean(WorkingRightY,'omitmissing'),'bo','markerfacecolor','b')
       
        % Calculate Two Dimensional Standard Deviation %%%%%%%%%%%%%%%%%%%%%%%
        figure(404); hold on;
    
        URightX = WorkingRightX - VortRightAvg(1);
        URightY = WorkingRightY - VortRightAvg(2);
        [STRightX_01,STRightY_01] = std2D(URightX,URightY,1);
        STRightX_01 = STRightX_01 + VortRightAvg(1);
        STRightY_01 = STRightY_01 + VortRightAvg(2);
        [STRightX_02,STRightY_02] = std2D(URightX,URightY,2);
        STRightX_02 = STRightX_02 + VortRightAvg(1);
        STRightY_02 = STRightY_02 + VortRightAvg(2);
        [STRightX_03,STRightY_03] = std2D(URightX,URightY,3);
        STRightX_03 = STRightX_03 + VortRightAvg(1);
        STRightY_03 = STRightY_03 + VortRightAvg(2);
        
        color = [0.1 0.1 0.1];
        plot(STRightX_01'/D,STRightY_01'/D,'-','color',color,'linewidth',2)
        plot(STRightX_02'/D,STRightY_02'/D,'-','color',color,'linewidth',2)
        % plot(STRightX_03'/D,STRightY_03'/D,'-','color',color,'linewidth',2)
    

        ULeftX = WorkingLeftX - VortLeftAvg(1);
        ULeftY = WorkingLeftY - VortLeftAvg(2);
        [STLeftX_01,STLeftY_01] = std2D(ULeftX,ULeftY,1);
        STLeftX_01 = STLeftX_01 + VortLeftAvg(1);
        STLeftY_01 = STLeftY_01 + VortLeftAvg(2);
        [STLeftX_02,STLeftY_02] = std2D(ULeftX,ULeftY,2);
        STLeftX_02 = STLeftX_02 + VortLeftAvg(1);
        STLeftY_02 = STLeftY_02 + VortLeftAvg(2);
        [STLeftX_03,STLeftY_03] = std2D(ULeftX,ULeftY,3);
        STLeftX_03 = STLeftX_03 + VortLeftAvg(1);
        STLeftY_03 = STLeftY_03 + VortLeftAvg(2);
    
        color = [0.1 0.1 0.1];
        plot(STLeftX_01'/D,STLeftY_01'/D,'-','color',color,'linewidth',2)
        plot(STLeftX_02'/D,STLeftY_02'/D,'-','color',color,'linewidth',2)
        % plot(STLeftX_03'/D,STLeftY_03'/D,'-','color',color,'linewidth',2)
        xlabel('x/D'); ylabel('y/D');

        save(sprintf('PSD_ReD_%iK_%ideg.mat',round(ReD(redI)/1000,-1),phi))
        figure(404); saveas(gcf,sprintf('VortTrack_%iK_%ideg.fig',round(ReD(redI)/1000,-1),phi));
        % figure(601); saveas(gcf,sprintf('VortTrack_%iK_Spectra_%ideg.fig',round(ReD(redI)/1000,-1),phi));

    end
%% POD Calculations
    if CompBool
        vMagComp= sqrt(vxComp.^2 + vyComp.^2+ vzComp.^2);
        vMagComp(isnan(vMagComp)) = 0;
        [U_POD, S_POD, V_POD] = pod(vMagComp);
        ModeEnergies=S_POD.^2;
%%
        ModeEnergyFraction=ModeEnergies/sum(ModeEnergies);
        figure('Color','w','Position',[146 620 403 357]);
        bar(1:length(ModeEnergies),ModeEnergyFraction,'k');
        title('Mode Energies');
        %Note only the first two modes have significant energy. The noise is spread out through all modes.
        
        %% Plots the mode shape m for visualization purposes:
        m=3;
        modeShape=squeeze(U_POD(m,:,:));
        figure(); imagesc(modeShape);
        caxis([-max(abs(modeShape(:))) max(abs(modeShape(:)))])
        colormap redblue
        % xlabel('x/L');
        % ylabel('y/H');
        % set(gca,'ydir','normal')
        % title(['Mode Shape ' num2str(m,'%0.0f')]);
    end

    %% Sanity Check! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % pause;

    %% Save Bool %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % if SaveBool
    %     saveStr = sprintf('Stereo_%ideg_ReD%iK_VortZ.fig',phi,round(ReD(redI)/1000,-1));
    %     figure(69); saveas(gcf,saveStr);
    % 
    %     saveStr = sprintf('Stereo_%ideg_ReD%iK_VMag.fig',phi,round(ReD(redI)/1000,-1));
    %     figure(70); saveas(gcf,saveStr);
    % 
    %     saveStr = sprintf('Stereo_%ideg_ReD%iK_Std.fig',phi,round(ReD(redI)/1000,-1));
    %     figure(71); saveas(gcf,saveStr);
    % end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Standard Dev. of Vort. Locs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [STX,STY] = std2D(UX,UY,stdN)

    % figure(1); clf; hold on;
    % plot(UX,UY,'k.')

    if stdN == 1; Pcnt = 2*0.341;
    elseif stdN == 2; Pcnt = 2*(0.341+0.136);
    else; Pcnt = 2*(0.341+0.136+0.023); 
    end

    binwidth = 4; 
    XMin = floor(min(UX)); XMax = ceil(max(UX));
    YMin = floor(min(UY)); YMax = ceil(max(UY));

    XBinCenters = XMin:binwidth:XMax; YBinCenters = YMin:binwidth:YMax;
    XBinIndices = XMin-binwidth/2:binwidth/2:XMax+binwidth/2;
    YBinIndices = YMin-binwidth/2:binwidth/2:YMax+binwidth/2;
    BinDepth = NaN(size(length(YBinCenters),length(XBinCenters)));

    [XBinCenters,YBinCenters] = meshgrid(XBinCenters,YBinCenters);
    
    xcount = 0; 
    for xbin = 2:2:length(XBinIndices)-1
    xcount = xcount+1;
    ycount = 0;
    for ybin = 2:2:length(YBinIndices)-1
        ycount = ycount+1;
        xmin = XBinIndices(xbin-1); xmax = XBinIndices(xbin+1);
        ymin = YBinIndices(ybin-1); ymax = YBinIndices(ybin+1);
        xIdx = (UX>xmin).*(UX<xmax);
        yIdx = (UY>ymin).*(UY<ymax);
        BinDepth(ycount,xcount) = sum(xIdx.*yIdx);
    end
    end

    % figure(2); clf; hold on;
    % contourf(XBinCenters,YBinCenters,BinDepth)

    R = linspace(0,50,100); stdR = NaN(size(1:365));
    for deg = 1:365
        xtemp = R*cosd(deg); ytemp = R*sind(deg);
        vtemp = interp2(XBinCenters,YBinCenters,BinDepth,xtemp,ytemp);
        percent = cumsum(vtemp)/sum(vtemp,'omitmissing');
        [~,stdI] = min(abs(percent-Pcnt));
        stdR(deg) = R(stdI);
    end
    STX = stdR.*cosd(1:365); STY = stdR.*sind(1:365);
end

% Plot Circle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = PlotCircle(xo,yo,R)
    phi = 0:0.5:365;
    x = R*cosd(phi) + xo;
    y = R*sind(phi) + yo;

    hold on; plot(x,y,'k--','LineWidth',2);
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

% Fills in NaN Holes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% ROTZ  Create a 3x3 rotation matrix about the Z-axis. %%%%%%%%%%%%%%%%%%%%
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
