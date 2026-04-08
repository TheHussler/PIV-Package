%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Determine Recirculation Region Bounds %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Rhylan A Huss %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% August 18th 2023 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This code serves as scratch paper for a delerious mind. 
% - The aforementioned delerious mind

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define Recirculation Region Bounds %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% File Path %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% path = {...
    % 'R:\01 - Lofted Cylinder\03 - LSWT HSPIV\ProcessedDataSetsMat\PlanarPIVData_ReD40K_nImg_1000_R03.mat'}
    % 'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\ProcessedDataSetsMat\PlanarPIVData_ReD110K_nImg_1000_R03.mat';
    % 'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\ProcessedDataSetsMat\PlanarPIVData_ReD120K_nImg_1000_R03.mat';
    % 'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\ProcessedDataSetsMat\PlanarPIVData_ReD130K_nImg_1000_R03.mat';
    % 'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\ProcessedDataSetsMat\PlanarPIVData_ReD140K_nImg_1000_R03.mat';
    % 'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\ProcessedDataSetsMat\PlanarPIVData_ReD150K_nImg_1000_R03.mat';
    % 'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\ProcessedDataSetsMat\PlanarPIVData_ReD160K_nImg_1000_R03.mat';
    % 'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\ProcessedDataSetsMat\PlanarPIVData_ReD170K_nImg_1000_R03.mat';
    % 'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\ProcessedDataSetsMat\PlanarPIVData_ReD180K_nImg_1000_R03.mat';
    % 'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\ProcessedDataSetsMat\PlanarPIVData_ReD190K_nImg_1000_R03.mat';
    % 'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\ProcessedDataSetsMat\PlanarPIVData_ReD200K_nImg_1000_R03.mat'};

path = {...
    'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\ProcessedDataSetsMat\PlanarPIVData_45_ReD40K_nImg_1000_R03.mat'};
%     'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\ProcessedDataSetsMat\PlanarPIVData_45_ReD120K_nImg_1000_R03.mat';
%     'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\ProcessedDataSetsMat\PlanarPIVData_45_ReD150K_nImg_1000_R03.mat';
%     'A:\01 - Lofted Cylinder\03 - LSWT HSPIV\ProcessedDataSetsMat\PlanarPIVData_45_ReD190K_nImg_1000_R03.mat'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define Recirculation Region Bounds %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for fileI = 1%:length(path)
    clearvars -except path fileI
    %% Defaults %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Colormaps
    mustang = [1 1 1; jet];
    temp = redblue(540);
    blue = [temp(1:270,:)];
    red = [temp(270:end,:)];
    
    % Variables
    D = 146.05; % mm 
    phi = 45;
    rho = 1.18; %kg/m^3
    dyn_vis = 1.81E-5; %kg/m^3
    nu = dyn_vis/rho;
    
    % Default Figure Formats
    set(groot,'defaultAxesFontName','Times New Roman')
    set(groot,'defaultAxesFontSize',14)
    set(groot,'defaultAxesFontWeight','Bold')
    set(groot,'defaultFigurePosition',[100 100 800 600])

    %% Read in Dataset %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    load(path{fileI}); clc;
    
    Winf = Data.Winf;
    ReD = Winf*(D/1000)/nu; fprintf('Reading Re_D: %iK...\n',round(ReD/1000,-1));
    z = Data.Coordinates.z; y = Data.Coordinates.y; [Z,Y] = meshgrid(z,y);
    Mask = Data.Mask;
    Vz = Data.Vec.Avg.Vz; Vy = Data.Vec.Avg.Vy; Vmag = Data.Vec.Avg.Vmag;

    % Gauss filter
    sigma = 1.25;
    Vz = imgaussfilt(Vz,sigma); Vy = imgaussfilt(Vy,sigma);
    Vmag = imgaussfilt(Vmag,sigma);
    
    % Rounded-Edge Model Vertices
    R = (18.487/60)*D;     %mm
    phj = (0:0.1:phi)';  %deg
    modelVert = [R*cosd(90-phj)-R*tand(phi/2),R*sind(90-phj)-R];
    whiteWash = [min(z)/D,max(y)/D; 0.07,max(y)/D; -0.27,0; modelVert/D; 1/tand(phi),-1; 1/tand(phi),min(y)/D; min(z)/D, min(y)/D];
    modelVert = [min(z)/D,0; modelVert/D; 1/tand(phi),-1; min(z)/D, -1];
    
    vortx = VortX(Z./D,Y./D,Vz./Winf,Vy./Winf);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PLOT PROFILES FOR SAVING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    skip = 4;
    
    % Plot a Visual for Sanity ... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plotVorticity(z/D,y/D,vortx,phi,0.5,Mask,modelVert,whiteWash,blue)
    hold on; quiver(Z(1:skip:end,1:skip:end)/D,Y(1:skip:end,1:skip:end)/D,...
        Vz(1:skip:end,1:skip:end),Vy(1:skip:end,1:skip:end));

    plotVelocity(z/D,y/D,Vmag/Winf,phi,0.5,Mask,modelVert,whiteWash,mustang)
    hold on; quiver(Z(1:skip:end,1:skip:end)/D,Y(1:skip:end,1:skip:end)/D,...
        Vz(1:skip:end,1:skip:end),Vy(1:skip:end,1:skip:end));

    % pause;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% EDGE PROFILES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % function [] = RecircBounds()

    skip = 4;

    % Plot a Visual for Sanity ... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plotVorticity(z/D,y/D,vortx,phi,0.5,Mask,modelVert,whiteWash,blue)
    hold on; q = quiver(Z(1:skip:end,1:skip:end)/D,Y(1:skip:end,1:skip:end)/D,...
        Vz(1:skip:end,1:skip:end),Vy(1:skip:end,1:skip:end)); q.Color = [0.75 0.75 0.75];

    % I need to find the vertices consistent with the PIV Data grid system.
    % What is the dx? <--- This is known already from the Q derivatives
    % Where is the start of the curve?

    %% For The Rounded-Edge
    % The Coords up to the curve:
    R = 18.487;      %mm
    phj = (0:phi)';  %deg
    modelVert = [R*cosd(90-phj)-R*tand(phi/2),R*sind(90-phj)-R];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Determine the Coordinates Prior to the Curve/Radius %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sze = size(Vz);
    [~,jt]=max(~isnan(flip(Mask)),[],1); jt = sze(1)-jt;

    % Determine how many points are prior to the curve
    sn = length(z(z <= modelVert(1,1))');
    % Preallocate sCurve Array
    sCurve = NaN(length(z),4);
    % sCurve = [zloc, yloc, Perpendicular Angle, y-index]
    sCurve(1:sn,1:3) = ...
        [z(z <= modelVert(1,1))',zeros(length(jt(1:sn)'),1),repmat(90,sn,1)];
    sCurve(:,4) = jt;
    % sCurve(1:sn,1:3) = ...
        % [z(z <= modelVert(1,1))',zeros(length(y(jt(1:sn))'),1),repmat(90,sn,1)];

    S(1:sn) = z(z <= modelVert(1,1))' - max(z(z <= modelVert(1,1))) + R*(-22.5)*pi/180;

    for it = find(z > modelVert(1,1) & z < modelVert(end,1))
        xt = z(it); yt = y(jt(it));
        sPhi = atand(abs((yt+R)/(xt+R*tand(phi/2))));
        xt = R*cosd(sPhi)-R*tand(phi/2); yt = R*sind(sPhi)-R;
        sCurve(it,1:3) = [xt,yt,sPhi];

        S(it) = R*(67.5 - sPhi)*pi/180;
    end

    sn = it+1 : it + length(z(z >= modelVert(end,1) & z <= D/tand(phi)));

    %What is this magic I perform? Perpindicular shit.
    b = y(jt(sn))' - (z(sn)'/tand(phi));
    zo = -b/(tand(phi) + 1/tand(phi));
    yo = -tand(phi)*zo;
    sCurve(sn,1:3) = ...
        [zo,yo,repmat(90-phi,length(sn),1)];
    S(sn) = sqrt((sCurve(sn,1) - (R*cosd(phi)-R*tand(phi/2))).^2 + (sCurve(sn,2) - (R*sind(phi)-R)).^2) + S(it);
    S = S';


    %% Take s-Perpendicualr Cuts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nLine = 200; ds = 0.6;
    sLocs = linspace(0.01,ds,nLine);
    sepcheck = 0; rSL = 1;

    % Preallocate Looped Variables
    Vmaginterp  = NaN(nLine+1,length(z(z <= D/tand(phi))));
    Vxinterp    = Vmaginterp;
    Vyinterp    = Vmaginterp;
    Vsinterp    = Vmaginterp;
    VsPERPinterp    = Vmaginterp;
    sPerpLinex  = Vmaginterp;
    sPerpLiney  = Vmaginterp;

    startx = -0*D;
    [~,ri] = min(abs(z-startx));

    [~,Edge] = min(abs(z-(R*cosd(90)-R*tand(phi/2))));

    % B = figure(102);
    grid on; grid minor; hold on;
    ranOnce = false; ShearVirgin = true;
    for r = ri:length(z(z <= D/tand(phi)))

        % Calculate and Store locations of the perpendicular cut line
        sPerpLine = rotz(sCurve(r,3)) * [sLocs*D;zeros(1,nLine);zeros(1,nLine)] + repmat([sCurve(r,1);sCurve(r,2);0],1,nLine);
        sPerpLine = sPerpLine(1:2,:)'/D;
        sPerpLinex(2:end,r) = sPerpLine(:,1); sPerpLiney(2:end,r) = sPerpLine(:,2);
        sPerpLiney(1,1:length(z(z <= D/tand(phi)))) = zeros(1,length(z(z <= D/tand(phi)))); sPerpLiney(1,1:length(z(z <= D/tand(phi)))) = zeros(1,length(z(z <= D/tand(phi))));

        % Interpolates the Slice!
        Vmaginterp(2:end,r) = interp2(z,y,Vmag,sPerpLine(:,1)*D,sPerpLine(:,2)*D,'linear');
        Vmaginterp(1,1:length(z(z <= D/tand(phi)))) = zeros(1,length(z(z <= D/tand(phi))));
        Vxinterp(2:end,r) = interp2(z,y,Vz,sPerpLine(:,1)*D,sPerpLine(:,2)*D,'linear');
        Vxinterp(1,1:length(z(z <= D/tand(phi)))) = zeros(1,length(z(z <= D/tand(phi))));
        Vyinterp(2:end,r) = interp2(z,y,Vy,sPerpLine(:,1)*D,sPerpLine(:,2)*D,'linear');
        Vyinterp(1,1:length(z(z <= D/tand(phi)))) = zeros(1,length(z(z <= D/tand(phi))));

        % Rotate to Find s-Directional Component!
        Vsinterp(:,r) = cosd(90-sCurve(r,3))*Vxinterp(:,r) - sind(90-sCurve(r,3))*Vyinterp(:,r);
        
        nanidx = find(~isnan(Vsinterp(:,r))==1);
        clear foot feet;
        foot = fittype("a*x^2 + b*x",dependent="y",independent="x",coefficients=["a","b"]);
        xfit = [1,nanidx(2):nanidx(2)+10];
        feet = fit(xfit',Vsinterp(xfit,r),foot,'StartPoint',[0,0]);
        % Fill in
        xfit = 2:nanidx(2)-1; yfit = feet.a.*xfit.^2 + feet.b.*xfit;
        Vsinterp(xfit,r) = yfit;
        figure(71); clf; plot(Vsinterp(:,r))
        hold on; plot(xfit,yfit)
        

        % Rotate to Find s-PERP-Directional Component!
        VsPERPinterp(:,r) = cosd(-sCurve(r,3))*Vxinterp(:,r) - sind(-sCurve(r,3))*Vyinterp(:,r);

        % Shows Where The Cut is on the Model
        figure(101); hold on;
        try delete(h); 
        catch; fprintf('Oopsies!\n');
        end
        h = plot(sPerpLine(:,1),sPerpLine(:,2),'r-','LineWidth', 1.5);

        % Calculate BL Max %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        EdgeMask = sPerpLiney < 0.02;
        VsTemp = EdgeMask.*Vsinterp(:,r);
        VsTemp(isnan(VsTemp)) = 0;

        neighbors = ones(size(VsTemp(:,r))); pad = 4;
        if ~ShearVirgin
            neighbors(1:ui-pad,ui+pad:end) = 0;
        end
        ShearVirgin = false;

        [UMax,ui] = max(neighbors.*VsTemp(:,r));

        delta(r,1) = sLocs(ui-1)*D/1000; %m
        delta(r,2) = UMax; % mps
        delta(r,3) = sPerpLinex(ui,r)*D/1000; %m
        delta(r,4) = sPerpLinex(2,r)*D/1000; %m


        figure(101)
        plot(sPerpLinex(ui,r),sPerpLiney(ui,r),'ko')

        % Calculate a Pohlhausen fit by Finite Difference Method %%%%%%%%%%
        % dUdx = ((11/6) * delta(r,2) + (-3) * delta(r-1,2) + (3/2) * delta(r-2,2) + (-1/3) * delta(r-3,2))/(delta(r,3) - delta(r-1,3));
        % Lambda = (delta(r,1)^2/nu) * dUdx;
        % dpdx = -Lambda * ((dyn_vis*UMax)/delta(r,2))/delta(r,2);
        % eta = linspace(0,1);
        % uPohl = UMax * ((2*eta - 2*eta.^3 + eta.^4) + (Lambda/6)*(eta - 3*eta.^2 + 3*eta.^3-eta.^4));

        % Forces a Pohlhausen fit to the data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        clear foot feet;
        foot = fittype("1-(1+x)*(1-x)^3 + (a/6)*x*(1-x)^3",dependent="y",independent="x",coefficients="a");
        eta = linspace(0,1); 

        nanI = isnan(Vsinterp(2:ui,r)/UMax);
        VsFit = rmmissing(Vsinterp(2:ui,r)/UMax);
        sFit =  (sLocs(1:ui-1)*D/1000)'/delta(r,1); sFit(nanI)=NaN;
        sFit =  rmmissing(sFit);

        feet = fit(sFit,VsFit,foot,'StartPoint',0);
        uPohl = 1-(1+eta).*(1-eta).^3 + (feet.a/6).*eta.*(1-eta).^3;

        if r == Edge
            deltaE = delta(r,1);
            dispThicknessE = trapz(eta*delta(r,1),1-uPohl); %delta1
            momThicknessE = trapz(eta*delta(r,1),uPohl.*(1-uPohl)); %delta2
            HE = dispThicknessE/momThicknessE;
            [maxShearE,maxShearLocE] = max(((uPohl(3:end) - uPohl(1:end-2))*UMax) ./ ((eta(3:end) - eta(1:end-2))*delta(r,1)));
            maxShearLocE = eta(maxShearLocE)*delta(r,1);

            % Shows the Cut %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % K. Flay
            figure(32); clf; hold on;
            xline(0,'k-'); xline(0.2,'r--'); xline(0.8,'r--'); 
            xline(1.0,'k--'); xline(delta(r,2)/Winf,'k-.');
            % title(sprintf('z/D Location = %f',z(r)/D));
            g(1) = plot(Vxinterp(:,r)/UMax,[0,sLocs]/(delta(r,1)/(D/1000)),'r--','LineWidth', 1.5);
            g(2) = plot(Vsinterp(:,r)/UMax,[0,sLocs]/(delta(r,1)/(D/1000)),'r-','LineWidth', 1.5);
            % g(3) = plot(Vyinterp(:,r)/Winf,[0,sLocs],'b--','LineWidth', 1.5);
            % g(4) = plot(VsPERPinterp(:,r)/Winf,[0,sLocs],'b-','LineWidth', 1.5);
            g(5) = plot(uPohl,eta,'b-'); g(5).LineWidth = 1.5;
            legend('','','','','','U_z','U_s','Pohlhausen Fit')
            xlim([-0.4 1.35]); ylim([0 1])
            xlabel('U/W_s|_{max}'), ylabel('\perps/\delta_e')

            saveFile = sprintf('RecircReg_%ideg_ReD%iK_EdgeBLProfile.fig',phi,round(ReD/1000,-1));
            savefig(gcf,saveFile)
            fprintf('Image %s\n saved.\n',saveFile)

        end


        if (uPohl(2)<0 || Vsinterp(2,r)<0) && ~ranOnce
            sepLocAll(r,:) = [sPerpLine(1,1),sPerpLine(1,2)];
            RSep = r;
            deltaS = delta(r,1);
            dispThicknessS = trapz(eta*delta(RSep,1),1-uPohl); %delta1
            momThicknessS  = trapz(eta*delta(RSep,1),uPohl.*(1-uPohl)); %delta2
            HS = dispThicknessS/momThicknessS;
            [maxShearS,maxShearLocS] = max(((uPohl(3:end) - uPohl(1:end-2))*UMax) ./ ((eta(3:end) - eta(1:end-2))*delta(RSep,1)));
            maxShearLocS = eta(maxShearLocS)*delta(RSep,1);
            ranOnce = true;

            % Shows the Cut %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % K. Flay
            figure(31); clf; hold on;
            xline(0,'k-'); xline(0.2,'r--'); xline(0.8,'r--'); 
            xline(1.0,'k--'); xline(delta(r,2)/Winf,'k-.');
            title(sprintf('z/D Location = %f',z(r)/D));
            g(1) = plot(Vxinterp(:,r)/UMax,[0,sLocs]/(delta(r,1)/(D/1000)),'r--','LineWidth', 1.5);
            g(2) = plot(Vsinterp(:,r)/UMax,[0,sLocs]/(delta(r,1)/(D/1000)),'r-','LineWidth', 1.5);
            % g(3) = plot(Vyinterp(:,r)/Winf,[0,sLocs],'b--','LineWidth', 1.5);
            % g(4) = plot(VsPERPinterp(:,r)/Winf,[0,sLocs],'b-','LineWidth', 1.5);
            g(5) = plot(uPohl,eta,'b-'); g(5).LineWidth = 1.5;
            legend('','','','','','U_z','U_s','Pohlhausen Fit')
            xlim([-0.4 1.35]); ylim([0 1])
            xlabel('U/W_s|_{max}'), ylabel('\perps/\delta_s')

            saveFile = sprintf('RecircReg_%ideg_ReD%iK_SepBLProfile.fig',phi,round(ReD/1000,-1));
            savefig(gcf,saveFile)
            fprintf('Image %s\n saved.\n',saveFile)

        end

        % Shows the Cut %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % K. Flay
        figure(102); clf; hold on;
        xline(0,'k-'); xline(0.2,'r--'); xline(0.8,'r--'); 
        xline(1.0,'k--'); xline(delta(r,2)/Winf,'k-.');
        title(sprintf('z/D Location = %f',z(r)/D));
        g(1) = plot(Vxinterp(:,r)/UMax,[0,sLocs]/sLocs(ui-1),'r--','LineWidth', 1.5);
        g(2) = plot(Vsinterp(:,r)/UMax,[0,sLocs]/sLocs(ui-1),'r-','LineWidth', 1.5);
        % g(3) = plot(Vyinterp(:,r)/Winf,[0,sLocs],'b--','LineWidth', 1.5);
        % g(4) = plot(VsPERPinterp(:,r)/Winf,[0,sLocs],'b-','LineWidth', 1.5);
        g(5) = plot(uPohl,eta,'b-'); g(5).LineWidth = 1.5;
        legend('','','','','','U_z','U_s')
        xlim([-0.4 1.35]); ylim([0 1])
        xlabel('U/W_{\infty}'), ylabel('\perps/D')

        % pause;

        if mean(Vsinterp(1:5,r),'omitnan')/Winf<0 && sepcheck ~= 1
            % figure(101); hold on;
            % plot(sPerpLine(1,1),sPerpLine(1,2),'ko','MarkerFaceColor','k')
            shearLayer = []; MeanDivSL = [];
            sepcheck = 1; reattachcheck = 0;
        end

        % Find Shear Layer %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        zeroIndex = find(Vsinterp(:,r)/Winf <= 0);
        pointTwoIndex = find(Vsinterp(:,r)/Winf <= 0.2);
        pointFourIndex = find(Vsinterp(:,r)/Winf <= 0.4);
        pointSixIndex = find(Vsinterp(:,r)/Winf <= 0.6);
        pointEightIndex = find(Vsinterp(:,r)/Winf <= 0.8);

        % zeroPoint(rSL,:) = [r,sLocs(zeroIndex(end)),sPerpLine(zeroIndex(end),1),sPerpLine(zeroIndex(end),2)];
        % PointTwo(r,:) = [r,sLocs(pointTwoIndex(end)),sPerpLine(pointTwoIndex(end),1),sPerpLine(pointTwoIndex(end),2)];
        % PointFour(r,:) = [r,sLocs(pointFourIndex(end)),sPerpLine(pointFourIndex(end),1),sPerpLine(pointFourIndex(end),2)];
        % PointSix(r,:) = [r,sLocs(pointSixIndex(end)),sPerpLine(pointSixIndex(end),1),sPerpLine(pointSixIndex(end),2)];
        % PointEight(r,:) = [r,sLocs(pointEightIndex(end)),sPerpLine(pointEightIndex(end),1),sPerpLine(pointEightIndex(end),2)];
        
        ShearEdge(r,:) = [r,delta(r,1)*1000/D,sPerpLine(ui-1,1),sPerpLine(ui-1,2)];
        ShearEdge(ShearEdge == 0) = NaN; ShearEdge = rmmissing(ShearEdge);

        % Plot For Visual
        % figure(101); hold on;
        % plot(zeroPoint(:,3),zeroPoint(:,4),'k--','LineWidth', 1)
        % plot(PointTwo(:,3),PointTwo(:,4),'k--','LineWidth', 1)
        % plot(PointFour(:,3),PointFour(:,4),'k--','LineWidth', 1)
        % plot(PointSix(:,3),PointSix(:,4),'k--','LineWidth', 1)
        % plot(PointEight(:,3),PointEight(:,4),'k--','LineWidth', 1)
        % plot(sPerpLine(ui-1,1),sPerpLine(ui-1,2),'k.')

        if sepcheck == 1% && reattachcheck ~= 1

            % Find the Mean Separating Streamline %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            meanDivPoint = NaN(nLine,1);
            for q = 1:nLine
                VsTemp = Vsinterp(1:q,r)/Winf;
                VsTemp(isnan(VsTemp)) = 0;
                meanDivPoint(q) = trapz(VsTemp);
            end
            try
                meanDivIndex = find(flip(meanDivPoint<0));
                meanDivIndex = (nLine+1)-meanDivIndex(1);

                meanDivS = interp1(meanDivPoint(meanDivIndex:meanDivIndex+1),sLocs(meanDivIndex:meanDivIndex+1),0);
                meanDivX = interp1(meanDivPoint(meanDivIndex:meanDivIndex+1),sPerpLine(meanDivIndex:meanDivIndex+1,1),0);
                meanDivY = interp1(meanDivPoint(meanDivIndex:meanDivIndex+1),sPerpLine(meanDivIndex:meanDivIndex+1,2),0);
                meanDivVs = interp1(meanDivPoint(meanDivIndex:meanDivIndex+1),Vsinterp(meanDivIndex:meanDivIndex+1,r),meanDivS);

                % MeanDivSL = [INDEX, Sloc, X, Y];
                MeanDivSL(rSL,:) = [r,meanDivS,meanDivX,meanDivY];
                MeanDivSL(MeanDivSL == 0) = NaN;

                figure(101);
                plot(MeanDivSL(:,3),MeanDivSL(:,4),'k-.','LineWidth', 2)

            catch
                try
                    if (sum(Vsinterp(:,r)<0) == 0 || isnan(MeanDivSL(rSL-1,3))) && sPerpLinex(2,r)>0.6
                        
                        if isnan(MeanDivSL(rSL,3))
                            MeanDivSL = MeanDivSL(1:rSL-1);
                        end
                        figure(101);
                        plot(MeanDivSL(:,3),MeanDivSL(:,4),'k-.','LineWidth', 2)
                        plot(sPerpLine(1,1),sPerpLine(1,2),'ko','MarkerFaceColor','k')
                        reAttLoc = [sPerpLine(1,1),sPerpLine(1,2)];
                        reattachcheck = 1;
                        break
                    end
                catch
                    figure(101);
                    plot(MeanDivSL(:,3),MeanDivSL(:,4),'k-.','LineWidth', 2)
                    plot(sPerpLine(1,1),sPerpLine(1,2),'ko','MarkerFaceColor','k')
                    reAttLoc = [sPerpLine(1,1),sPerpLine(1,2)];
                    reattachcheck = 1;
                    break
                end
            end
            rSL = rSL+1;
        end
    end

    pause;

    %%
    if exist('sepLocAll', 'var') && exist('reAttLoc', 'var') 
        % close 101 
        sepLocAll(sepLocAll==0)=NaN;
        sepLoc = rmmissing(sepLocAll); sepLoc = sepLoc(1,:);

        sepBubbleLength = reAttLoc(1) - sepLoc(1);
    
        %% Shear Layer and MDSL Fit %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % foot = fittype("a*x+b",dependent="y",independent="x",coefficients=["a","b"]);
        % foot = fittype("a*x+b*x^2+c*x^3+d*x^4",dependent="y",independent="x",coefficients=["a","b","c","d"]);
        
        figure(22); clf; hold on

        % Shear Layer Fit %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        clear foot;
        foot = fittype("a+b*x+c*x^2",dependent="y",independent="x",coefficients=["a","b","c"]);
        temp = (sPerpLinex(2,ShearEdge(:,1))-sepLoc(1))>0;
        Why = temp'.*ShearEdge(:,2); Why(Why==0)=NaN; Why = rmmissing(Why); 
        Ex = temp.*(sPerpLinex(2,ShearEdge(:,1))-sepLoc(1));
        Ex(Ex==0)=NaN; Ex = rmmissing(Ex); 
        Qx = linspace(0,Ex(end)); 
        ShearFit = fit(Ex',Why,foot);

        % Remove Outliers and Re-Fit
        clear inliners foot
        residuals = abs(Why - ShearFit(Ex'));
        threshold = 2 * std(residuals);
        Shearinliers = residuals < threshold;
        Ex = Ex(Shearinliers);
        Why = Why(Shearinliers);
        foot = fittype("a+b*x+c*x^2",dependent="y",independent="x",coefficients=["a","b","c"]);
        Qx = linspace(0,Ex(end)); 
        ShearFit = fit(Ex',Why,foot);
    
        plot(Ex,Why,'rx')
        plot(ShearFit,'r-')

        % Shear Layer Save Specifics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        RecircReg.Sep.Shear.Fit = ShearFit;
        RecircReg.Sep.Shear.X = Ex;
        RecircReg.Sep.Shear.Y = Why;

        % MDSL Fit %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        clear Ex Why foot;
        foot = fittype("a*x+b*x^2+c*x^3+d*x^4",dependent="y",independent="x",coefficients=["a","b","c","d"]);
        temp = (sPerpLinex(2,rmmissing(MeanDivSL(:,1)))-sepLoc(1))>0;
        Why = temp'.*rmmissing(MeanDivSL(:,2)); Why(Why==0)=NaN; Why = rmmissing(Why); 
        Ex = temp.*(sPerpLinex(2,rmmissing(MeanDivSL(:,1)))-sepLoc(1));
        Ex(Ex==0)=NaN; Ex = rmmissing(Ex); Qx = linspace(0,Ex(end)); 
        MDSLFit = fit(Ex',Why,foot);

        % Remove Outliers and Re-Fit
        clear inliners foot
        
        residuals = abs(Why - MDSLFit(Ex'));
        threshold = 2 * std(residuals);
        MDSLinliers = residuals < threshold;
        Ex = Ex(MDSLinliers);
        Why = Why(MDSLinliers);
        foot = fittype("a*x+b*x^2+c*x^3+d*x^4",dependent="y",independent="x",coefficients=["a","b","c","d"]);
        Qx = linspace(0,Ex(end)); 
        MDSLFit = fit(Ex',Why,foot);

        MDSLArea = integral(@(x) MDSLFit(x), 0, fzero(@(x) MDSLFit(x), 1), 'ArrayValued', true);
    
        plot(Ex,Why,'kx')
        plot(MDSLFit,'k-')

        % MDSL Save Specifics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        RecircReg.Sep.MDSL.Fit = ShearFit;
        RecircReg.Sep.MDSL.X = Ex;
        RecircReg.Sep.MDSL.Y = Why;
    
        % xline(sepLoc(1),'k--')
        xline(reAttLoc(1),'k--')
    
        L = legend('Shear Layer Upper Bound',sprintf('Shear Layer Fit = (%.1f)+(%.1f)x+(%.1f)x^2',ShearFit.a,ShearFit.b,ShearFit.c),...
            'M.D.S.L.',sprintf('M.D.S.L. Fit = (%.1f)x+(%.1f)x^2+(%.1f)x^3+(%.1f)x^4',MDSLFit.a,MDSLFit.b,MDSLFit.c,MDSLFit.d));
        L.NumColumns = 2; L.Location = "southoutside"; L.Orientation = 'horizontal';
        
        xlabel('(z-z_s)/D'); ylabel('y/D')
        xlim([0 1.2]); ylim([0 0.5]); daspect([1 1 1]);

        saveFile = sprintf('RecircReg_%ideg_ReD%iK_MDSLProfile.fig',phi,round(ReD/1000,-1));
        savefig(gcf,saveFile)
        fprintf('Image %s\n saved.\n',saveFile)
    
        %% Shear Layer Growth %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        IndAlign = zeros(size(ShearEdge(:,1))); MeanDivSLComp = rmmissing(MeanDivSL(:,1:2));
        for I = 1:length(MeanDivSLComp(:,1))
            IndAlign = IndAlign+(ShearEdge(:,1) == MeanDivSLComp(I,1));
            if sum((ShearEdge(:,1) == MeanDivSLComp(I,1)))==0
                MeanDivSLComp(I,1) = NaN;
            end
        end
        MeanDivSLComp = rmmissing(MeanDivSLComp);
        ShearEdgeComp = IndAlign.*ShearEdge(:,2); ShearEdgeComp(ShearEdgeComp==0) = NaN;
        ShearEdgeComp = rmmissing(ShearEdgeComp);
        ShearLayerSize = ShearEdgeComp-MeanDivSLComp(:,2);

        figure(23); clf; hold on;
        plot(sPerpLinex(2,rmmissing(MeanDivSLComp(:,1)))-sepLoc(1),ShearLayerSize,'kx')
        Qy = -MDSLFit.d * Qx.^4 + -MDSLFit.c * Qx.^3 + (-MDSLFit.b+ShearFit.c) * Qx.^2 + (-MDSLFit.a+ShearFit.b) * Qx + ShearFit.a;
        plot(Qx,Qy,'k-')

        xlabel('(z-z_s)/D'); ylabel('\delta_{SL}/D')
        xlim([0 1.2]); ylim([0 0.5]); daspect([1 1 1]);

        saveFile = sprintf('RecircReg_%ideg_ReD%iK_ShearGrowthProfile.fig',phi,round(ReD/1000,-1));
        savefig(gcf,saveFile)
        fprintf('Image %s\n saved.\n',saveFile)
        
        % Shear Layer Growth Specifics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        RecircReg.Sep.Shear.Growth.FitX = Qx;
        RecircReg.Sep.Shear.Growth.FitY = Qy;
        RecircReg.Sep.Shear.Growth.X = sPerpLinex(2,rmmissing(MeanDivSLComp(:,1)))-sepLoc(1);
        RecircReg.Sep.Shear.Growth.Size = ShearLayerSize;

        % Shear Layer Growth Rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        dx = sPerpLinex(2,rmmissing(MeanDivSLComp(:,1)))-sepLoc(1);
        DSLdx = (ShearLayerSize(2:end) - ShearLayerSize(1:end-1))./(dx(2:end) - dx(1:end-1))';

        figure(24); clf; hold on;
        plot(dx(2:end),DSLdx,'kx')
        Qy = -4*MDSLFit.d * Qx.^3 + -3*MDSLFit.c * Qx.^2 + 2*(-MDSLFit.b+ShearFit.c) * Qx + (-MDSLFit.a+ShearFit.b);
        plot(Qx,Qy,'k-')

        xlabel('(z-z_s)/D'); ylabel('\partial\delta_{SL}/\partialz')
        xlim([0 1.2]); ylim([-3 5]);

        saveFile = sprintf('RecircReg_%ideg_ReD%iK_ShearGrowthRateProfile.fig',phi,round(ReD/1000,-1));
        savefig(gcf,saveFile)
        fprintf('Image %s\n saved.\n',saveFile)

        % Shear Layer Growth Rate Save Specifics %%%%%%%%%%%%%%%%%%%%%%%%%%
        RecircReg.Sep.Shear.GrowthRate.FitX = Qx;
        RecircReg.Sep.Shear.GrowthRate.FitY = Qy;
        RecircReg.Sep.Shear.GrowthRate.X = dx(2:end);
        RecircReg.Sep.Shear.GrowthRate.Size = DSLdx;

        %% Plot Model Display %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        modelVert = [R*cosd(90-phj)-R*tand(phi/2),R*sind(90-phj)-R];
        whiteWash = [min(z)/D,max(y)/D; 0.07,max(y)/D; -0.27,0; modelVert/D; 1/tand(phi),-1; 1/tand(phi),min(y)/D; min(z)/D, min(y)/D];
        modelVert = [min(z)/D,0; modelVert/D; 1/tand(phi),-1; min(z)/D, -1];
        
        plotVorticity(z/D,y/D,vortx,phi,0.5,Mask,modelVert,whiteWash,blue)
        hold on; q = quiver(Z(1:skip:end,1:skip:end)/D,Y(1:skip:end,1:skip:end)/D,...
            Vz(1:skip:end,1:skip:end),Vy(1:skip:end,1:skip:end)); q.Color = [0.75 0.75 0.75];
        plot(sepLoc(1),sepLoc(2),'ko','MarkerFaceColor','k')
        plot(reAttLoc(1),reAttLoc(2),'ko','MarkerFaceColor','k')
    
        Ex = (sPerpLinex(2,ShearEdge(:,1))-sepLoc(1))>0;
        Why = Ex'.*ShearEdge(:,4); Why = Shearinliers.*Why(Why~=0);
        Ex =  Ex'.*ShearEdge(:,3); Ex =  Shearinliers.*Ex(Ex~=0);
        plot(Ex',Why,'rx')
        Ex = (sPerpLinex(2,rmmissing(MeanDivSL(:,1)))-sepLoc(1))>0;
        Why = Ex'.*rmmissing(MeanDivSL(:,4)); Why = MDSLinliers.*Why(Why~=0);
        Ex =  Ex'.*rmmissing(MeanDivSL(:,3)); Ex =  MDSLinliers.*Ex(Ex~=0);
        plot(MeanDivSL(:,3),MeanDivSL(:,4),'kx')

        % Model Plot Save Specifics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        saveFile = sprintf('RecircReg_%ideg_ReD%iK_ModelOverview.fig',phi,round(ReD/1000,-1));
        savefig(gcf,saveFile)
        fprintf('Image %s\n saved.\n',saveFile)

        %% Vorticity Thickness Calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        temp = ((sPerpLinex(2,ShearEdge(:,1))-sepLoc(1))>0)';
        temp = repmat(temp,[1,length(Vsinterp(:,1))])';
        VeeEss = temp.*Vsinterp(:,ShearEdge(:,1));
        idx = find(temp(1,:)>0); idx = idx(1);
        VeeEss = VeeEss(:,idx:end); % mps

        temp = ((sPerpLinex(2,ShearEdge(:,1))-sepLoc(1))>0)';
        VortThicknessEx = temp.*(sPerpLinex(2,ShearEdge(:,1))-sepLoc(1))';
        VortThicknessEx = VortThicknessEx(idx:end); %x/D
        % These are SLocations
        MDSLWhy = MDSLFit(VortThicknessEx); ShearWhy = ShearFit(VortThicknessEx); %x/D

        BLMask = double((delta(ShearEdge(:,1),4)-sepLoc(1)*D/1000)>0);
        BLMask(BLMask==0) = NaN;
        BLThickness = rmmissing(BLMask.*delta(ShearEdge(:,1),1)*1000/D);
        UMax = rmmissing(BLMask.*delta(ShearEdge(:,1),2));
        
        clear ShearLoc MDSLLoc
        for I = 1:length(VortThicknessEx)

            [~,ShearLoc(I)] = min(abs(ShearWhy(I)-sLocs));
            [~,MDSLLoc(I)]  = min(abs(MDSLWhy(I)-sLocs));
            
            figure(34); clf; hold on
            plot(VeeEss(:,I)/UMax(I),[0,sLocs]/BLThickness(I),'r--')
            plot(VeeEss(MDSLLoc(I)+1:ShearLoc(I)+1,I)/UMax(I),sLocs(MDSLLoc(I):ShearLoc(I))'/BLThickness(I),'b-')
            
            Vee = VeeEss(MDSLLoc(I)+1:ShearLoc(I)+1,I)/UMax(I); %W/Winf
            % Ess = sLocs(MDSLLoc(I):ShearLoc(I))'/BLThickness(I); dVeedEss = NaN(size(Ess));
            Ess = sLocs(MDSLLoc(I):ShearLoc(I))'*D/1000; %m
            dVeedEss = NaN(size(Ess));

            % Correct Up To Here
            dVeedEss(2:end-1) = (Vee(3:end)-Vee(1:end-2))./(Ess(3:end)-Ess(1:end-2)); %1/m

            [vortThickness(I),vortThicknessI] = max(dVeedEss);
            vortThickness(I) = vortThickness(I)^(-1); %m
            plot(Vee(vortThicknessI),(Ess(vortThicknessI)*1000/D)/BLThickness(I),'ko')

            xlim([-0.4 1.35]); ylim([0 1])
            xlabel('U/W_{\infty}(z)'), ylabel('\perps/D')
            % pause()

        end

        clear Ex Why foot;
        foot = fittype("a+b*x+c*x^2+d*x^3+e*x^4",dependent="y",independent="x",coefficients=["a","b","c","d","e"]); 
        Ex = VortThicknessEx*D/1000; %m
        Why = vortThickness'; %m
        vortThickFit = fit(Ex,Why,foot);

        figure(25); clf; hold on;
        plot(vortThickFit,'r-')

        % Remove Outliers and Re-Fit
        clear inliners foot

        residuals = abs(Why - vortThickFit(Ex));
        threshold = 2 * std(residuals);
        vortThickinliers = residuals < threshold;
        Ex = Ex(vortThickinliers);
        Why = Why(vortThickinliers);
        foot = fittype("a+b*x+c*x^2+d*x^3+e*x^4",dependent="y",independent="x",coefficients=["a","b","c","d","e"]); 
        vortThickFit = fit(Ex,Why,foot);

        figure(25);hold on;
        plot(Ex,Why,'k.')
        plot(vortThickFit,'k-')
        xlim([min(Ex),max(Ex)]); ylim([min(Why),max(Why)]);

        saveFile = sprintf('RecircReg_%ideg_ReD%iK_VortThickness.fig',phi,round(ReD/1000,-1));
        savefig(gcf,saveFile)
        fprintf('Image %s\n saved.\n',saveFile)

        xlim([0 1]); ylim([0 1])
        xlabel('z/D'), ylabel('delta`')

        clear foot
        temp = Ex<0.15*sepBubbleLength;
        Exx = temp.*Ex; Exx(Exx==0) = NaN; Exx = rmmissing(Exx);
        Whyy = temp.*Why; Whyy(Whyy==0) = NaN; Whyy = rmmissing(Whyy);
        foot = fittype("a+b*x",dependent="y",independent="x",coefficients=["a","b"]); 
        vortThickGrowthRateFit = fit(Exx,Whyy,foot);

        plot(vortThickGrowthRateFit,'k-')

        vortThickGrowthRateLinear = vortThickGrowthRateFit.b;
        vortThickGrowthRateFunc = vortThickFit.b+2*vortThickFit.c.*Ex+3*vortThickFit.d*Ex.^2+4*vortThickFit.e*Ex.^3;

        figure(26); clf; hold on;
        yline(vortThickGrowthRateLinear)
        plot(Ex,vortThickGrowthRateFunc,'k-')
        xlim([0 1]); ylim([-3.5 1])

        saveFile = sprintf('RecircReg_%ideg_ReD%iK_VortThicknessGrowthRate.fig',phi,round(ReD/1000,-1));
        savefig(gcf,saveFile)
        fprintf('Image %s\n saved.\n',saveFile)

        % Vort Thickness Save Specifics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        RecircReg.Sep.Shear.VortThickness.X = Ex;
        RecircReg.Sep.Shear.VortThickness.Y = Why;
        RecircReg.Sep.Shear.VortThickness.Fit = vortThickFit;
        RecircReg.Sep.Shear.VortThickness.GrowthRate.Linear = vortThickGrowthRateLinear;
        RecircReg.Sep.Shear.VortThickness.GrowthRate.X = Ex;
        RecircReg.Sep.Shear.VortThickness.GrowthRate.Func = vortThickGrowthRateFunc;

        %% Sep Bubble Save Specifics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        RecircReg.Sep.SepLoc = sepLoc;
        RecircReg.Sep.ReAttachLoc = reAttLoc;
        RecircReg.Sep.SepBubLength = sepBubbleLength;
        RecircReg.Sep.SepBubArea = MDSLArea;
        RecircReg.Sep.delta = deltaS;
        RecircReg.Sep.disThickness = dispThicknessS;
        RecircReg.Sep.momThickness = momThicknessS;
        RecircReg.Sep.H = HS;
        RecircReg.Sep.maxShear = maxShearS;
        RecircReg.Sep.maxShearLoc = maxShearLocS;
        RecircReg.Sep.ReDelta = Winf*deltaS/nu;


    else
        fprintf('\nNo Detected Separation Region\n\n')
    end

    % Beginning Edge BL Save Specifics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    RecircReg.Edge.delta = deltaE;
    RecircReg.Edge.disThickness = dispThicknessE;
    RecircReg.Edge.momThickness = momThicknessE;
    RecircReg.Edge.H = HE;
    RecircReg.Edge.maxShear = maxShearE;
    RecircReg.Edge.maxShearLoc = maxShearLocE;
    RecircReg.Edge.ReDelta = Winf*deltaE/nu;

    % Main Variable Save Specifics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    RecircReg.Variables.deltaX = delta;
    RecircReg.Variables.sPerpLine = sPerpLine;
    RecircReg.Variables.sPerpLinex = sPerpLinex;
    RecircReg.Variables.sPerpLiney = sPerpLiney;
    RecircReg.Variables.Vmaginterp = Vmaginterp;
    RecircReg.Variables.Vsinterp = Vsinterp;
    RecircReg.Variables.Vxinterp = Vxinterp;
    RecircReg.Variables.Vyinterp = Vyinterp;
    RecircReg.Variables.Vx = Vz;
    RecircReg.Variables.Vy = Vy;
    RecircReg.Variables.X = Z;
    RecircReg.Variables.Y = Y;

    % Global Save Specifics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    RecircReg.Global.D = D;
    RecircReg.Global.phi = phi;
    RecircReg.Global.ReD = ReD;
    RecircReg.Global.Winf = Winf;


    %% Save Output Variable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % close all
    saveFile = sprintf('RecircReg_%ideg_ReD%iK_R01.mat',phi,round(ReD/1000,-1));
    fprintf('\nSaving %s... ',saveFile)
    save(saveFile,'RecircReg');
    fprintf('Saved.\n\n')

    if exist('deltaS', 'var')
        fprintf('Finished.\nReD %f\nReDeltaE %f\nReDeltaS %f\n\n',ReD,Winf*deltaE/nu,Winf*deltaS/nu)
    else
        fprintf('Finished.\nReD %f\nReDeltaE %f\n\n',ReD,Winf*deltaE/nu)
    end

    % pause;

    % stop;
end
clc; fprintf('Completed %ideg All Cases.\n',phi)

stop;
%% Smooth Recirc Region
Rotted = rotz(45)*[MeanDivSL(:,3:4),repmat(0,[length(MeanDivSL(:,1)),1])]';

% Parameters
sigma = 5;      % standard deviation of Gaussian
kernel_size = 6*sigma + 1;  % rule of thumb for kernel width
% Create normalized Gaussian kernel
x = linspace(-3*sigma, 3*sigma, kernel_size);
g = exp(-x.^2 / (2*sigma^2));
g = g / sum(g);  % normalize

Rotted(2,:) = conv(Rotted(2,:), g, 'same'); 
% Rotted(2,:) = medfilt1(Rotted(2,:),25);
ReRott = rotz(-45)*Rotted; ReRott = ReRott(1:2,:)'; 

plot(ReRott(:,1),ReRott(:,2),'rp')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot Vorticity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = plotVorticity(z,y,vortx,phi,sigma,mask,modelVert,whiteWash,coloringBook)
    
    vortx(isnan(vortx)) = 0;

    figure(101); clf; hold on; 
    set(gca,'Visible','off'); set(gca,'XTick',[]); set(gca,'YTick',[]);
    ax1 = axes;
    contourf(z,y, mask.*imgaussfilt(vortx,sigma),[-10,0:2.5:30], 'LineStyle','none');
    xlim([-0.2 1.75]); ylim([-1.10 0.25]); daspect([1 1 1]);
    xlabel('z/D'); ylabel('y/D');
    clim([0 30])
    
    hold on; fill(whiteWash(:,1),whiteWash(:,2),[1 1 1],'EdgeColor','none');
    
    view(2); ax2 = axes;
    patch(modelVert(:,1),modelVert(:,2),modelVert(:,2),'LineWidth',2,'EdgeColor',[0.5,0.5,0.5])
    xlim([-0.2 1.75]); ylim([-1.10 0.25]); daspect([1 1 1]);
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

    if phi == 45
        xlim([-0.2 1.2])
        set(gcf, 'Position', [229.5714 379.2857 760.5714 600])
    else
        set(gcf, 'Position', [773 275.8571 987.4286 600]);
    end

end

% Plot Vorticity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = plotVelocity(z,y,vmag,phi,sigma,mask,modelVert,whiteWash,coloringBook)
    
    vmag(isnan(vmag)) = 0;

    figure(102); clf; hold on; 
    set(gca,'Visible','off'); set(gca,'XTick',[]); set(gca,'YTick',[]);
    ax1 = axes;
    contourf(z,y, mask.*imgaussfilt(vmag,sigma),[-10,0:0.05:1.25,5], 'LineStyle','none');
    xlim([-0.2 1.75]); ylim([-1.10 0.25]); daspect([1 1 1]);
    xlabel('z/D'); ylabel('y/D');
    clim([0 1.25])
    
    hold on; fill(whiteWash(:,1),whiteWash(:,2),[1 1 1],'EdgeColor','none');
    
    view(2); ax2 = axes;
    patch(modelVert(:,1),modelVert(:,2),modelVert(:,2),'LineWidth',2,'EdgeColor',[0.5,0.5,0.5])
    xlim([-0.2 1.75]); ylim([-1.10 0.25]); daspect([1 1 1]);
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
    
    if phi == 45
        xlim([-0.2 1.2])
        set(gcf, 'Position', [229.5714 379.2857 760.5714 600])
    else
        set(gcf, 'Position', [773 275.8571 987.4286 600]);
    end

end

% Calculate Vorticity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vortx] = VortX(Z,Y,Vz,Vy)
    % Calculate Vorticity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dx = Z(1,2)-Z(1,1); dy = Y(2,1)-Y(1,1);
    % Calculate Gradients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dwdy = zeros(size(Vz)); dvdx = zeros(size(Vz));
    dwdy(2:end-1,:) = (Vz(3:end,:)-Vz(1:end-2,:))./(2*dy);
    dvdx(:,2:end-1) = (Vy(:,3:end)-Vy(:,1:end-2))./(2*dx);
    % Output Vorticity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vortx = dwdy - dvdx;
end

% Create a 3x3 rotation matrix about the Z-axis %%%%%%%%%%%%%%%%%%%%%%%%%%%
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


