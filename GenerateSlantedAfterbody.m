%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate Slanted Afterbody for X-Foil %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Rhylan A Huss %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% August 18th 2023 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This code serves as scratch paper for a delerious mind. 
% - The aforementioned delerious mind

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;
D = 146.05;
phi = 20;

% Nose Vertices
yNose = D*linspace(-1,0);
xNose = -sqrt(D^2-4*(yNose+0.5*D).^2)-2*D;

% Rounded-Edge Vertices
R = (18.487/60)*D;     %mm
phj = (0:0.1:phi)';  %deg
roundedVert = [R*cosd(90-phj)-R*tand(phi/2),R*sind(90-phj)-R];

% Top Side Vertices
xtop = linspace(-2*D,roundedVert(1,1)); xtop = xtop(2:end-1);
ytop = zeros(size(xtop));

% Bottom Side Vertices
xbot = linspace(-2*D,D/tand(phi)); xbot = xbot(2:end-1);
ybot = -D * ones(size(xbot));

% Slant Vertices
xslant = linspace(roundedVert(end,1), D/tand(phi)); xslant = xslant(2:end);
yslant = linspace(roundedVert(end,2),-D); yslant = yslant(2:end);

xpoints = [flip(xslant)'; flip(roundedVert(:,1)); flip(xtop)'; flip(xNose)'; xbot'] - D/tand(phi) + D*3 + D/tand(phi);
ypoints = [flip(yslant)'; flip(roundedVert(:,2)); flip(ytop)'; flip(yNose)'; ybot'] + D;

xpoints = xpoints / (D*3 + D/tand(phi));
ypoints = ypoints / (D*3 + D/tand(phi));

figure(1); clf;
plot(xpoints,ypoints); daspect([1 1 1])

[~,i] = min(xpoints);

fileid = fopen('20Deg.txt','w');
fprintf(fileid,'%1.6f  %1.6f\n',[xpoints(1:i),ypoints(1:i)]');
fprintf(fileid,'\n');
fprintf(fileid,'%1.6f  %1.6f\n',[xpoints(i+1:end-1),ypoints(i+1:end-1)]');
fclose(fileid);