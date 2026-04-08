clear; clc; close all;

%This code generates user-selected masks for the SSPIV processing.
d=dir('*.im7');

startImage=1; %First image of the set to process
startMotionImage=1; %[Run 5409] 
endImage=700; %[Run 5409] 

skipMask=50; %Number of images to skip between manual masks
cAxisMin=350; cAxisMax=1300;

%% Opens the images in descending order, using the last frame as the polygon basis
filesToOpen=endImage:-skipMask:startMotionImage;
if filesToOpen(end)~=startMotionImage
    filesToOpen(end)=startMotionImage;
end

for i=1:length(filesToOpen)
    I=readimx(fullfile(d(filesToOpen(i)).folder,d(filesToOpen(i)).name));
    %One of the versions of readimx uses I.Frames.Components.Planes
%     Camera1=I.Frames{1,1}.Components{1,1}.Planes{1}; Camera1=Camera1';
%     Camera2=I.Frames{3,1}.Components{1,1}.Planes{1}; Camera2=Camera2';
    
    %There's another version of readimx which uses the "I.Data" to store
    %the camera frame instead:
    Camera1=I.Data(:,1:I.Ny)';
    Camera2=I.Data(:,(I.Ny*2+1):(I.Ny*3))';

    if i==1
        imagesc(Camera1); caxis([cAxisMin cAxisMax]); title('Select Mask for CAMERA 1')
        [M, xM1, yM1]=roipoly();
        close all;
        
        imagesc(Camera2); caxis([cAxisMin cAxisMax]); title('Select Mask for CAMERA 2')
        [M2, xM2, yM2]=roipoly();
        close all;

        Mask1(i).x=xM1;Mask1(i).y=yM1;
        Mask2(i).x=xM2;Mask2(i).y=yM2;
    else
        imagesc(Camera1); caxis([cAxisMin cAxisMax]); title(['Review Mask for CAMERA 1 for Frame ' num2str(filesToOpen(i)) ', then ENTER']);
        h = drawpolygon('Position',[Mask1(i-1).x Mask1(i-1).y]);
        pause
        Mask1(i).x=h.Position(:,1);Mask1(i).y=h.Position(:,2);
        close all;

        imagesc(Camera2); caxis([cAxisMin cAxisMax]); title(['Review Mask for CAMERA 2 for Frame ' num2str(filesToOpen(i)) ', then ENTER']);
        h = drawpolygon('Position',[Mask2(i-1).x Mask2(i-1).y]);
        pause
        Mask2(i).x=h.Position(:,1);Mask2(i).y=h.Position(:,2);
        close all;        
    end
end

%% Shows the polygon shrinking around the model as an animation for confirmation purposes
skipImg=1;
saveImages=1; %Whether to save the images or not

if saveImages    
    Maskdir=fullfile(pwd,'Mask');
    %Creates a mask folder if it doesn't exist
    if ~exist(Maskdir, 'dir')
        mkdir(Maskdir);
    else
        %Deletes and re-creates the mask folder if it already exists to
        %start afresh
        rmdir(Maskdir);
        mkdir(Maskdir);
    end
end

figure('position',[770 418 1075 536])
for i=startImage:skipImg:endImage
    I=readimx(fullfile(d(i).folder,d(i).name));
    %One of the versions of readimx uses I.Frames.Components.Planes
%     Camera1=I.Frames{1,1}.Components{1,1}.Planes{1}; Camera1=Camera1';
%     Camera2=I.Frames{3,1}.Components{1,1}.Planes{1}; Camera2=Camera2';
    
    %There's another version of readimx which uses the "I.Data" to store
    %the camera frame instead:
    Camera1=I.Data(:,1:I.Ny)';
    Camera2=I.Data(:,(I.Ny*2+1):(I.Ny*3))';

    subplot(1,2,1)
    imagesc(Camera1); daspect([1 1 1]); caxis([cAxisMin cAxisMax]); title('Camera 1')
    subplot(1,2,2)
    imagesc(Camera2); daspect([1 1 1]); caxis([cAxisMin cAxisMax]); title('Camera 2')

    if i<startMotionImage
        %Does not interpolate; use the mask where the motion starts
        subplot(1,2,1)
        h1 = drawpolygon('Position',[Mask1(end).x Mask1(end).y]);
        subplot(1,2,2)
        h2 = drawpolygon('Position',[Mask2(end).x Mask2(end).y]);

        %Saves the mask if flag requests
        if saveImages
            M1=roipoly(Camera1,Mask1(end).x,Mask1(end).y);
            M2=roipoly(Camera2,Mask2(end).x,Mask2(end).y);
            FullMask=[M1; M1; M2; M2];
            imwrite(~FullMask,fullfile(Maskdir, ['B' num2str(i,'%05.f') '.jpg']),'jpg');
        end
    else
        %Interpolates the mask from the manual sets
        AllM1x=[Mask1(:).x]'; AllM1y=[Mask1(:).y]';
        AllM2x=[Mask2(:).x]'; AllM2y=[Mask2(:).y]';

        ThisM1x = interp1(filesToOpen,AllM1x,i); ThisM1y = interp1(filesToOpen,AllM1y,i);
        ThisM2x = interp1(filesToOpen,AllM2x,i); ThisM2y = interp1(filesToOpen,AllM2y,i);

        subplot(1,2,1)
        h1 = drawpolygon('Position',[ThisM1x' ThisM1y']);
        subplot(1,2,2)
        h2 = drawpolygon('Position',[ThisM2x' ThisM2y']);

        %Saves the mask if flag requests
        if saveImages
            M1=roipoly(Camera1,ThisM1x,ThisM1y);
            M2=roipoly(Camera2,ThisM2x,ThisM2y);
            FullMask=[M1; M1; M2; M2];
            imwrite(~FullMask,fullfile(Maskdir, ['B' num2str(i,'%05.f') '.jpg']),'jpg');
        end
    end

    drawnow;
    pause(0.1);

    %Removes the polygons from the figure to prevent drawing more polygons
    delete(h1);
    delete(h2);    
end







