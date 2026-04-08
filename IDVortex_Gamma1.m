function varargout = IDVortex_Gamma1(x,y,Vx,Vy)%,U_Vx, U_Vy)
    % Gamma 1 Vortex finder
    %Returns coordinates  [x y] of the peak gamma 1 in the region 
    %If nargout is 1:
    %vortexCoords = IDVortex_Gamma1(x,y,Vx,Vy)
    
    %If nargout is 2:
    %[vortexCoords, G1] = IDVortex_Gamma1(x,y,Vx,Vy)
    
    if nargout>2
        error('Too many output arguments (max is 2)');
    end
    
    
    n=5; %n defines the region size
    s=(n-1)/2; %s defines the shift to left and right
    G1=zeros(size(Vx)); %Store Gamma 1 in G1
    dx=-s:s;
    dy=-s:s;
    [PMx, PMy]=meshgrid(dx,dy);
    PMmag=sqrt(PMx.^2+PMy.^2);
    Vmag=sqrt(Vx.^2+Vy.^2);
    
    Mask=Vmag==0;
    se = strel('disk',floor(s*1.5));
    Mask=imerode(~Mask,se);
    Mask(:,1:2) = 0;
    Mask(:,end-1:end) = 0;
    Mask(1:2,:) = 0;
    Mask(end-1:end,:) = 0;
    
    
    %%
    %==================METHOD 1 - 4 NESTED FOR LOOPS====================
% tic
%     for i=1:size(Vx,1)
%         for j=1:size(Vx,2)
%             g1=zeros(size(PMx));
%             for a=1:length(dx)
%                 for b=1:length(dy)
%                     c=i+a-s-1;
%                     d=j+b-s-1;
%                     
%                     inSquare=(c>=1 & c<=size(Vx,1)) & (d>=1 & d<=size(Vx,2));
%                     
%                     if inSquare
%                         PM=[PMx(a,b); PMy(a,b); 0];
%                         U=[Vx(c,d); Vy(c,d); 0];
%                         
%                         
%                         
% %                         g1(a,b)=(PMx(a,b)*Vy(c,d)-PMy(a,b)*Vx(c,d))/(PMmag(a,b)*Vmag(c,d));
%                         
%                         gector=cross(PM,U)/(norm(PM)*norm(U));      
%                         g1(a,b)=gector(3);
%                         
%                     else
%                         g1(a,b)=nan;
%                     end
%                 end
%             end
%             if all(isnan(g1))
%                 G1(i,j)=nan;
%             else
%                 %G1(i,j)=nanmean(g1(:));
%                 G1(i,j)=nansum(g1(:))/numel(PMx); %By dividing by numel(PMx) we have the same effect as the conv2 method below.
%             end
%         end
%     end
%     
%     G1(~Mask)=nan; %Msks out the effect of borders (quite dramatic)
%     
% toc
    %%
    %==================METHOD 2 - CONVOLUTION====================  
% tic
    
    PMxNorm=PMx./PMmag;
    PMyNorm=PMy./PMmag;
    
    VxNorm=Vx./Vmag;
    VyNorm=Vy./Vmag;
    
    PMxNorm(isnan(PMxNorm))=0;PMyNorm(isnan(PMyNorm))=0;
    VxNorm(isnan(VxNorm))=0;VyNorm(isnan(VyNorm))=0;
    
    %See notes on Matlab exchange
    ConvResult1=conv2(VyNorm,PMxNorm,'same');
    ConvResult2=conv2(VxNorm,PMyNorm,'same');
    
    G1_conv=(ConvResult1-ConvResult2)/numel(PMx);
    G1_conv=abs(G1_conv);
        
    G1_conv(~Mask)=0; %Masks out the effect of borders (quite dramatic)
    
    
% toc
%%
%Benchmark against a 212x438 vector field in an intel i7 4720HQ is:
    %4.11s - 4 nested for loops
    %5.92ms - conv2 method
    
    %Results agreed to numerical roundoff (1e-15)
    
    %%
    %Finds the vortex by finding the max of G1 absolute value and
    %interpolating in a small grid around it

%     figure;
%     imagesc(Mask.*G1_conv)
%     figure;
%     imagesc(Mask)

    [row,col] = find(Mask.*G1_conv==max(Mask.*G1_conv,[],'all'));
    
%     disp(row);
%     disp(col)

    gSize=n;
    colInterval=(col-gSize):(col+gSize);
    rowInterval=(row-gSize):(row+gSize);
    
    xCrop=x(colInterval);
    yCrop=y(rowInterval);
    
    GammaCrop=G1_conv(rowInterval,colInterval);
    
    [Xcrop,Ycrop]=meshgrid(xCrop,yCrop);
    Ggrid=griddedInterpolant(Xcrop',Ycrop',GammaCrop','spline');
    
    xc=xCrop(gSize+1); yc=yCrop(gSize+1);
    gc=Ggrid(xc,yc);
    
    %Solves the optimization problem in a spline grid and returns
    options = optimoptions('fmincon','Display','off');
    vortexCoords = fmincon(@(X)-Ggrid(X(1),X(2)),[xc yc],[0 0],0,[0 0],0,[xCrop(1) yCrop(1)],[xCrop(end) yCrop(end)],[],options);
    
    disp('Vortex Found')

    if nargout==1
        varargout={vortexCoords};
    elseif  nargout==2
        varargout={vortexCoords,G1_conv};
    end
    
    %imagesc(abs(G1_conv)); colormap inferno; drawnow; pause(0.01);%If you
    %want to see how the gamma map look like
end