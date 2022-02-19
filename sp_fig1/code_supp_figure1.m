clear all;clc
load('mycolor.mat');
%%%%%%%%%%%%%%%%  For Doppler shift of topographic roughness
filename=strcat('S1A_IW_OCN__2SDV_20211009T231409_20211009T231434_040049_04BDB9_A87E.nc');
 lonimg =ncread(filename,'lon');
 latimg =ncread(filename,'lat');
lonimg=(rot90(lonimg,1));
latimg=(rot90(latimg,1));
lon1=mean(lonimg)';
lat1=mean(latimg,2);
RadVel =ncread(filename,'vv_001_owiRadVel');
RadVel=fliplr(rot90(RadVel,3));
 h1=fspecial('average',[3,3]);
 RadVel=imfilter(RadVel,h1);
Fdca=RadVel/(0.056*3.1416);
Fdca(:,any(isnan(Fdca))) = [];
      

y=Fdca(120,:);y(find(isnan(y)))=[];
x=Fdca(:,200);y(find(isnan(y)))=[];



figure('Color',[1 1 1]);
h=imagesc(lon1,lat1,Fdca);set(gca,'YDir','normal');
i=colorbar;
 colormap(rwb2);
set(h,'alphadata',~isnan(Fdca));
 set(gca,'XTickLabel',{'79.5W','79W','78.5W','78W','77.5W','77W'}) 
 set(gca,'YTickLabel',{'33.4N','33.4N','33.6N','33.8N','34N','34.2N','34.4N','34.6N','34.8N','35N'}) 
box on;
hold on
plot([-79.96,-76.82],[34.8,34.8],'k--','linewidth',5);
hold on
plot([-77.5,-77.5],[33.18,35.14],'r--','linewidth',5);
ylabel(colorbar,'Doppler centroid anomaly (Hz)');
title('(a)','position',[-79.85,35.00],'FontSize',16);


figure('Color',[1 1 1]);
plot(x,'r');
hold on
plot(y,'k');
xlabel('Grid');ylabel('Dca (Hz)');
set(gca,'XLim',[0 255]); 
title('(b)','position',[10,0.5],'FontSize',16);

 
 
 
 
 
 
 
 
 
 
 
 
 