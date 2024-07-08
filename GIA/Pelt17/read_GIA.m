% Main purpose: Test the unit of GIA is yearly data or not

clc;clear;

% change the 'IFILES' for your route
setenv('IFILES','./');
getenv('IFILES')

addpath(fullfile(getenv('IFILES'),'slepian_abd-master'))

% doi: 10.1002/2016JB013844
fid=fopen('ICE-6G_D_VM5a_O512_GRACE.txt');
C=textscan(fid,'%f%f%f%f');

lmcosi=cell2mat(C);

a=fralmanac('a_EGM96','Earth');

% Change from geopotential coefficients to the rate of change of 
% surface mass density (kg/m*2) per year as required by 'correct4gia.m'
% Note: the value of surface density (kg/m*2) is also equivalent to EWH (mm)
lmcosiM=plm2pot([lmcosi(:,1:2) lmcosi(:,3:4)*a],[],[],[],4); 

% you can also check wther the spatial pattern of GIA is consistent with
% previous work
% note that m_map package is needed for plot

c11cmn=[0 90 360 -90];

[r,lon,lat,Plm,degres]=plm2xyz(lmcosiM,1,[0 90 360 -90]);

figure

[lon1,lat1]=meshgrid([c11cmn(1):c11cmn(3)],[c11cmn(2):-1:c11cmn(4)]);

m_proj('robinson','long',[c11cmn(1), c11cmn(3)],'lat',[c11cmn(4), c11cmn(2)]);
m_pcolor(lon1,lat1,r); 
% m_coast('patch',[1 .85 .7]);
m_coast
m_grid('box','fancy','tickdir','in');
% m_line(XY(:,1),XY(:,2),'color','k','linewidth',1);

colormap(flipud (jet))
hc=colorbar('southoutside');
caxis([-25 25]);
ylabel(hc, 'GIA (mm)');

% Results:
% By comparing the figure with Figure 2 of B.D. Vishwakarma et al. (2022), 
% we confirm that the data is values per year
% B.D. Vishwakarma et al., 2022, https://doi.org/10.1093/gji/ggab464

% then save the the value of surface density (kg/m*2)

save('Pelt17_SD',"lmcosiM");
