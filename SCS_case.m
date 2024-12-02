% Main purpose: use slepian functions to calculate the sea level anomaly in
%   a regional ocean.
% This code is based on the Free Software from the Simons Laboratories
%   (https://geoweb.princeton.edu/people/simons/software.html).
%   with some changes for oceanic application (e.g., GAD, GIA, IB correction).
% It is strongly recommended to familiarize yourself with the original code
%   before using this add-on code.
%
% Required software:
%   slepian_alpha (by csdms-contrib)
%   slepian_bravo (by csdms-contrib)
%   slepian_delta (by csdms-contrib)
%   slepian_zero (actually only guyotphysics.m) (by csdms-contrib)
%   m_map (by https://www.eoas.ubc.ca/~rich/map.html)

% Main Reference:
% Ref1: Harig, C., & Simons, F. J. (2012). Mapping GreenlandAs mass loss in space and time. Proceedings of the National Academy of Sciences, 109(49), 19934-19937.
% Ref2: Zhongtian, Ma et al., A Novel Slepian Approach for Determining Mass-term Sea Level from GRACE over the
%           South China Sea

% Last modified by zhongtian.ma@connect.polyu.hk, 2/12/2024

clc;clear;close all

setenv('IFILES','./');
getenv('IFILES')
ddir1=fullfile(getenv('IFILES'),'Results');

% where the required software is
addpath(fullfile(getenv('IFILES'),'slepian_abd-master_untouched'))
addpath(fullfile(getenv('IFILES'),'slepian_abd-master_untouched','REGIONS'))
addpath(fullfile(getenv('IFILES'),'slepian_abd-master_untouched','ICEREGIONS'));

% where this add-on code is
addpath(fullfile(getenv('IFILES'),'updated_code'))
addpath(fullfile(getenv('IFILES'),'updated_code','REGIONS'));

%% The main procedures below follow the README.MD files in 'slepian_delta' software

% changeable coefficients
Dataproduct={'CSR','RL06',60,'ocean','Pelt17'}; % dataproduct
phi=0; theta=0; omega=0; % coefficients of tapers (see grace2slept.m)
Lwindow=Dataproduct{3}; % maximum degree of useful gravity coefficient
data_period=[2005,2015]; % the start year and end year of the data set

Radius=500; % radius of additional Gaussian smooth (unit: km) 
TH='SCSpTH'; Area=TH; % study area (take SCS as example)

% The degree of the buffer zones
%  positive values (for land application) make buffer zone outside the study area
%  negative values (for ocean application) make buffer zone inside the study area
XY_buffer=-1; % unit: degree

% study boundary
c11cmn=[95.5 29.5 124.5 -4.5];

if XY_buffer>=0
    buffer_str=num2str(XY_buffer);
    buffer_str(buffer_str=='.')='p';
else
    buffer_str=['neg' num2str(abs(XY_buffer))];
    buffer_str(buffer_str=='.')='p';
end

FIG_Attach=sprintf('%s_%s_%s_%s_%s_%s',...
    Dataproduct{1},Dataproduct{2}(3:4),Area,num2str(Lwindow),buffer_str,num2str(Radius));

radius=6370000; % be consistent with code

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot for the study area (blue line) and that with buffer 
% zone (yellow shaded area)

f0=figure;
XY=eval(sprintf('%s(%i,%f)',TH,10,XY_buffer));
BasinArea=spharea(XY)*4*pi*radius^2;
geoshow(XY(:,2), XY(:,1), 'DisplayType', 'polygon', ...
    'FaceColor', 'yellow')

hold on

XY_shp=eval(sprintf('%s(%i,%f)',TH,10,0));
BasinArea0=spharea(XY_shp)*4*pi*radius^2;
plot(XY_shp(:,1),XY_shp(:,2),'b','linewidth',3);

title('Study area');

tif_name0=['Fig0_',FIG_Attach,'.tif'];
print(f0,'-dtiff','-r125',fullfile(getenv('IFILES'),'Figure',tif_name0));
%%%%%%%%%%%%%%%%%%%%%%%%%%%

save(fullfile(getenv('IFILES'),'Results',['XY_',FIG_Attach,'.mat']),'XY')

%%  1. Set up your data and directory structure, and then with GRACE2PLMT
% read in the GRACE data files into a matrix for use in Matlab. This
% function does corrections for C2,0 from the SLR values of Cheng and
% Tapley, [2004] and degree 1 values from Swenson, [2008].

% Please be aware to check the release of your CSR data and revise the
%   grace2plmt_up accordingly.
% Please also note that you may need to download the corresponding degree 1 gravity
%   coefficient file and save it into ./GRACE/Degree1/ 
%   (e.g., TN-13_GEOC_CSR_RL0602.txt from https://podaac.jpl.nasa.gov/gravity/grace-documentation#TechnicalNotes) 

% Updated function
% Additional input:
%   land_or_ocean   Whether the data was used on land or ocean. This will
%                   change the process on degree 1
%   GIA             Whether you want to remove GIA (e.g., 'Pelt17' for RL06) 
%                   or not [default].

land_or_ocean='ocean';
GIA='Pelt17';
[potcoffs,thedates]=grace2plmt_up('CSR','RL06',60,'SD',0,land_or_ocean,GIA);

%% 2. Decide on your choice of basis, depending on your region of interest
% and the bandwidth you want. Using GRACE2SLEPT, project the results of
% GRACE2PLMT into your chosen basis. We recommend that your basis is
% chosen based on a set of synthetic experiments which estimate the
% leakage/recovery tradeoffs.

% Updated function suitable for ocean application
[slepcoffs,calerrors,thedates,TH,G,CC,V]=grace2slept_up(Dataproduct,...
    TH,XY_buffer,Lwindow,phi,theta,omega,[],'SD',1);

% selection of truncation number of Spherical Slepian function

% Total N used (concentration radio > 0.1)
N=sum(V>0.1);
% First threhold N1 (Shannon Number)
N1=round((Lwindow+1)^2*spharea(XY)); 
% Second threhold N2 (concentration radio > 0.1)
N2=sum(V>0.1);

FIG_Attach=[FIG_Attach '_N' num2str(N1) 't' num2str(N2)]

%% Additional procedure for finding leakaging date
% manually adjustment is somehow needed for specific dataset

num_year_str=char(Dataproduct(1));
num_month=(data_period(2)-data_period(1)+1)*12;

given_mon=(str2num(datestr(thedates,'yyyy'))-data_period(1))*12+ ...
    str2num(datestr(thedates,'mm'));

% manually adjustment 
if data_period(1)==2003
    if length(given_mon)>=104
        given_mon(104)=107; % Only for JPL
    end
    if length(given_mon)>=138
        given_mon(138)=149; % for CSR, GFZ, JPL 05
    end
elseif data_period(1)==2005
    if length(given_mon)>=81
        given_mon(81)=83; % Only for JPL
    end
    if length(given_mon)>=115
        given_mon(115)=125; % for CSR, GFZ, JPL 05
    end
else

end

% check if all leakage months are calculated without double-counted
if length(unique(given_mon))~=length(given_mon')
    uu=unique(given_mon);
    repeat_num=length(given_mon')-length(unique(given_mon));

    compare=[unique(given_mon); ones(repeat_num,1)]-given_mon;

    error(['Repeat number exists. Repeat number: ' ...
        num2str(repeat_num)])
else
    fprintf('No repeat number');
end

EST_mon=[1:num_month]';

leakage=setdiff(EST_mon,given_mon)';

%fill the time epoch by linear interpolate
thedates_int = interp1(given_mon,thedates,leakage);
ESTthedates = interp1(given_mon,thedates,EST_mon)';

Sleptthedates=thedates;
if any(leakage)
    for i=1:numel(leakage)
        numleak=leakage(i);
        %select the mean time epoch between the adjacent months and put
        %them behind the raw data array
        Sleptthedates(length(given_mon)+i)=round(thedates_int(i));
    end
end

%% 3. Next run SLEPT2RESID to fit a choice of functions (e.g. lines,
% quadratics, etc.) to the Slepian coefficients. Note: if you want to
% remove a model of GIA then you should do that before this step, using
% a function like CORRECT4GIA.

fitwhat=[1 365.25 365.25/2]; %linear trend and an annual and a semiannual harmonic term

[ESTsignal,ESTresid,ftests,extravalues,total,alphavarall,totalparams,...
    totalparamerrors,totalfit,functionintegrals,alphavar]=slept2resid_up(slepcoffs,...
    Sleptthedates,fitwhat,[],[],CC,TH,[N1,N2],Radius);

%% 4. If you want to examine the total mass trend in a region, this
% information is calculated and returned as a result from SLEPT2RESID.
% To summarize, each Slepian coefficient up to the Shannon truncation is
% multiplied by the integral (found using INTEGRATEBASIS) of the corresponding
% function over the region. This total data is then fit with TIMESERIESFIT
% which allows fitting with a custom data variance.

% Mass sea level (MSL)= Ocean Bottom Pressure (OBP) - Inverted Barometer (IB)

given_nmonths=length(thedates);EST_nmonths=length(ESTthedates);

[Cab] = slepresid2cov(ESTresid);

% Make the coefficients with reference to some mean
% If they already are, then this won't matter
sleptdelta = slepcoffs(1:given_nmonths,:) - repmat(mean(slepcoffs(1:given_nmonths,:),1),given_nmonths,1);

% load Inverted Barometer correction
load(fullfile(getenv('IFILES'),'IB','IB0315_global_NCL'));
IB=IB0315_global_NCL/10; % cm

% only use the fit procedure in this functionto get the residual of IB
[IBsignal,IBresid,IB_ftests,IB_extravalues]=slept2resid_up(IB(given_mon+12*(data_period(1)-2003))*10,...
    Sleptthedates,fitwhat); % from cm to mm (values equivalent to kg/m^2)

% filled leakage values for estiamted signals and residuals
if any(leakage)
    ESTsignal_com=zeros(length(ESTthedates),size(ESTsignal,2));
    ESTresid_com=zeros(length(ESTthedates),size(ESTresid,2));
    IBsignal_com=zeros(length(ESTthedates),size(IBsignal,2));
    IBresid_com=zeros(length(ESTthedates),size(IBresid,2));
    for i=1:size(ESTsignal,1)
        ESTsignal_com(given_mon(i),:)=ESTsignal(i,:);
        ESTresid_com(given_mon(i),:)=ESTresid(i,:);
        IBsignal_com(given_mon(i),:)=IBsignal(i,:);
        IBresid_com(given_mon(i),:)=IBresid(i,:);
    end

    % calculated leakage values
    extravalues_resid=zeros(length(leakage),size(ESTresid,2));
    extravalues_IBresid=zeros(length(leakage),size(IBresid,2));
    for j=1:size(ESTresid_com,2)
        ESTresid_leak = interp1(thedates,ESTresid(:,j)',thedates_int,'spline');
        extravalues_resid(:,j)=ESTresid_leak';
    end
    IBresid_leak = interp1(thedates,IBresid(:,1)',thedates_int,'spline');
    extravalues_IBresid(:,1)=IBresid_leak';

    for i=1:size(leakage,2)
        ESTsignal_com(leakage(i),:)=extravalues(i,:);
        %             ESTresid_com(leakage(i),:)=extravalues_resid(i,:);
        ESTresid_com(leakage(i),:)=0; % here we assume zero for leakage residuals
        IBsignal_com(leakage(i),:)=IB_extravalues(i,:);
        %             IBresid_com(leakage(i),:)=extravalues_IBresid(i,:);
        IBresid_com(leakage(i),:)=0; % here we assume zero for leakage residuals
    end
    ESTsignaldelta=ESTsignal_com(1:EST_nmonths,:) - repmat(mean(ESTsignal(1:given_nmonths,:),1),EST_nmonths,1);
    ESTresiddelta=ESTresid_com(1:EST_nmonths,:) - repmat(mean(ESTresid_com(1:given_nmonths,:),1),EST_nmonths,1);
    IBsignaldelta=IBsignal_com(1:EST_nmonths,:) - repmat(mean(IBsignal(1:given_nmonths,:),1),EST_nmonths,1);
    IBresiddelta=IBresid_com(1:EST_nmonths,:) - repmat(mean(IBresid_com(1:given_nmonths,:),1),EST_nmonths,1);

else
    ESTsignaldelta=ESTsignal(1:given_nmonths,:) - repmat(mean(ESTsignal(1:given_nmonths,:),1),given_nmonths,1);
    ESTresiddelta=ESTresid(1:given_nmonths,:) - repmat(mean(ESTresid(1:given_nmonths,:),1),given_nmonths,1);
    IBsignaldelta=IBsignal(1:given_nmonths,:) - repmat(mean(IBsignal(1:given_nmonths,:),1),given_nmonths,1);
    IBresiddelta=IBresid(1:given_nmonths,:) - repmat(mean(IBresid(1:given_nmonths,:),1),given_nmonths,1);
end
ESTfilldelta=ESTsignaldelta+ESTresiddelta;
IBfilldelta=IBsignaldelta+IBresiddelta;

%%%%%%%%%%%%%%%%%%%%%%%%%
% Total mass without considering IB (i.e., OBP)

% Multiple the Slepian coefficient (signal) by the integral (found using INTEGRATEBASIS) 
% of the corresponding function over the region

% raw OBP time series without leakage filling
[total_CC]=integral_fit(N,functionintegrals,sleptdelta);

% signal OBP time series with leakage filling
[ESTsigtotal_CC]=integral_fit(N,functionintegrals,ESTsignaldelta);

OBP_ESTsigtotal=sum(ESTsigtotal_CC);

% residual OBP time series with leakage filling
[ESTrestotal_CC]=integral_fit(N,functionintegrals,ESTresiddelta);

OBP_ESTrestotal=sum(ESTrestotal_CC);

% total OBP time series with leakage filling
OBP_ESTreconsttotal=OBP_ESTsigtotal+OBP_ESTrestotal;
%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%
% Total mass considering IB (i.e., MSL)

% substitude the N+1 slepian function as IB, because we only used the
% truncated slepian functions up to N, the N+1 was used to correct the IB
% mannualy

% consider IB in error propagation
ESTresid_IB=ESTresid;
ESTresid_IB(:,N+1)=IBresid;
[Cab_IB] = slepresid2cov(ESTresid_IB);

% functionintegrals_IB=[functionintegrals -1]; % -1 represents subtracting IB
% correction on 2/12/2024 (the unit of IB is wrong in previous version)
functionintegrals_IB=[functionintegrals -1/10^3*BasinArea*10^3/10^3/10^9]; % represents subtracting IB (change unit from mm to Gt)

alphavarall_IB=functionintegrals_IB*Cab_IB(1:N+1,1:N+1)*functionintegrals_IB';
alphavarall_EST=functionintegrals*Cab(1:N,1:N)*functionintegrals';

sleptdelta_IB=sleptdelta;
sleptdelta_IB(:,N+1)=IB(given_mon+12*(data_period(1)-2003))*10;

% Multiple the Slepian coefficient (signal) by the integral (found using INTEGRATEBASIS) 
% of the corresponding function over the region

% total mass considering IB (i.e., MSL)
[MSL_total_CC]=integral_fit(N+1,functionintegrals_IB,sleptdelta_IB);

MSL_total=sum(MSL_total_CC);

[MSL_fit,MSL_delta,MSL_totalparams,MSL_paramerrors] = timeseriesfit([thedates' MSL_total'],...
    alphavarall_IB,1,1);
% Make a matrix for the line, and 95% confidence in the fit
MSL_totalfit = [thedates' MSL_fit MSL_delta];
MSL_totalparamerrors = MSL_paramerrors*365.25;
%%%%%%%%%%%%%%%%%%%%%%%%%

% Multiple the Slepian coefficient (signal) by the integral (found using INTEGRATEBASIS) 
% of the corresponding function over the region
% Note that here Cab could not be used actually.

% signal MSL time series with leakage filling
MSL_ESTsigtotal=sum(ESTsigtotal_CC)-IBsignaldelta';

% residual MSL time series with leakage filling
MSL_ESTrestotal=sum(ESTrestotal_CC)-IBresiddelta';

% total MSL time series with leakage filling
MSL_ESTreconsttotal=MSL_ESTsigtotal+MSL_ESTrestotal;

%% 5. If you want the total map of mass change, multiply the difference 
% between estimated signal coefficients by the corresponding Slepian 
% function, and add them up. Remember the appropriate units. This sum can
% then be expanded to space using PLM2XYZ By limiting the set of 
% coefficients you use, you can instead make this for any date span, 
% such as for several years or a single year.

[r_example,lon,lat,Plm,degres]=plm2xyz(CC{1},1,c11cmn,60);

% J could be one or two
% J=[N]; % example 1
J=[N1,N2]; % example 2

% time wasting procedure, please wait
if length(J)==1
    if Radius>0
        for j=1:J
            CC_smo=plm_Gausssmooth(CC{j},Radius)
            [r_record(j,:,:),lon,lat,Plm,degres]=plm2xyz(CC_smo,1,c11cmn,Lwindow);
        end

    else
        for j=1:J
            [r_record(j,:,:),lon,lat,Plm,degres]=plm2xyz(CC{j},1,c11cmn,Lwindow);
        end
    end

elseif length(J)==2
    for j=1:J(1)
        [r_record(j,:,:),lon,lat,Plm,degres]=plm2xyz(CC{j},1,c11cmn,60);
    end
    for j=J(1)+1:J(2)
        CC_smo=plm_Gausssmooth(CC{j},Radius);
        [r_record(j,:,:),lon,lat,Plm,degres]=plm2xyz(CC_smo,1,c11cmn,60);
    end

else
    eror('Wong input J.')
end

MSL_ESTsignal=zeros(EST_nmonths,size(r_example,1),size(r_example,2));
MSL_ESTresidual=zeros(EST_nmonths,size(r_example,1),size(r_example,2));
OBP_ESTsignal=zeros(EST_nmonths,size(r_example,1),size(r_example,2));
OBP_ESTresidual=zeros(EST_nmonths,size(r_example,1),size(r_example,2));
for i=1:EST_nmonths
    sp_msl_signal=0;sp_msl_residual=0;
    for j=1:N

        r=squeeze(r_record(j,:,:));

        sp_msl_signal=sp_msl_signal+r*ESTsignaldelta(i,j)' ... % (you can manually adjust this period)
            /1000*1000; % change from kg/m^2 to mm (/ 1000 kg/m*3 * 1000)
        sp_msl_residual=sp_msl_residual+r*ESTresiddelta(i,j)' ... % (you can manually adjust this period)
            /1000*1000; % change from kg/m^2 to mm (/ 1000 kg/m*3 * 1000)
    end
    % substract IB
    MSL_ESTsignal(i,:,:)=sp_msl_signal-IBsignaldelta(i)';
    MSL_ESTresidual(i,:,:)=sp_msl_residual-IBresiddelta(i)';
    OBP_ESTsignal(i,:,:)=sp_msl_signal;
    OBP_ESTresidual(i,:,:)=sp_msl_residual;
end

MSL_ESTreconst=MSL_ESTsignal+MSL_ESTresidual;
OBP_ESTreconst=OBP_ESTsignal+OBP_ESTresidual;

% signal + residual
MSL_slept=zeros(given_nmonths,size(r_example,1),size(r_example,2));
OBP_slept=zeros(given_nmonths,size(r_example,1),size(r_example,2));
for i=1:given_nmonths
    sp_msl=0;
    for j=1:N

        r=squeeze(r_record(j,:,:));

        sp_msl=sp_msl+r*sleptdelta(i,j)' ... % (you can manually adjust this period)
            /1000*1000; % change from kg/m^2 to mm (/ 1000 kg/m*3 * 1000)
    end
    % add IB
    MSL_slept(i,:,:)=sp_msl-IB(given_mon(i)+12*(data_period(1)-2003))'*10; % from cm to mm
    OBP_slept(i,:,:)=sp_msl;
end

%% (Case Figure 1) Temporal mass sea level 
% change mass to EWH (transform Gt to kg/m^2 to cm: Gt * (10^9 * 10^3 / BasinArea) * (/10^3 *10^2))
Total_MSL_EST=MSL_total*10^9*10^3/BasinArea/10;
Total_MSL_EST_fit=MSL_totalfit*10^9*10^3/BasinArea/10;
Total_MSL_EST_slope=MSL_totalparams(2)*10^9*10^3/BasinArea/10*365.25;
Total_MSL_EST_signal=MSL_ESTsigtotal*10^9*10^3/BasinArea/10;
Total_MSL_EST_resid=MSL_ESTrestotal*10^9*10^3/BasinArea/10;

figure
ttpp=floor(given_mon/12)+data_period(1)+mod(given_mon,12)/12-1/24;
ttff=floor(EST_mon/12)+data_period(1)+mod(EST_mon,12)/12-1/24;

p1=plot(ttpp,Total_MSL_EST,'color','black','linewidth',2);
hold on
p2=plot(ttpp,Total_MSL_EST_fit(:,2),'r','linewidth',1);
hold on
p3=plot(ttpp,(Total_MSL_EST_fit(:,2)+Total_MSL_EST_fit(:,3)),'m','linestyle','--','linewidth',0.7);
hold on
plot(ttpp,(Total_MSL_EST_fit(:,2)-Total_MSL_EST_fit(:,3)),'m','linestyle','--','linewidth',0.7);
title('GRACE SSF','FontWeight','bold')
hold on
p4=plot(ttff,Total_MSL_EST_signal,'linewidth',2);
ylabel('Mass sea level (cm)','FontWeight','bold')

xlim(data_period)
%% (Case Figure 2) coefficients eigenfunctions (CC)

[lonlon,latlat]=meshgrid(lon,lat);

f5=figure;

t=tiledlayout('flow'); % nexttile method

for i=1:N

    % choose nexttile or subplot (based on your matlab version containing tiledlayout or not)
    nexttile % nexttile method

    r=squeeze(r_record(i,:,:));

    r_column=reshape(r,[],1);

    r_res=rescale(r_column,-1,1);

    r_res=reshape(r_res,size(r));

    r_ewh=r*( ESTsignaldelta(end-11,i) - ESTsignaldelta(1,i) )/(ceil(length(EST_mon)/12) - 1) ... % (you can manually adjust this period)
        /1000*100; % change from kg/m^2 to cm (/ 1000 kg/m*3 * 100)

    m_proj('mercator','long',[c11cmn(1)-1.5, c11cmn(3)+0.5],'lat',[c11cmn(4)-0.5, c11cmn(2)+1.5]);
    m_pcolor(lonlon,latlat,r_res);
    % m_coast('patch',[1 .85 .7]);
    m_coast
    m_grid('box','fancy','tickdir','in');
    m_line(XY(:,1),XY(:,2),'color','k','linewidth',1);

    title(['CC(' num2str(i) ') \lambda=' num2str(V(i),'%.4f')],'color','red','fontsize',12)

end

%%% nexttile method
t.TileSpacing='Compact';
t.Padding='Compact';
%     colormap(jet)
colormap(flipud (jet))

cb = colorbar;
cb.Layout.Tile = 'south';
title(cb,'magnitude','fontsize',14);
%%%


% You'd better adjust this manually
set (gcf,'Position',[100,100,1000,800])

%% (Case Figure 3) PLOT for spatial total, siganl and residual mass (Fig.2 in Ref1)
% From the first month of the first year to the first month of the last year (default without missing month)

[lon1,lat1]=meshgrid([c11cmn(1)-0.5:c11cmn(3)-0.5],[c11cmn(2)+0.5:-1:c11cmn(4)+0.5]);

%     [r_example,lon,lat,Plm,degres]=plm2xyz(CC{1},1,c11cmn,60);

str_begdate=[num2str(floor(EST_mon(1)/12)+data_period(1)) '.' num2str(mod(EST_mon(1)-1,12)+1)];
str_enddate=[num2str(floor(EST_mon(end-11)/12)+data_period(1)) '.' num2str(mod(EST_mon(end-11)-1,12)+1)];

f2=figure;

subplot(1,3,1)
m_proj('mercator','long',[c11cmn(1)-1.5, c11cmn(3)+0.5],'lat',[c11cmn(4)-0.5, c11cmn(2)+1.5]);
m_pcolor(lon1,lat1, squeeze(MSL_slept(111,:,:) - MSL_slept(1,:,:) )/10 );
% m_coast('patch',[1 .85 .7]);
m_coast
m_grid('box','fancy','tickdir','in');
m_line(XY(:,1),XY(:,2),'color','k','linewidth',1);

title({['total mass variance between ' str_begdate '-' str_enddate],['Int=' num2str( ( MSL_ESTreconsttotal(end-11) - MSL_ESTreconsttotal(1) ), '%.2f') ' Gt']},'color','red')

colormap(flipud (jet))
hc=colorbar('southoutside');
caxis([-100 100]);
ylabel(hc, 'surface density (cm)');

subplot(1,3,2)
m_proj('mercator','long',[c11cmn(1)-1.5, c11cmn(3)+0.5],'lat',[c11cmn(4)-0.5, c11cmn(2)+1.5]);
m_pcolor(lon1,lat1, squeeze(MSL_ESTsignal(121,:,:) - MSL_ESTsignal(1,:,:) )/10 );
% m_coast('patch',[1 .85 .7]);
m_coast
m_grid('box','fancy','tickdir','in');
m_line(XY(:,1),XY(:,2),'color','k','linewidth',1);

title({['signal mass variance between ' str_begdate '-' str_enddate],['Int=' num2str( ( MSL_ESTsigtotal(end-11) - MSL_ESTsigtotal(1) ), '%.2f') ' Gt']},'color','red')

colormap(flipud (jet))
hc=colorbar('southoutside');
caxis([-20 20]);
ylabel(hc, 'surface density (cm)');

subplot(1,3,3)
m_proj('mercator','long',[c11cmn(1)-1.5, c11cmn(3)+0.5],'lat',[c11cmn(4)-0.5, c11cmn(2)+1.5]);
m_pcolor(lon1,lat1, squeeze(MSL_ESTresidual(121,:,:) - MSL_ESTresidual(1,:,:) )/10 );
% m_coast('patch',[1 .85 .7]);
m_coast
m_grid('box','fancy','tickdir','in');
m_line(XY(:,1),XY(:,2),'color','k','linewidth',1);

title({['Residual mass variance between ' str_begdate '-' str_enddate],['Int=' num2str( ( MSL_ESTrestotal(end-11) - MSL_ESTrestotal(1) ), '%.2f') ' Gt']},'color','red')

colormap(flipud (jet))
hc=colorbar('southoutside');
caxis([-100 100]);
ylabel(hc, 'surface density (cm)');

set (gcf,'Position',[100,100,1200,400])
