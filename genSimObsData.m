function [obs,rate_bg,phi_tes,theta_tes,psi_tes,simobs]=...
    genSimObsData(varargin)
%GENSIMOBSDATA Generate simulated observation data

% Satellite orbital demo (pointing vector and velocity vector)
% orbital sampling interval. default: 3.9s

global config;loadConfig
global axfontsize
global fontname
global tlfontsize
global gpsf

filesuffix=datetimefilename;

tau=3.9;
modelfile='integral7yr_data.csv';
modelname='INTEGRAL 7-year survey';
marginsz=2.9*pi/180;
fontname='Helvetica';
axfontsize=20;
tlfontsize=axfontsize+4;
fwaxfontsize=14;
fwtlfontsize=fwaxfontsize+2;
bgflux=0.01;

interactive=false;

%% Parsing user input
if isempty(varargin)
    interactive=true;
    userinput=inputdlg({...
        'Sampling interval of orbital simulation, in seconds:',...
        'Model filename:',...
        'Model name:',...
        'EPS fontname:',...
        'EPS axes fontsize:',...
        'EPS title fontsize:',...
        'EPS axes fontsize for full-line-width figures:',...
        'EPS title fontsize for full-line-width figures:',...
        'Number of pixels along each side of a tessella:',...
        'Margin size, in degree:',...
        'Tessella size, in degree:',...
        'Background flux, in mCrabs:'...
        },'User input',1,{...
        num2str(3.9,3),...
        'integral7yr_data.csv',...
        'INTEGRAL 7-year survey',...
        'Helvetica',...
        '20',...
        '24',...
        '14',...
        '16',...
        int2str(config.SAMPLING_FREQ),...
        '2.9',...
        num2str(config.WINDOW_SIZE*180/pi),...
        '0.01'...
        });
    if isempty(userinput)
        errordlg('Canceled by user.')
        error('Canceled by user.')
    end

    if ~isempty(userinput{1}),  tau=str2double(userinput{1});end
    if ~isempty(userinput{2}),  modelfile=userinput{2};end
    if ~isempty(userinput{3}),  modelname=userinput{3};end
    if ~isempty(userinput{4}),  fontname=userinput{4};end
    if ~isempty(userinput{5}),  axfontsize=str2double(userinput{5});end
    if ~isempty(userinput{6}),  tlfontsize=str2double(userinput{6});end
    if ~isempty(userinput{7}),  fwaxfontsize=str2double(userinput{7});end
    if ~isempty(userinput{8}),  fwtlfontsize=str2double(userinput{8});end
    if ~isempty(userinput{9}),  config.SAMPLING_FREQ=str2double(userinput{9});end
    if ~isempty(userinput{10}), marginsz=str2double(userinput{10})*pi/180;end
    if ~isempty(userinput{11}), config.WINDOW_SIZE=str2double(userinput{11})*pi/180;end
    if ~isempty(userinput{12}), bgflux=str2double(userinput{12});end
else
    k=1;
    while k<length(varargin)
        switch varargin{k}
            case 'modelfile'
                modelfile=varargin{k+1};
                k=k+2;
            case 'modelname'
                modelname=varargin{k+1};
                k=k+2;
            case 'bgflux'
                bgflux=varargin{k+1};
                k=k+2;
            case 'exptime'
                tau=varargin{k+1};
                k=k+2;
            case 'ra'
                phi_sel=varargin{k+1}*pi/180;
                k=k+2;
            case 'dec'
                theta_sel=varargin{k+1}*pi/180;
                k=k+2;
            case 'margin'
                marginsz=varargin{k+1}*pi/180;
                k=k+2;
            otherwise
                error('Unsupported parameter.')
        end
    end
end

%% Orbital simulation
if exist('./orbitsim.mat','file')==2
    orbitsim=load('./orbitsim.mat');
    L=length(orbitsim.clock);
    clock=0:tau:orbitsim.clock(L);
    pt_sate=zeros(3,length(clock));
    pt_sate(1,:)=interp1(orbitsim.clock,orbitsim.pt_sate(1,:),clock);
    pt_sate(2,:)=interp1(orbitsim.clock,orbitsim.pt_sate(2,:),clock);
    pt_sate(3,:)=interp1(orbitsim.clock,orbitsim.pt_sate(3,:),clock);
    clock=clock(1:length(clock)-1);
    az_sate=zeros(3,length(clock));
    az_sate(1,:)=interp1(orbitsim.clock(1:L-1),orbitsim.az_sate(1,:),clock);
    az_sate(2,:)=interp1(orbitsim.clock(1:L-1),orbitsim.az_sate(2,:),clock);
    az_sate(3,:)=interp1(orbitsim.clock(1:L-1),orbitsim.az_sate(3,:),clock);
    pt_sate=pt_sate./(ones(3,1)*sqrt(sum(pt_sate.^2,1)));
    [phi_sate,theta_sate]=xyz2ptr(pt_sate(1,:),pt_sate(2,:),pt_sate(3,:));
else
    disp(['orbitsim.mat does not exist.',...
        ' Run satedemo.m to simulate satellite orbital data.'])
    [~,~,~,orbitsim]=satedemo(tau,[],'','','orbitsim.mat');
    pt_sate=orbitsim.pt_sate;
    az_sate=orbitsim.az_sate;
    phi_sate=orbitsim.phi_sate;
    theta_sate=orbitsim.theta_sate;
end
L=size(az_sate,2);
psi_sate=getAzimuth(pt_sate(:,1:L),az_sate);


%% Import model & select the tessella to process

% load model from file:
obj_model=dlmread(modelfile,',',1,0);
RA_obj=obj_model(:,1)*pi/180;
RA_obj(RA_obj>pi)=RA_obj(RA_obj>pi)-2*pi;
Dec_obj=obj_model(:,2)*pi/180;
mCrab_obj=obj_model(:,3);

clmt_cm2=300*5;
photon_keV=20;
mCrab2keV=1.43*6.2415e-3;
rate_obj=((mCrab_obj)*mCrab2keV*clmt_cm2/photon_keV);
rate_bg=(bgflux*mCrab2keV*clmt_cm2/photon_keV);


fgallsky=figure;
axallsky=axes('fontname',fontname,'fontsize',fwaxfontsize);hold all
gridmollweide(30,30,'--k')
[xlim,ylim]=mollweideproj(pi,max(theta_sate(:)));
plot([-xlim,xlim],[ylim,ylim],'-.r')
plot([-xlim,xlim],[-ylim,-ylim],'-.r')
xlabel(axallsky,'R.A., rad','fontsize',fwaxfontsize,...
    'fontname',fontname)
ylabel(axallsky,'Dec., rad','fontsize',fwaxfontsize,...
    'fontname',fontname)

title(axallsky,modelname,'fontsize',fwtlfontsize,'fontname',fontname)
set(axallsky,'fontname',fontname,'fontsize',fwaxfontsize,'Box','on')


[x,y]=mollweideproj(RA_obj,Dec_obj);
scatter(x(:),y(:),4,log10(rate_obj(:)),'filled'),hold all
cbax=colorbar;
set(cbax,'fontname',fontname,'fontsize',fwaxfontsize,'Box','on')
set(get(cbax,'YLabel'),'string','pts/s, log to base 10',...
    'fontname',fontname,...
    'fontsize',fwaxfontsize)

axis xy
axis image
axis tight
drawnow

if interactive
    button=questdlg('How would you like to select the tessella?',...
        'Select a tessella','Input the coordinates',...
        'Click on the map','Use default','Use default');
    switch button
        case 'Input the coordinates'
            coordinates=inputdlg({'R.A., in degrees [-80]: ',...
                'Dec., in degrees [0]: '},...
                'Input the coordinates of the center of the tessella');
            if isempty(coordinates{1})
                phi_sel=-80/180*pi;
            else
                phi_sel=str2double(coordinates{1})/180*pi;
            end
            if isempty(coordinates{2})
                theta_sel=0;
            else
                theta_sel=str2double(coordinates{2})/180*pi;
            end
        case 'Click on the map'
            axis([-pi pi -ylim+config.WINDOW_SIZE/2 ylim-config.WINDOW_SIZE/2])
            drawnow
            [xsel,ysel]=ginput(1);
            [phi_sel,theta_sel]=imollweideproj(xsel,ysel);
        case 'Use default'
            phi_sel=-80/180*pi;
            theta_sel=0;
        otherwise
            errordlg('Unknown user input.')
    end
end

if (theta_sel+config.WINDOW_SIZE/2)>max(theta_sate(:)) ||...
        (theta_sel-config.WINDOW_SIZE/2)<min(theta_sate(:))
    errordlg('Selected tessella exceeds simulated all-sky survey range.')
    error('Selected tessella exceeds simulated all-sky survey range.')
end

filesuffix=[num2str(phi_sel*180/pi,3),'d_',...
    num2str(theta_sel*180/pi,3),'d_',...
    filesuffix];

quat=getQuatDetStatus(phi_sel,theta_sel);
showTessellaBorderOnMollweide(quat,[],'-r');
axis xy
axis image
axis tight
drawnow
figure(fgallsky)
axes(axallsky)
if interactive
    epsname=inputdlg('*.eps filename:','Save to EPS file');
else
%     epsname={['allsky_',filesuffix,'.eps']};
    epsname={};
end

if ~isempty(epsname)
    if ~isempty(epsname{1})
        if exist('epswrite.m','file')==2
            epswrite(fgallsky,epsname{1},'BoundingBox','loose')
        else
            print(fgallsky,'-depsc2',epsname{1},'-loose')
        end
    end
end

obj_tes=isPixelOnTessella(quat,[],RA_obj,Dec_obj);
[~,~,~,phi_tes,theta_tes]=genQuadProj(quat,[],config.SAMPLING_FREQ);
cts_tes=list2map(RA_obj(obj_tes),Dec_obj(obj_tes),rate_obj(obj_tes),...
    phi_tes,theta_tes,0);
u=(((1:config.SAMPLING_FREQ)-1)/(config.SAMPLING_FREQ-1)-0.5)*config.WINDOW_SIZE;
screen=ones(size(u'))*(abs(u)<=(config.WINDOW_SIZE*0.5-marginsz));
screen=screen.*screen';
cts_tes=cts_tes.*screen+rate_bg;
gpsf=fspecial('gauss',config.SAMPLING_FREQ,4);
fgtes=figure;
axtes=axes;
im_tes=imconv(cts_tes,gpsf);
imagesc(u,u,log10(im_tes)),colormap('hot')
cbax=colorbar;
xlabel(axtes,['x rad, centered at ',num2str(phi_sel*180/pi,3),'^\circ R.A.'],...
    'fontsize',axfontsize,'fontname',fontname)
ylabel(axtes,['y rad, centered at ',num2str(theta_sel*180/pi,3),'^\circ Dec.'],...
    'fontsize',axfontsize,'fontname',fontname)
titlestr='Model image of selected region';
title(axtes,titlestr,'fontsize',tlfontsize,'fontname',fontname)
set(axtes,'fontname',fontname,'fontsize',axfontsize,'Box','on')
set(cbax,'fontname',fontname,'fontsize',axfontsize,'Box','on')
set(get(cbax,'YLabel'),'string','pts/s, log to base 10',...
    'fontname',fontname,...
    'fontsize',axfontsize)
figure(fgtes)
axes(axtes)
axis xy
axis tight
axis square
drawnow
if interactive
    epsname=inputdlg('*.eps filename:','Save to EPS file');
else
    epsname={['model_',filesuffix,'.eps']};
end

if ~isempty(epsname)
    if ~isempty(epsname{1})
        if exist('epswrite.m','file')==2
            epswrite(fgtes,epsname{1},'BoundingBox','loose')
        else
            print(fgtes,'-depsc2',epsname{1},'-loose')
        end
    end
end


sate_tes=isPixelOnTessella(quat,[],phi_sate,theta_sate);
azInterp=TriScatteredInterp(phi_sate(sate_tes)',theta_sate(sate_tes)',...
    psi_sate(sate_tes)','nearest');
psi_tes=azInterp(phi_tes,theta_tes);
fgaz=figure;
axaz=axes;
imagesc(u,u,psi_tes*180/pi)
cbax=colorbar;
xlabel(axaz,['x rad, centered at ',num2str(phi_sel*180/pi,3),'^\circ R.A.'],...
    'fontsize',axfontsize,'fontname',fontname)
ylabel(axaz,['y rad, centered at ',num2str(theta_sel*180/pi,3),'^\circ Dec.'],...
    'fontsize',axfontsize,'fontname',fontname)
titlestr='Position angle map of selected region';
title(axaz,titlestr,'fontsize',tlfontsize,'fontname',fontname)
set(axaz,'fontname',fontname,'fontsize',axfontsize,'Box','on')
set(cbax,'fontname',fontname,'fontsize',axfontsize,'Box','on')
set(get(cbax,'YLabel'),'string','position angle, in degree','fontname',fontname,...
    'fontsize',axfontsize)
figure(fgaz)
axes(axaz)
axis xy
axis tight
axis square
drawnow
if interactive
    epsname=inputdlg('*.eps filename:','Save to EPS file');
else
    epsname={['posangle_map_',filesuffix,'.eps']};
end
if ~isempty(epsname)
    if ~isempty(epsname{1})
        if exist('epswrite.m','file')==2
            epswrite(fgaz,epsname{1},'BoundingBox','loose')
        else
            print(fgaz,'-depsc2',epsname{1},'-loose')
        end
    end
end


fgazhist=figure;
axazhist=axes;
hist(psi_tes(:)*180/pi,100);
set(findobj(axazhist,'Type','patch'),...
    'FaceColor',[0.2 0.2 0.2],'EdgeColor','black')
set(axaz,'fontname',fontname,'fontsize',axfontsize)
xlabel(axazhist,'position angle, in degree','fontsize',axfontsize,...
    'fontname',fontname)
ylabel(axazhist,'frequency','fontsize',axfontsize,...
    'fontname',fontname)
titlestr=['Position angle distribution, R.A.: ',num2str(phi_sel*180/pi,3),...
    '^\circ, Dec.: ',num2str(theta_sel*180/pi,3),'^\circ'];
title(axazhist,titlestr,'fontsize',tlfontsize,'fontname',fontname)
figure(fgazhist)
axes(axazhist)
drawnow
if interactive
    epsname=inputdlg('*.eps filename:','Save to EPS file');
else
    epsname={['posangle_dist_',filesuffix,'.eps']};
end
if ~isempty(epsname)
    if ~isempty(epsname{1})
        if exist('epswrite.m','file')==2
            epswrite(fgazhist,epsname{1},'BoundingBox','loose')
        else
            print(fgazhist,'-depsc2',epsname{1},'-loose')
        end
    end
end


%% Simulate observation
[PSF,PHI,THETA]=genPSFClmt(0);
if ~isempty(regexpi(config.SCAN_ACCELERATION,'FFT'))
    obs=getScanDataFFT(PSF,cts_tes*tau,psi_tes);
    if config.POISNOISE_ON
        obs=poissrnd(obs);
    end
    if config.GAUSSNOISE_SIGMA>0
        obs=abs(normrnd(obs,config.GAUSSNOISE_SIGMA));
    end
else
    obs=getScanData(PSF,PHI,THETA,cts_tes*tau,phi_tes,theta_tes,psi_tes);
end

% SNR=10*log10(var(obs(:))/mean(obs(:)));

fgobs=figure;
axobs=axes;
imagesc(u,u,log10(obs)),colormap('hot')
cbax=colorbar;
xlabel(axobs,['x rad, centered at ',num2str(phi_sel*180/pi,3),'^\circ R.A.'],...
    'fontsize',axfontsize,'fontname',fontname)
ylabel(axobs,['y rad, centered at ',num2str(theta_sel*180/pi,3),'^\circ Dec.'],...
    'fontsize',axfontsize,'fontname',fontname)
titlestr='Observed data in selected region';

title(axobs,titlestr,'fontsize',tlfontsize,'fontname',fontname)
set(axobs,'fontname',fontname,'fontsize',axfontsize,'Box','on')
set(cbax,'fontname',fontname,'fontsize',axfontsize,'Box','on')
set(get(cbax,'YLabel'),'string','counts, log to base 10',...
    'fontname',fontname,...
    'fontsize',axfontsize)
figure(fgobs)
axes(axobs)
axis xy
axis tight
axis square
drawnow
if interactive
    epsname=inputdlg('*.eps filename:','Save to EPS file');
else
    epsname={['obs_',filesuffix,'.eps']};
end

if ~isempty(epsname)
    if ~isempty(epsname{1})
        if exist('epswrite.m','file')==2
            epswrite(fgobs,epsname{1},'BoundingBox','loose')
        else
            print(fgobs,'-depsc2',epsname{1},'-loose')
        end
    end
end

simobs=struct('PSF',PSF,...
    'ExpTime',tau,...
    'Obs',obs,...
    'Model',cts_tes,...
    'Bg',rate_bg,...
    'AzObs',psi_tes,...
    'PhiObs',phi_tes,...
    'ModelName',modelname,...
    'ModelFile',modelfile,...
    'RA',phi_sel,...
    'Dec',theta_sel,...
    'ThetaObs',theta_tes,...
    'LocalRad',u);
% save
return
