function [obs,rate_bg,phi_tes,theta_tes,psi_tes,simobs]=...
    genStdObsData(varargin)
%CALHXMTDDM HXMT DDM calibration.

% config
global config;loadConfig
global axfontsize
global fontname
global tlfontsize
global gpsf

filesuffix=datetimefilename;
interactive=false;
tau=3.9;
crab_sig=20;
fontname='Helvetica';
axfontsize=20;
tlfontsize=axfontsize+4;
fwaxfontsize=14;
fwtlfontsize=fwaxfontsize+2;

if isempty(varargin)
    interactive=true;
    userinput=inputdlg({...
        'Sampling interval of orbital simulation, in seconds:',...
        'Significance of standard point source:',...
        'EPS fontname:',...
        'EPS axes fontsize:',...
        'EPS title fontsize:',...
        'EPS axes fontsize for full-line-width figures:',...
        'EPS title fontsize for full-line-width figures:',...
        'Number of pixels along each side of a tessella:',...
        'Tessella size, in degree:'...
        },'User input',1,{...
        '3.9',...
        '20',...
        'Helvetica',...
        '20',...
        '24',...
        '14',...
        '16',...
        int2str(config.SAMPLING_FREQ),...
        num2str(config.WINDOW_SIZE*180/pi)...
        });
    if isempty(userinput)
        errordlg('Canceled by user.')
        error('Canceled by user.')
    end
    if ~isempty(userinput{1}),  tau=str2double(userinput{1});end
    if ~isempty(userinput{2}),  crab_sig=str2double(userinput{2});end
    if ~isempty(userinput{3}),  fontname=userinput{3};end
    if ~isempty(userinput{4}),  axfontsize=str2double(userinput{4});end
    if ~isempty(userinput{5}),  tlfontsize=str2double(userinput{5});end
    if ~isempty(userinput{6}),  fwaxfontsize=str2double(userinput{6});end
    if ~isempty(userinput{7}),  fwtlfontsize=str2double(userinput{7});end
    if ~isempty(userinput{8}),  config.SAMPLING_FREQ=str2double(userinput{8});end
    if ~isempty(userinput{9}),  config.WINDOW_SIZE=str2double(userinput{9})*pi/180;end
else
    k=1;
    while k<length(varargin)
        switch varargin{k}
            case 'sig'
                crab_sig=varargin{k+1};
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
            otherwise
                error('Unsupported parameter.')
        end
    end
end

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
    az_sate=orbitsim.pt_sate;
    phi_sate=orbitsim.phi_sate;
    theta_sate=orbitsim.theta_sate;
end
L=size(az_sate,2);
psi_sate=getAzimuth(pt_sate(:,1:L),az_sate);

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
title(axallsky,'Sky region selection',...
    'fontsize',fwtlfontsize,'fontname',fontname)
set(axallsky,'fontname',fontname,'fontsize',fwaxfontsize,'Box','on')

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

quat=getQuatDetStatus(phi_sel,theta_sel);
showTessellaBorderOnMollweide(quat,[],'-r');
axis xy
axis image
axis tight
drawnow
figure(fgallsky)
axes(axallsky)
[~,~,~,phi_tes,theta_tes]=genQuadProj(quat,[],config.SAMPLING_FREQ);

[rate_tes,rate_bg]=genStdModel(crab_sig,tau);
u=(((1:config.SAMPLING_FREQ)-1)/(config.SAMPLING_FREQ-1)-0.5)*config.WINDOW_SIZE;
gpsf=fspecial('gauss',config.SAMPLING_FREQ,1);
fgtes=figure;
axtes=axes;
im_tes=imconv(rate_tes,gpsf);
imagesc(u,u,log10(im_tes)),colormap('hot')
cbax=colorbar;
xlabel(axtes,['x rad, centered at ',num2str(phi_sel*180/pi,3),'^\circ R.A.'],...
    'fontsize',axfontsize,'fontname',fontname)
ylabel(axtes,['y rad, centered at ',num2str(theta_sel*180/pi,3),'^\circ Dec.'],...
    'fontsize',axfontsize,'fontname',fontname)
titlestr=['Standard point source in selected region, ',num2str(crab_sig),'\sigma'];
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
    epsname={['crab',num2str(crab_sig),'sig_model_',filesuffix,'.eps']};
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
[PSF,PHI,THETA]=genPSFClmt(0);
if ~isempty(regexpi(config.SCAN_ACCELERATION,'FFT'))
    obs=getScanDataFFT(PSF,rate_tes*tau,psi_tes);
    if config.POISNOISE_ON
        obs=poissrnd(obs);
    end
    if config.GAUSSNOISE_SIGMA>0
        obs=abs(normrnd(obs,config.GAUSSNOISE_SIGMA));
    end
else
    obs=getScanData(PSF,PHI,THETA,rate_tes*tau,phi_tes,theta_tes,psi_tes);
end

obs_smooth=imconv(fspecial('gauss',config.SAMPLING_FREQ,10),obs);
mu_signal=obs_smooth(round(config.SAMPLING_FREQ/2),round(config.SAMPLING_FREQ/2))...
    - tau*sum(PSF(:))*rate_bg;
snr=10*log10(mu_signal/sqrt(tau*sum(PSF(:))*rate_bg));

fgobs=figure;
axobs=axes;
imagesc(u,u,log10(obs)),colormap('hot')
cbax=colorbar;
xlabel(axobs,['x rad, centered at ',num2str(phi_sel*180/pi,3),'^\circ R.A.'],...
    'fontsize',axfontsize,'fontname',fontname)
ylabel(axobs,['y rad, centered at ',num2str(theta_sel*180/pi,3),'^\circ Dec.'],...
    'fontsize',axfontsize,'fontname',fontname)
titlestr={['Observed data, ',num2str(crab_sig),'\sigma, ',num2str(snr,3),' dB']};

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
    epsname={['crab',num2str(crab_sig),'sig_obs_',filesuffix,'.eps']};
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
    'Model',rate_tes,...
    'Bg',rate_bg,...
    'AzObs',psi_tes,...
    'PhiObs',phi_tes,...
    'ThetaObs',theta_tes,...
    'LocalRad',u,...
    'SNR',snr,...
    'RA',phi_sel,...
    'Dec',theta_sel,...
    'Sig',crab_sig);

return
