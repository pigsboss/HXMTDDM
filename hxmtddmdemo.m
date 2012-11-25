function [I,t,obs,bg,psi_obs,simobs,filename]=hxmtddmdemo(varargin)
%HXMTDDMDEMO

% config
rng('shuffle')
global config
global axfontsize
global fontname
global tlfontsize
global gpsf
interactive=false;
datetimestr=datetimefilename;
NUMIT=100;
contour_zoom=1.5;

% generate simulated observed data
if isempty(varargin)
    interactive=true;
    config=configHXMTDDM;
    button=questdlg('Please select a model type for the HXMT DDM demo:',...
        'Select a model type',...
        'Standard point source',...
        'Others',...
        'Others');
    switch button
        case 'Others'
            demo_mode='spec';
        case 'Standard point source'
            demo_mode='std';
        otherwise
            errordlg('Unknown user input.');
    end
else
    loadConfig;
    k=1;
    while k<length(varargin)
        switch varargin{k}
            case 'mode'
                demo_mode=varargin{k+1};
                k=k+2;
            case 'numit'
                NUMIT=varargin{k+1};
                k=k+2;
            case 'window_size'
                config.WINDOW_SIZE=varargin{k+1}*pi/180;
                k=k+2;
            case 'sampling_freq'
                config.SAMPLING_FREQ=varargin{k+1};
                k=k+2;
            case 'cuda_on'
                config.CUDA_ON=varargin{k+1};
                k=k+2;
            case 'contour_zoom'
                contour_zoom=varargin{k+1};
                k=k+2;
            otherwise
                break
        end
    end
end

switch demo_mode
    case 'spec'
        if interactive
            [obs,bg,~,~,psi_obs,simobs]=genSimObsData;
        else
            [obs,bg,~,~,psi_obs,simobs]=...
                genSimObsData(varargin{k:length(varargin)});
        end
        modelstr=strrep(simobs.ModelName,' ','_');
    case 'std'
        if interactive
            [obs,bg,~,~,psi_obs,simobs]=genStdObsData;
        else
            [obs,bg,~,~,psi_obs,simobs]=...
                genStdObsData(varargin{k:length(varargin)});
        end
        modelstr=['crab',int2str(simobs.Sig),'sig'];
    otherwise
        errordlg('Unknown demo mode.');
end
regionstr=[num2str(simobs.RA*180/pi),'d_',num2str(simobs.Dec*180/pi),'d_'];
obs=obs/sum(simobs.PSF(:))/simobs.ExpTime;
I=obs;
PSF=simobs.PSF/sum(simobs.PSF(:));
t=0;
fgim=figure;
axim=axes;
fgct=figure;
axct=axes;
promptstr='Number of iterations to run:';
while NUMIT
    if interactive
        userinput=inputdlg({promptstr},...
            'Check point',1,{'1000'});
        if ~isempty(userinput)
            NUMIT=str2double(userinput{1});
        else
            NUMIT=0;
        end
    end
    if isempty(NUMIT)
        NUMIT=0;
    end
    if NUMIT>0
        tic
        I=addmaz(PSF,obs,psi_obs,I,bg,[],ceil(NUMIT));
        timecost=toc;
        promptstr=['Number of iterations to run, ',...
            num2str(timecost),' seconds for ',...
            int2str(NUMIT),' iterations.'];
        t=t+ceil(NUMIT);
        u=simobs.LocalRad;
        figure(fgim)
        imagesc(u,u,log10(imconv(gpsf,I))),colormap('hot')
        cbax=colorbar;
        xlabel(axim,['x rad, centered at ',...
            num2str(simobs.RA*180/pi,3),'^\circ R.A.'],...
            'fontsize',axfontsize,'fontname',fontname)
        ylabel(axim,['y rad, centered at ',...
            num2str(simobs.Dec*180/pi,3),'^\circ Dec.'],...
            'fontsize',axfontsize,'fontname',fontname)
        titlestr=['Reconstructed image, ',int2str(t),' iterations'];
        title(axim,titlestr,'fontsize',tlfontsize,'fontname',fontname)
        set(axim,'fontname',fontname,'fontsize',axfontsize,'Box','on')
        set(cbax,'fontname',fontname,'fontsize',axfontsize,'Box','on')
        set(get(cbax,'YLabel'),'string','pts/s, log to base 10',...
            'fontname',fontname,...
            'fontsize',axfontsize)
        figure(fgim)
        axis xy
        axis tight
        axis square
        drawnow
        if interactive
            epsname=inputdlg('*.eps filename:','Save to EPS file');
        else
            if strcmp(demo_mode,'std')
                epsname={['crab',num2str(simobs.Sig),...
                    'sig_ddm_image_',num2str(t),'itr_',datetimestr,'.eps']};
            else
                epsname={['ddm_image_',num2str(t),'itr_',datetimestr,'.eps']};
            end
        end
        if ~isempty(epsname)
            if ~isempty(epsname{1})
                if exist('epswrite.m','file')==2
                    epswrite(fgim,epsname{1})
                else
                    print(fgim,'-depsc2',epsname{1},'-loose')
                end
            end
        end
        if interactive
            zoomin=inputdlg('Zoom-in factor for contour plot:',...
                'Draw contour plot');
        else
            zoomin={num2str(contour_zoom)};
        end
        if ~isempty(zoomin)
            if ~isempty(zoomin{1})
                figure(fgct)
                zoomin=max(1,str2double(zoomin{1}));
                u_length=length(u)/zoomin;
                u_start=max(1,round(1+0.5*(length(u)-u_length)));
                u_end=min(length(u),round(u_start+u_length));
                Ig=imconv(gpsf,I);
                contour(axct,u(u_start:u_end)*180*60/pi,...
                    u(u_start:u_end)*180*60/pi,...
                    log10(Ig(u_start:u_end,u_start:u_end)));
                colormap('Jet')
                cbax=colorbar;
                xlabel(axct,['x arcmin, centered at ',...
                    num2str(simobs.RA*180/pi,3),'^\circ R.A.'],...
                    'fontsize',axfontsize,'fontname',fontname)
                ylabel(axct,['y arcmin, centered at ',...
                    num2str(simobs.Dec*180/pi,3),'^\circ Dec.'],...
                    'fontsize',axfontsize,'fontname',fontname)
                if strcmp(demo_mode,'std')
                    titlestr=['Reconstructed image, ',...
                        num2str(simobs.Sig),'\sigma, ',...
                        num2str(simobs.SNR,3),' dB'];
                else
                    titlestr=['Reconstructed image, ',...
                        int2str(t),' iterations'];
                end
                title(axct,titlestr,'fontsize',tlfontsize,...
                    'fontname',fontname)
                set(axct,'fontname',fontname,...
                    'fontsize',axfontsize,'Box','on')
                set(cbax,'fontname',fontname,...
                    'fontsize',axfontsize,'Box','on')
                set(get(cbax,'YLabel'),'string','pts/s, log to base 10',...
                    'fontname',fontname,...
                    'fontsize',axfontsize)
                grid on
                axis xy
                drawnow
                if interactive
                    epsname=inputdlg('*.eps filename:','Save to EPS file');
                else
                    if strcmp(demo_mode,'std')
                        epsname={['crab',num2str(simobs.Sig),...
                            'sig_ddm_contour_',num2str(t),'itr_',...
                            datetimestr,'.eps']};
                    else
                        epsname={['ddm_contour_',num2str(t),'itr_',...
                            datetimestr,'.eps']};
                    end
                end
                if ~isempty(epsname)
                    if ~isempty(epsname{1})
                        if exist('epswrite.m','file')==2
                            epswrite(fgct,epsname{1})
                        else
                            print(fgct,'-depsc2',epsname{1},'-loose')
                        end
                    end
                end
            end
        end
    end
    filename=[modelstr,'_',regionstr,'_',datetimestr,'.mat'];
    save(filename)
    if ~interactive
        break
    end
end
return
