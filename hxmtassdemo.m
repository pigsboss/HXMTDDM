function hxmtassdemo(M,NUMIT,matlistfile)
%HXMTASSDEMO HXMT all-sky survey demo

global config;loadConfig
if ischar(M)
    M=round(str2double(M));
end

if ischar(NUMIT)
    NUMIT=floor(str2double(NUMIT));
end

MgnSz=2.9;
CnvSz=90/M;
config.WINDOW_SIZE=CnvSz*pi/180;
[~,phi,theta]=genQuadProjTessella(M);
display(['The celestial sphere is divided into ',...
    int2str(length(phi)),' regions.'])
WinSz=CnvSz+MgnSz*2;
ThetaMax=42-WinSz/2;
phi=phi*180/pi;
theta=theta*180/pi;
NumRgn=sum(double(abs(theta)<=ThetaMax));
display([int2str(NumRgn),' regions will be observed.'])
RaRgn=phi(abs(theta)<=ThetaMax);
DecRgn=theta(abs(theta)<=ThetaMax);
fid=fopen(matlistfile,'a');
t=0;
for r=1:NumRgn
    tStart=tic;
    close all;
    [~,~,~,~,~,~,matfilename]=hxmtddmdemo('mode','spec',...
        'contour_zoom',1.0,...
        'window_size',WinSz,...
        'sampling_freq',512,...
        'numit',NUMIT,...
        'ra',RaRgn(r),...
        'dec',DecRgn(r),...
        'modelname','I7YR');
    t=t+toc(tStart);
    fprintf(fid,'%s\r\n',matfilename);
    display([int2str(r),' of ',int2str(NumRgn),...
        ' regions finished. Time elapsed: ',num2str(t),...
        ' seconds. Time remaining: ',num2str(t/r*(NumRgn-r)),' seconds.'])
end
fclose(fid);
return