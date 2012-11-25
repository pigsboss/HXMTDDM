function hxmtassviewer(filename)
%HXMTASSVIEWER

load(filename)

% draw image
fgim=figure;
axim=axes;
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

% draw contour diagram
fgct=figure;
axct=axes;
figure(fgct)
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
return