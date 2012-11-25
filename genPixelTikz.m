function genPixelTikz(distH,distV,azH,azV,tikzname)
%GENPIXELTIKZ generate pixel diagram in tikz code (latex code).
x_init=0;
y_init=0;
draw_color='black';
fill_color='lightgray';

if nargin==4
  tikzname=[];
end

if isempty(tikzname)
  tikzname=[datetimefilename,'.tex'];
end

fid=fopen(tikzname,'w');

[~,Nx]=size(distH);
[Ny,~]=size(distV);

l=mean([distH(:);distV(:)]);
distH=distH/l;
distV=distV/l;

x_step=max(distH(:));
y_step=max(distV(:));

for ky=1:Ny
  for kx=1:Nx
    x_pix=x_init+(kx-1)*x_step;
    y_pix=y_init+(ky-1)*y_step;
    tilt=azH(ky,kx)-pi/2;

    x_r=distH(ky,kx)+x_pix;
    y_r=y_pix;

    x_bot=sin(azV(ky,kx)-tilt)*distV(ky,kx)+x_pix;
    y_bot=cos(azV(ky,kx)-tilt)*distV(ky,kx)+y_pix;

    x_botr=sin(azH(ky+1,kx)-tilt)*distH(ky+1,kx)+x_bot;
    y_botr=cos(azH(ky+1,kx)-tilt)*distH(ky+1,kx)+y_bot;
    x_rbot=sin(azV(ky,kx+1)-tilt)*distV(ky,kx+1)+x_r;
    y_rbot=cos(azV(ky,kx+1)-tilt)*distV(ky,kx+1)+y_r;

    x_botr=0.5*(x_botr+x_rbot);
    y_botr=0.5*(y_botr+y_rbot);

    tikzcmd=['\filldraw [fill=', fill_color,',draw=',draw_color,'] (',num2str(x_pix,'%f'),',',num2str(y_pix,'%f'),') -- (',...
        num2str(x_r,'%f'),',',num2str(y_r,'%f'),') -- (',...
        num2str(x_botr,'%f'),',',num2str(y_botr,'%f'),') -- (',...
        num2str(x_bot,'%f'),',',num2str(y_bot,'%f'),') -- cycle;'];
    disp(tikzcmd);
    fprintf(fid,'%s\r\n',tikzcmd);
  end
end
fclose(fid);

disp('Tikz saved to: ');
disp(tikzname);
return
