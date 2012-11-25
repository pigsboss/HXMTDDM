function demoQuadProjTessella(M,N)
%DEMOQUADPROJTESSELLA Demonstrate quadrilateral projection pixelization based
% tessellation for spherical surface.
%
% N is the base resolution of pixelization of each window.
% M is the number of windows along each quarter meridian.
%
% window size:
a=pi/2/M;

% quaternion of each tessella:
q=genQuadProjTessella(M);
[NumWin,~]=size(q);

% color name wheel:
colorNameWheel='ykmcrgb';

for k=1:NumWin
    [x,y,z]=genQuadProj(q(k,:),a,N);
    scatter3(x(:),y(:),z(:),[colorNameWheel(mod(k,7)+1),'.']);
%     scatter3(x(:),y(:),z(:),'b.');
    hold all
end
return