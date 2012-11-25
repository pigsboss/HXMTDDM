function pixel=anaPixel(phi,theta,output)
%ANAPIXEL analyze pixelization on spherical surface
%
% phi and theta are in radians
% -pi < phi <= pi, -pi/2 <= theta <= pi/2
%
% output is a string describing the output of the program.
%   supported substrings:
%     'tikz' - will generate latex code to render a pixel diagram by tikz engine.
%     'hist' - will draw histograms of pixelization properties.
%

if nargin==2
  output='';
end

[NumTheta,NumPhi]=size(phi);
HalfNumTheta=ceil(NumTheta/2.0);
HalfNumPhi=ceil(NumPhi/2.0);
phi=phi/pi*180; % convert phi to degrees
theta=theta/pi*180; % convert theta to degrees
distPixH=zeros((HalfNumTheta),(HalfNumPhi));
distPixV=zeros((HalfNumTheta),(HalfNumPhi));
azPixH=zeros((HalfNumTheta),(HalfNumPhi));
azPixV=zeros((HalfNumTheta),(HalfNumPhi));
thetaRight=theta(:,2:NumPhi);
phiRight=phi(:,2:NumPhi);
thetaDown=theta(2:NumTheta,:);
phiDown=phi(2:NumTheta,:);
parfor p=1:HalfNumPhi
  for t=1:HalfNumTheta
    [distPixH(t,p), azPixH(t,p)]=distance(theta(t,p),phi(t,p),thetaRight(t,p),phiRight(t,p));
    [distPixV(t,p), azPixV(t,p)]=distance(theta(t,p),phi(t,p),thetaDown(t,p),phiDown(t,p));
  end
end
% t=NumTheta;
% for p=1:NumPhi-1
%   [distPixH(t,p),azPixH(t,p)]=distance(theta(t,p),phi(t,p),theta(t,p+1),phi(t,p+1));
% end
% p=NumPhi;
% for t=1:NumTheta-1
%   [distPixV(t,p),azPixV(t,p)]=distance(theta(t,p),phi(t,p),theta(t+1,p),phi(t+1,p));
% end
distPixH=distPixH/180*pi;
distPixV=distPixV/180*pi;
azPixH=azPixH/180*pi;
azPixV=azPixV/180*pi;

distPixH=[distPixH,fliplr(distPixH(:,1:(NumPhi-HalfNumPhi)));...
  flipud(distPixH(1:(NumTheta-HalfNumTheta),:)),rot90(distPixH(1:(NumTheta-HalfNumTheta),1:(NumPhi-HalfNumPhi)),2)];
distPixV=[distPixV,fliplr(distPixV(:,1:(NumPhi-HalfNumPhi)));...
  flipud(distPixV(1:(NumTheta-HalfNumTheta),:)),rot90(distPixV(1:(NumTheta-HalfNumTheta),1:(NumPhi-HalfNumPhi)),2)];
azPixH=[azPixH,fliplr(azPixH(:,1:(NumPhi-HalfNumPhi)));...
  flipud(azPixH(1:(NumTheta-HalfNumTheta),:)),rot90(azPixH(1:(NumTheta-HalfNumTheta),1:(NumPhi-HalfNumPhi)),2)];
azPixV=[azPixV,fliplr(azPixV(:,1:(NumPhi-HalfNumPhi)));...
  flipud(azPixV(1:(NumTheta-HalfNumTheta),:)),rot90(azPixV(1:(NumTheta-HalfNumTheta),1:(NumPhi-HalfNumPhi)),2)];

distPixVRight=distPixV(:,2:NumPhi);
distPixHDown=distPixH(2:NumTheta,:);
azPixVRight=azPixV(:,2:NumPhi);
azPixHDown=azPixH(2:NumTheta,:);
aUL=distPixV(1:HalfNumTheta,1:HalfNumPhi).*...
  distPixH(1:HalfNumTheta,1:HalfNumPhi).*...
  abs(sin(azPixV(1:HalfNumTheta,1:HalfNumPhi)-...
  azPixH(1:HalfNumTheta,1:HalfNumPhi)))*0.5;
aDR=distPixVRight(1:HalfNumTheta,1:HalfNumPhi).*...
  distPixHDown(1:HalfNumTheta,1:HalfNumPhi).*...
  abs(sin(azPixVRight(1:HalfNumTheta,1:HalfNumPhi)-...
  azPixHDown(1:HalfNumTheta,1:HalfNumPhi)))*0.5; % area of lower right triangle.
areaPix=aUL+aDR;
areaPix=[areaPix,fliplr(areaPix(:,1:(NumPhi-HalfNumPhi)));...
  flipud(areaPix(1:(NumTheta-HalfNumTheta),:)),rot90(areaPix(1:(NumTheta-HalfNumTheta),1:(NumPhi-HalfNumPhi)),2)];

if ~isempty(strfind(output,'hist'))
  figure;hist(distPixV(:),30);title('distances between adjacent pixels, vertical');xlabel('distance, in radian');ylabel('frequency');drawnow
  figure;hist(distPixH(:),30);title('distances between adjacent pixels, horizontal');xlabel('distance, in radian');ylabel('frequency');drawnow
  figure;hist(areaPix(:),30);title('pixel area');xlabel('area, in steradian');ylabel('frequency');drawnow
end

if ~isempty(strfind(output,'tikz'))
  genPixelTikz(distPixH,distPixV,azPixH,azPixV)
end

if ~isempty(strfind(output,'mat='))
  namestart=strfind(output,'mat=')+4;
  namestr=['pixel_',output(namestart:length(output)),'.mat'];
else
  namestr='pixel.mat';
end

pixel=struct('area',areaPix,...
  'distRight',distPixH,...
  'distDown',distPixV,...
  'azRight',azPixH,...
  'azDown',azPixV);
save(namestr,'-struct','pixel','-v7.3')
return
