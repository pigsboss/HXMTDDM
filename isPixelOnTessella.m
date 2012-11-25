function TFList=isPixelOnTessella(q,alpha,PhiList,ThetaList)
%ISPIXELONTESSELLA find if given pixels are on the tessella defined by a
%quaternion as well as an angle alpha. The tessella is about an alpha x alpha
%spherical square. four edges are all arcs of great circles on the sphere.
%the quaternion q defines the orientation of the tessella. PhiList and
%ThetaList are the coordinates of the centers of the pixels.

global config;loadConfig
if isempty(q)
    q=[1 0 0 0];
end
if isempty(alpha)
    alpha=config.WINDOW_SIZE;
end
% find the four vertices of the tessella:
% 1: top left; 2: top right; 3: bottom right; 4: bottom right; 5: center
y=[-1,1,1,-1,0];
z=[1,1,-1,-1,0];
switch config.PIXELIZATION
    case {'radial','r'}
        a=tan(alpha/2);
    case {'parallel','p'}
        a=sin(alpha/2);
    otherwise
        error('Undefined pixelization scheme.')
end
y=y*a;
z=z*a;
x=sqrt(1-y.^2-z.^2);
n=quatARotate(q,[x;y;z],0);
[phi,theta]=xyz2ptr(n(1,:),n(2,:),n(3,:));

% find the outer sphere
R=sqrt(sum((n(:,1)-n(:,5)).^2));
% find the inner sphere
r=R/1.5;

[XList,YList,ZList]=ptr2xyz(PhiList,ThetaList);
TFList=logical(sqrt((XList-n(1,5)).^2+(YList-n(2,5)).^2+(ZList-n(3,5)).^2)<=r);
sz=size(PhiList);
PhiList=PhiList(:);
ThetaList=ThetaList(:);
TFList=TFList(:);
flags=logical(sqrt((XList-n(1,5)).^2+(YList-n(2,5)).^2+(ZList-n(3,5)).^2)<=R);
flags=flags(:);

topisos=find(flags.*(~TFList));
L=length(topisos);
h=waitbar(0,'Please wait...');
prog=0.1;
tic
for l=1:L
    TFList(topisos(l))=pisos(PhiList(topisos(l)),ThetaList(topisos(l)),...
        phi(1),theta(1),phi(2),theta(2),phi(3),theta(3),phi(4),theta(4),...
        phi(5),theta(5));
    if (l/L)>prog
        lasttime=toc;
        etastr=num2str((1-l/L)/(l/L)*lasttime,2);
        waitbar(l/L,h,['Please wait...',etastr,'s left.'])
        prog=prog+0.1;
    end
end
close(h)

TFList=logical(reshape(TFList,sz));
return