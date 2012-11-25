function [xTes,yTes,zTes,RATes,DecTes,xPix,yPix,zPix,RAPix,DecPix]=...
    genQuadCubeTessellation(NumPixTes,Margin,THETA)
%GENQUADCUBETESSELLATION Generate quadrilateralized spherical cube
%tessellation of the celestial sphere.
%
%  THETA: square tessellation size in degrees.
%  Margin: margin of each tessellation in degrees. the tessellation 
%  generated is larger than THETA by 2*Margin along each side.
%  NumPixTes: number of pixels along a lateral side of the square 
%    tessellation.

if nargin<2
    THETA=11.25;
    Margin=0.0;
end
if nargin<3
    THETA=11.25;
end
Margin=Margin/THETA;
THETA=THETA/180*pi;

% number of tessellations along each meridian:
NumTesTheta=ceil(pi/THETA)+1;

% Dec. of tessellations:
ThetaTes=(((1:NumTesTheta)-1)/(NumTesTheta-1))*pi-pi/2;

% number of tessellations along each latitude:
NumTesPhi=ceil(2*pi/THETA);

PhiTes=(((1:NumTesPhi)-1)/NumTesPhi)*2*pi-pi;

[RATes,DecTes]=meshgrid(PhiTes,ThetaTes);
xTes=cos(DecTes).*cos(RATes);
yTes=cos(DecTes).*sin(RATes);
zTes=sin(DecTes);

NumPixTesOuter=ceil(NumPixTes*(1+2.0*Margin));
THETAOUTER=THETA*NumPixTesOuter/NumPixTes;
xPix=zeros(NumPixTesOuter,NumPixTesOuter,NumTesTheta,NumTesPhi);
yPix=xPix;
zPix=yPix;
RAPix=xPix;
DecPix=xPix;
for k=1:NumTesPhi
    for l=1:NumTesTheta
        [x,y,z]=genQuadCube(getQuatDetStatus(PhiTes(k),ThetaTes(l),0),...
            THETAOUTER,NumPixTesOuter);
        xPix(:,:,l,k)=x;
        yPix(:,:,l,k)=y;
        zPix(:,:,l,k)=z;
        RAPix(:,:,l,k)=angle(x+y*1i);
        DecPix(:,:,l,k)=asin(z);
    end
end
return