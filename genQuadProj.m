function [x,y,z,phi,theta,rho,pixInfo]=genQuadProj(varargin)
%GENQUADCUBE Generate quadrilateral projection on unit spherical surface.
%
% The cube has 4 essential properties. The first is its location,
% represented by R.A. and Dec. of its center.
% The second is its orientation, represented by R.A. and
% Dec. of the direction of its axis of symmetry.
% The third is its span, represented by the radian of projection of its
% axis of symmetry.
% The last is its sampling frequency on the spherical surface, represented
% by number of pixels along one of its sides.
%
%
% Define a reference cube is such a cube that its center is on (1,0,0) and
% its axis of symmetry points to (0,0,1).
%  q is the quaternion represents both the position and the orientation of
%  the cube. In fact the quaternion q represents the rotation from
%  reference cube to the current cube.
%
%  alpha is the radian of the projection of its axis.
%
%  N is the number of pixels along one of its lateral sides.
%
%  METHOD
%    'parallel' projecting pixels parallelly from a quadrilateral
%    'radial' projecting pixels from an tangent quadrilateral
%               towards center of sphere.
%
% Returns:
% (x, y, z) - xyz coordinates of centers of pixels generated.
% (phi, theta, rho) - celestial sphere coordinates of centers of pixels.
% phi the r.a., theta the dec., and rho the radius.
%
% given the area of the original quadrilateral, a square, to be a x a, with
% a being the length of its side, the area of each pixel shall be (a/N)^2.
% this can make calculation of numerical integrations on the spherical surface
% easier.
%
% syntax:
%  genQuadProj(N)
%  genQuadProj(alpha,N)
%  genQuadProj(alpha,METHOD)
%  genQuadProj(q,alpha,N)
%  genQuadProj(q,alpha,N,METHOD)

global config;loadConfig

q=[1 0 0 0];
alpha=config.WINDOW_SIZE;
N=config.SAMPLING_FREQ;
METHOD=config.PIXELIZATION;

switch length(varargin)
  case 0
  case 1
    N=varargin{1};
  case 2
    alpha=varargin{1};
    if ischar(varargin{2})
      METHOD=varargin{2};
    else
      N=varargin{2};
    end
  case 3
    q=varargin{1};
    alpha=varargin{2};
    N=varargin{3};
  case 4
    q=varargin{1};
    alpha=varargin{2};
    N=varargin{3};
    METHOD=varargin{4};
  otherwise
end

if isempty(METHOD)
  METHOD=config.PIXELIZATION;
end

if isempty(q)
  q=[1 0 0 0];
end

if isempty(alpha)
  alpha=config.WINDOW_SIZE;
end

if isempty(N)
  N=config.SAMPLING_FREQ;
end

[y,z]=meshgrid(1:N);

switch METHOD
  case {'parallel','p','Parallel','P'}
    a=2*sin(alpha/2);
    y=(((y-0.5)/(N))-0.5)*a;
    z=(((z-0.5)/(N))-0.5)*a;
    x=sqrt(1-y.^2-z.^2);
    sz=size(x);
    n=[x(:)';y(:)';z(:)'];
    n=quatARotate(q,n);
    x=reshape(n(1,:),sz);
    y=reshape(n(2,:),sz);
    z=reshape(n(3,:),sz);
    [phi,theta,rho]=xyz2ptr(x,y,z);
    matname='p';
  case {'radial','r','Radial','R'}
    a=2*tan(alpha/2);
    y=(((y-0.5)/N)-0.5)*a;
    z=(((z-0.5)/N)-0.5)*a;
    x=ones(size(y));
    sz=size(x);
    n=[x(:)';y(:)';z(:)'];
    n=quatARotate(q,n);
    x=reshape(n(1,:),sz);
    y=reshape(n(2,:),sz);
    z=reshape(n(3,:),sz);
    [phi,theta,~]=xyz2ptr(x,y,z);
    rho=ones(size(phi));
    [x,y,z]=ptr2xyz(phi,theta,rho);
    matname='r';
  otherwise
    error('unsupported method.')
end
matname=[matname,num2str(alpha*180/pi),'d',num2str(N)];

if exist(['./pixel_',matname,'.mat'],'file')~=2
  pixInfo=anaPixel(phi,theta,['mat=',matname]);
else
  pixInfo=load(['./pixel_',matname,'.mat']);
end
return
