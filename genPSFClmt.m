function [P,X,Y,Z,phi,theta]=genPSFClmt(ClmtType,varargin)
%GENPSFCLMT Generate PSF for collimator of type ClmtType.
% ClmtType 0: rectangular collimator on HE.
% ClmtType 1: square collimator on HE.
%
% syntax:
%
% genPSFClmt(ClmtType)
% 
% genPSFClmt(ClmtType, N)
%
% genPSFClmt(ClmtType, WINDOW_SIZE)
%
% genPSFClmt(ClmtType, PIXELIZATION)
%
% genPSFClmt(ClmtType, N, WINDOW_SIZE)
%
% genPSFClmt(ClmtType, N, WINDOW_SIZE, PIXELIZATION)
%
% genPSFClmt(ClmtType, PHI, THETA)


global config;loadConfig

% FOV size (in degrees);
switch ClmtType
    case 0 % rectangular collimator on HE:
        FOV_THETA=5.7*pi/180;
        FOV_PHI=1.1*pi/180;
    case 1 % square collimator on HE:
        FOV_THETA=5.7*pi/180;
        FOV_PHI=5.7*pi/180;
    otherwise
        error('unknown collimator type.')
end

switch length(varargin)
  case 0
    WINDOW_SIZE=config.WINDOW_SIZE;
    N=config.SAMPLING_FREQ;
    Pixel=config.PIXELIZATION;
    [~,~,~,phi,theta,~]=genQuadProj([],WINDOW_SIZE,N,Pixel);
  case 1
    N=config.SAMPLING_FREQ;
    WINDOW_SIZE=config.WINDOW_SIZE;
    Pixel=config.PIXELIZATION;
    switch class(varargin{1})
      case 'double'
        if varargin{1} > 4 % generally window_size should not be greater than 4 (in radian).
          N=varargin{1};
        else
          WINDOW_SIZE=varargin{1};
        end
      case 'char'
        Pixel=varargin{1};
      otherwise
        error(['unsupported argument type: (',class(varargin{1}),')',varargin{1}])
    end
    [~,~,~,phi,theta,~]=genQuadProj([],WINDOW_SIZE,N,Pixel);
  case 2
    if isscalar(varargin{1}) && isscalar(varargin{2})
      WINDOW_SIZE=varargin{1};
      N=varargin{2};
      Pixel=config.PIXELIZATION;
      [~,~,~,phi,theta,~]=genQuadProj([],WINDOW_SIZE,N,Pixel);
    else
      phi=varargin{1};
      theta=varargin{2};
    end
  case 3
    WINDOW_SIZE=varargin{1};
    N=varargin{2};
    Pixel=varargin{3};
    [~,~,~,phi,theta,~]=genQuadProj([],WINDOW_SIZE,N,Pixel);
  otherwise
    error('syntax error: too many input arguments')
end

P=(double((0.5*FOV_PHI)>abs(phi)).*(0.5*FOV_PHI-abs(phi))).*...
    (double((0.5*FOV_THETA)>abs(theta)).*(0.5*FOV_THETA-abs(theta)));

% calculate coordinates in x,y,z
X=cos(theta).*cos(phi);
Y=cos(theta).*sin(phi);
Z=sin(theta);


% normalize P.
% DeltaTheta=WINDOW_THETA/N*pi/180;
% DeltaPhi=WINDOW_PHI/N*pi/180;
% P=P/sum(P(:).*cos(theta(:)))/DeltaTheta/DeltaPhi;
% [d1,d2]=size(P);
% dw=distance(phi(round(d1/2),round(d2/2))*180/pi,theta(round(d1/2),round(d2/2))*180/pi,...
%   phi(round(d1/2),1+round(d2/2))*180/pi,theta(round(d1/2),1+round(d2/2))*180/pi)/180*pi;
% dh=distance(phi(round(d1/2),round(d2/2))*180/pi,theta(round(d1/2),round(d2/2))*180/pi,...
%   phi(1+round(d1/2),round(d2/2))*180/pi,theta(1+round(d1/2),round(d2/2))*180/pi)/180*pi;
% da=dw*dh;
pixSz=getPixelSize(WINDOW_SIZE,N,Pixel);
P=P.*pixSz/sum(P(:).*pixSz(:));
P=P/max(P(:));
return
