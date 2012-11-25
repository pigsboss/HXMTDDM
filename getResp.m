function [Resp,PhiSrcLocal,ThetaSrcLocal]=getResp(PSF,PHI,THETA,...
    quatDetectorStatus,...
    PhiSrc,ThetaSrc,...
    METHOD)
%GETRESP Get response to the sources on (PhiSrc,ThetaSrc) while the
%detector is in status characterized by quaternion quatDetectorStatus
%according to the given PSF P(PHI,THETA).
%
% PHI and THETA are generated through meshgrid.
% PSF has the same dimension as PHI's and THETA's.
% PhiSrc and ThetaSrc have the same dimension.
% Returned Resp should have the same dimension as PhiSrc's and ThetaSrc's.

% syntax:
% [Resp,PhiSrcLocal,ThetaSrcLocal] = getResp(PSF,PHI,THETA,quatDetStat,PhiSrc,ThetaSrc,METHOD)
%   default syntax
%
% [Resp,PhiSrcLocal,ThetaSrcLocal] = getResp(PSF,PHI,THETA,quatDetStat,PhiSrc,ThetaSrc)
%   METHOD is omitted. autodetect method using monotonicity of PHI and
%   THETA
%
% [Resp,PhiSrcLocal,ThetaSrcLocal] = getResp(PSF,[],[],quatDetStat,PhiSrc,ThetaSrc)
%   In this mode PSF is interpolant generated by calling function of
%   getResp through triScatteredInterp. 'triangle' method is enforced.

global config;loadConfig

if nargin==6
  monotest=ismonotonic(PHI)*ismonotonic(THETA);
  if monotest==true
    METHOD='meshgrid';
  else
    METHOD='triangle';
  end
end

if isempty(PHI) || isempty(THETA)
  METHOD='triangle'; % enforcing 'triangle' interpolation
end

% record the input dimension of PhiSrc.
szSrc=size(PhiSrc);
% [NumDet,~]=size(quatDetectorStatus);

% serialize PhiSrc and ThetaSrc.
PhiSrc=PhiSrc(:)';
ThetaSrc=ThetaSrc(:)';

% convert R.A. and Dec. into x, y, and z.
r=[cos(ThetaSrc).*cos(PhiSrc);cos(ThetaSrc).*sin(PhiSrc);sin(ThetaSrc)];

% rotate volume of sources into detector's local coordinate system.
n=quatARotate(quatconj(quatDetectorStatus),r,config.CUDA_ON);

% convert x, y, and z into R.A. and Dec. R.A. is from -pi to pi, Dec is
% from -pi/2 to pi/2.
ThetaSrcLocal=asin(n(3,:));
PhiSrcLocal=angle(n(1,:)+n(2,:)*1i);

% compute response through interpolation:
switch METHOD
  case {'triangle','T','tri','t'}
    if ~isempty(PHI) && ~isempty(THETA)
      PSF=TriScatteredInterp(PHI(:),THETA(:),PSF(:),'linear');
    end
    Resp=PSF(PhiSrcLocal,ThetaSrcLocal);
    tmp=zeros(size(Resp));
    tmp(~isnan(Resp))=Resp(~isnan(Resp));
    Resp=tmp;
  case {'meshgrid','grid','g','G','mesh'}
    Resp=interp2(PHI,THETA,PSF,PhiSrcLocal,ThetaSrcLocal,'linear',0);
  otherwise
    error(['unsupported interpolation method: ',METHOD])
end
Resp=Resp.*double(Resp>0);

% recover the shape of Resp.
if numel(PhiSrc)>1
    Resp=reshape(Resp,szSrc);
    PhiSrcLocal=reshape(PhiSrcLocal,szSrc);
    ThetaSrcLocal=reshape(ThetaSrcLocal,szSrc);
end
return
