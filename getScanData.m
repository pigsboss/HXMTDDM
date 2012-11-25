function datScan=getScanData(PSF,PHI,THETA,...
  Obj,ObjPhi,ObjTheta,...
  varargin)
%GETSCANDATA Calculate scan data of the given object.
%
%PSF: normalized PSF defined on (PHI,THETA) of local coordinate system of
%the detector.
%
%Obj is defined on (ObjPhi,ObjTheta).
%
%ScanPtPhi and ScanPtTheta are two parameters of pointing of the detector.
%ScanAzimuth is the parameter of the orientation of the detector.
%
%SimulatePoisNoise: logical, true or false.
%
%GaussNoiseSigma: standard deviation of simulated gaussian noise.
%
% ACCEL_METHOD: acceleration methods.
%   'none': no acceleration at all (default).
%
%   'parallel': parallel acceleration.
%
%   'shift': using shifting to accelerate calculation of response matrices
%     sharing the same azimuth angle. to employ this acceleration method
%     the scan points, i.e., the grids of detector status, must have the
%     same sampling interval as the grids of object space. Or the sampling
%     interval of the detector status space must be N times of the sampling
%     interval of the object space. N must be integer. If not the response
%     matrices calculated with this acceleration method would be
%     inaccurate.
%
%   'fft': using fft accelerated circular convolution.
%
%   'pretri': pre-triangulation. using pre-built triangulation scheme for
%     all scattered interpolation following. this scheme can be used
%     together with other accelerations.
%
%
%Returns:
%  datScan. defined on (ScanPtPhi, ScanPtTheta, ScanAzimuth).
%

% declare global variables:
global config;loadConfig
global PSFInterp
global azGroupInfo

% default values of optional arguments:
ScanPtPhi=ObjPhi;
ScanPtTheta=ObjTheta;
ScanAzimuth=zeros(size(ScanPtTheta)); % no azimuth component in scan point vectors
SimulatePhotonNoise=config.POISNOISE_ON;
GaussianNoiseSigma=config.GAUSSNOISE_SIGMA;
ACCEL_METHOD=[config.INTERP_METHOD,config.SCAN_ACCELERATION];

% parsing input arguments:
SAMEGRID=true;
switch length(varargin)
  case 0
  case 1
    if ischar(varargin{1})
      ACCEL_METHOD=varargin{1};
    else
      ScanAzimuth=varargin{1};
    end
  case 2
    ScanAzimuth=varargin{1};
    if ischar(varargin{2})
      ACCEL_METHOD=varargin{2};
    else
      SimulatePhotonNoise=varargin{2};
    end
  case 3
    ScanAzimuth=varargin{1};
    SimulatePhotonNoise=varargin{2};
    if ischar(varargin{3})
      ACCEL_METHOD=varargin{3};
    else
      GaussianNoiseSigma=varargin{3};
    end
  case 4
    ScanAzimuth=varargin{1};
    SimulatePhotonNoise=varargin{2};
    GaussianNoiseSigma=varargin{3};
    ACCEL_METHOD=varargin{4};
  case 6
    ScanPtPhi=varargin{1};
    ScanPtTheta=varargin{2};
    ScanAzimuth=varargin{3};
    SimulatePhotonNoise=varargin{4};
    GaussianNoiseSigma=varargin{5};
    ACCEL_METHOD=varargin{6};
    if numel(ScanPtPhi)~=numel(ObjPhi) || numel(ScanPtTheta)~=numel(ObjTheta)
      SAMEGRID=false;
    elseif sum(abs(ScanPtPhi(:)-ObjPhi(:)))+sum(abs(ScanPtTheta(:)-ObjTheta(:)))~=0
      SAMEGRID=false;
    end
  otherwise
    error('syntax error')
end

if isempty(ACCEL_METHOD)
  ACCEL_METHOD='pretrifft';
end
if isempty(ScanAzimuth)
  ScanAzimuth=0;
end
if isempty(SimulatePhotonNoise)
  SimulatePhotonNoise=false;
end
if isempty(GaussianNoiseSigma)
  GaussianNoiseSigma=0;
end
if isempty(ScanPtPhi)
  ScanPtPhi=ObjPhi;
end
if isempty(ScanPtTheta)
  ScanPtTheta=ObjTheta;
end
if isempty(ScanAzimuth)
  ScanAzimuth=0;
end
if isscalar(ScanAzimuth)
  ScanAzimuth=ones(size(ScanPtTheta))*ScanAzimuth;
end

NumScanPts=numel(ScanPtPhi);
datScan=zeros(size(ScanPtPhi));
[NumTheta,NumPhi]=size(Obj);
[NumPSFTheta,NumPSFPhi]=size(PSF);
PSFCtrIndexPhi=ceil((NumPSFPhi+1)/2);
PSFCtrIndexTheta=ceil((NumPSFTheta+1)/2);
PSFCtrPhi=PHI(PSFCtrIndexTheta,PSFCtrIndexPhi);
PSFCtrTheta=THETA(PSFCtrIndexTheta,PSFCtrIndexPhi);
quatPSF=getQuatDetStatus(PSFCtrPhi,PSFCtrTheta,0);
GroupInfoLoaded=false;
ROTATE_ON=false;

if ~isempty(regexpi(ACCEL_METHOD,'rotate'))
  ROTATE_ON=true;
  ACCEL_METHOD=regexprep(ACCEL_METHOD,'rotate','','ignorecase');
end

if ~isempty(regexpi(ACCEL_METHOD,'pretri'))
  if ROTATE_ON
    error('ROTATING acceleration conflicts with PRETRI acceleration.');
  end
  ACCEL_METHOD=regexprep(ACCEL_METHOD,'pretri','','ignorecase');
  PSFVal=PSF;
  if config.LAZYMODE_ON
    if ~isempty(PSFInterp)
      PSF=PSFInterp;
    else
      PSF=TriScatteredInterp(PHI(:),THETA(:),PSF(:),'linear');
      PSFInterp=PSF;
    end
  else
    PSF=TriScatteredInterp(PHI(:),THETA(:),PSF(:),'linear');
    PSFInterp=PSF;
  end
  PHI=[];
  THETA=[];
end






switch ACCEL_METHOD
  case {'parallel','par','Parallel','Par','p','P'}
    parfor k=1:NumScanPts
      quatObj=getQuatDetStatus(ScanPtPhi(k),ScanPtTheta(k),ScanAzimuth(k));
      quat=quatAMultiply(quatconj(quatPSF),quatObj); % calculate the central quaternion of the current group
      R=getResp(PSF,PHI,THETA,quat,ObjPhi,ObjTheta);
      datScan(k)=sum(sum(Obj.*R));
    end

  case {'CONVOLUTION','CONV','conv','convolution','shift','shifting','Shift','Shifting','s','S'}
    if ~SAMEGRID
      error('Shift acceleration requires the scan grid to be the same as the object grid.')
    end
    if config.LAZYMODE_ON
      if ~isempty(azGroupInfo)
        if numel(azGroupInfo{5})==numel(ScanAzimuth) && azGroupInfo{6}==config.AZIMUTH_CLUSTER
          if sum(abs(azGroupInfo{5}(:)-ScanAzimuth(:)))<1e-8
            iPt=azGroupInfo{1};
            azimuth=azGroupInfo{2};
            groupSize=azGroupInfo{3};
            numGroups=azGroupInfo{4};
            GroupInfoLoaded=true;
          end
        end
      end
    end
    if ~GroupInfoLoaded
      [iPt,azimuth,groupSize,numGroups,azGroupInfo]=groupScalar(ScanAzimuth, config.AZIMUTH_CLUSTER); % group azimuth data
    end
    CtrIndexPhi=ceil((NumPhi+1)/2);
    CtrIndexTheta=ceil((NumTheta+1)/2);
    ObjCtrPhi=ObjPhi(CtrIndexTheta, CtrIndexPhi);
    ObjCtrTheta=ObjTheta(CtrIndexTheta, CtrIndexPhi); % calculate the global central RA and Dec (pointing vector).
    for g=1:numGroups
      quatObj=getQuatDetStatus(ObjCtrPhi,ObjCtrTheta,azimuth(g));
      quat=quatAMultiply(quatconj(quatPSF),quatObj); % calculate the central quaternion of the current group
      if sum(abs(quat(2:4)))<1e-12
        if isfloat(PSF)
          resp=PSF;
        else
          resp=PSFVal;
        end
      else
        if ROTATE_ON
          resp=imrotate(PSF,-1*azimuth(g)*180/pi,'bilinear','crop');
        else
          resp=getResp(PSF,PHI,THETA,quat,ObjPhi,ObjTheta); % calculate the central response matrix of the current group
        end
      end
      resp=circshift(rot90(resp,2),[2*CtrIndexTheta-NumTheta-1, 2*CtrIndexPhi-NumPhi-1]);
      for k=1:groupSize(g)
        IndexTheta=mod(iPt{g}(k)-1,NumTheta)+1;
        IndexPhi=((iPt{g}(k)-IndexTheta)/NumTheta)+1;
        PhiShift=IndexPhi-CtrIndexPhi;
        ThetaShift=IndexTheta-CtrIndexTheta;
        R=zeroshift2(resp,[ThetaShift PhiShift]);
        datScan(iPt{g}(k))=sum(Obj(:).*R(:));
      end
    end

  case {'fft','FFT','f','F'}
    if ~SAMEGRID
      error('FFT acceleration requires the scan grid to be the same as the object grid.')
    end
    if config.LAZYMODE_ON
      if ~isempty(azGroupInfo)
        if numel(azGroupInfo{5})==numel(ScanAzimuth) && azGroupInfo{6}==config.AZIMUTH_CLUSTER
          if sum(abs(azGroupInfo{5}(:)-ScanAzimuth(:)))<1e-8
            iPt=azGroupInfo{1};
            azimuth=azGroupInfo{2};
            groupSize=azGroupInfo{3};
            numGroups=azGroupInfo{4};
            GroupInfoLoaded=true;
          end
        end
      end
    end
    if ~GroupInfoLoaded
      [iPt,azimuth,groupSize,numGroups,azGroupInfo]=groupScalar(ScanAzimuth, config.AZIMUTH_CLUSTER); % group azimuth data
    end
    CtrIndexPhi=ceil((NumPhi+1)/2);
    CtrIndexTheta=ceil((NumTheta+1)/2);
    ObjCtrPhi=ObjPhi(CtrIndexTheta, CtrIndexPhi);
    ObjCtrTheta=ObjTheta(CtrIndexTheta, CtrIndexPhi); % calculate the global central RA and Dec (pointing vector).
    parfor g=1:numGroups
      quatObj=getQuatDetStatus(ObjCtrPhi,ObjCtrTheta,azimuth(g));
      quat=quatAMultiply(quatconj(quatPSF),quatObj); % calculate the central quaternion of the current group
      if ROTATE_ON
        resp=imrotate(PSF,-1*azimuth(g)*180/pi,'bilinear','crop');
      else
        resp=getResp(PSF,PHI,THETA,quat,ObjPhi,ObjTheta); % calculate the central response matrix of the current group
      end
      groupMask=zeros(size(ScanAzimuth));
      groupMask(iPt{g}(1:groupSize(g)))=ones(groupSize(g),1);
      datScan=datScan+imconv(Obj,resp).*groupMask;
    end

  otherwise % no acceleration at all
    for k=1:NumScanPts
      quatObj=getQuatDetStatus(ScanPtPhi(k),ScanPtTheta(k),ScanAzimuth(k));
      quat=quatAMultiply(quatconj(quatPSF),quatObj); % calculate the central quaternion of the current group
      R=getResp(PSF,PHI,THETA,quat,ObjPhi,ObjTheta);
      datScan(k)=sum(sum(Obj.*R));
    end
end

if SimulatePhotonNoise==true
  datScan=round(normrnd(datScan,sqrt(datScan)));
end
if GaussianNoiseSigma>1e-10
  datScan=(normrnd(datScan,GaussianNoiseSigma));
end
return
