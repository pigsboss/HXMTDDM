function Im=addm(PSF,PHI,THETA,...
  ScanData,ScanPtPhi,ScanPtTheta,ScanAzimuth,...
  InitIm,ImLL,ImUL,...
  NUMIT)
%ADDM Accelerated direct demodulation method.
%
%PSF is defined on (PHI,THETA).
%
%ScanData is obtained as ScanPtPhi, ScanPtTheta, ScanAzimuth.
% ScanPtPhi and ScanPtTheta are indicating the pointing of the involved
% detector(s) while ScanAzimuth indicates the orientation of the
% detector(s) with respect to its(their) null position(s). (optional)
%
%DDM parameters:
% InitIm -- initial image for DDM iterations. (optional)
% ImLL -- lower limites of images, pixelwise. (optional)
% ImUL -- upper limites of images, pixelwise. (optional)
%
%NUMIT is the number of iterations to be performed. (optional)
%
%Returns:
%  Im -- the reconstructed image



if isempty(ScanAzimuth)
  ScanAzimuth=0;
end
if isempty(InitIm)
  InitIm=ScanData;
end
if isempty(ImUL)
  ImUL=sum(ScanData(:))*ones(size(ScanData));
end
if isempty(ImLL)
  ImLL=min(ScanData(:))*ones(size(ScanData));
end
if isempty(NUMIT)
  NUMIT=10;
  disp('No number of iterations specified. Perform 10 iterations.');
end

if isscalar(ScanAzimuth)
  ScanAzimuth=ones(size(ScanData))*ScanAzimuth;
end
if isscalar(ImUL)
  ImUL=ones(size(ScanData))*ImUL;
end
if isscalar(ImLL)
  ImLL=ones(size(ScanData))*ImLL;
end
if isscalar(InitIm)
  InitIm=ones(size(ScanData))*InitIm;
end

% Employ R-L iterations for DDM:
Im=itrRL(PSF,PHI,THETA,ScanData,ScanPtPhi,ScanPtTheta,ScanAzimuth,InitIm,ImLL,ImUL,NUMIT);
return



function u=itrRL(PSF,PHI,THETA,d,ScanPtPhi,ScanPtTheta,ScanAzimuth,u0,ul,uu,NUMIT)
u=u0;
IPSF=rot90(PSF,2);
for t=1:NUMIT
  tic;
  disp(['iteration ',num2str(t),' of ',num2str(NUMIT)]);

  % save the previously calculated image:
  uold=u;

  % the R-L iteration:
  M=getScanData(PSF,PHI,THETA,u,ScanPtPhi,ScanPtTheta,ScanAzimuth);
  u=getScanData(IPSF,PHI,THETA,d./M,ScanPtPhi,ScanPtTheta,ScanAzimuth).*u;

  % lower limit constraint
  u=double(u>ul).*(u-ul)+ul;

  % upper limit constraint
  u=double(u<uu).*(u-uu)+uu;

  % display the RMS:
  disp(['RMS: ',num2str(std(u(:)-uold(:)))])
  toc
end
return
