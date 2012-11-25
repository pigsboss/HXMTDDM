function datScan=getScanDataFFT(PSF,Obj,ScanAzimuth)

global config;loadConfig
global azGroupInfo

datScan=zeros(size(PSF));
GroupInfoLoaded=false;
if isempty(ScanAzimuth)
  ScanAzimuth=0;
end
if isscalar(ScanAzimuth)
  ScanAzimuth=ones(size(PSF))*ScanAzimuth;
end

if ~isempty(azGroupInfo)
    if numel(azGroupInfo{5})==numel(ScanAzimuth) &&...
        azGroupInfo{6}==config.AZIMUTH_CLUSTER
        if sum(abs(azGroupInfo{5}(:)-sort(ScanAzimuth(:)))) < ...
                (config.AZIMUTH_CLUSTER)*0.5
            iPt=azGroupInfo{1};
            azimuth=azGroupInfo{2};
            groupSize=azGroupInfo{3};
            numGroups=azGroupInfo{4};
            GroupInfoLoaded=true;
        end
    end
end

if ~GroupInfoLoaded
  [iPt,azimuth,groupSize,numGroups,azGroupInfo]=groupScalar(ScanAzimuth,...
      config.AZIMUTH_CLUSTER); % group azimuth data
end

for g=1:numGroups
    resp=imrotate(PSF,-1*azimuth(g)*180/pi,'bilinear','crop');
    groupMask=zeros(size(ScanAzimuth));
    groupMask(iPt{g}(1:groupSize(g)))=ones(groupSize(g),1);
    datScan=datScan+imconv(Obj,resp).*groupMask;
end
return
