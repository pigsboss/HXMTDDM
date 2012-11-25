function Im=addmaz(PSF,Obs,AzObs,ImInit,ImLL,ImUL,NUMIT)

global config;loadConfig
global azGroupInfo
GroupInfoLoaded=false;

if isempty(AzObs),   AzObs=0;                           end
if isempty(ImInit),  ImInit=Obs;                        end
if isscalar(AzObs),  AzObs=ones(size(PSF))*AzObs;       end
if isscalar(ImUL),   ImUL=ones(size(Obs))*ImUL;         end
if isscalar(ImLL),   ImLL=ones(size(Obs))*ImLL;         end
if isscalar(ImInit), ImInit=ones(size(Obs))*ImInit;     end

if ~isempty(azGroupInfo)
    if numel(azGroupInfo{5})==numel(AzObs) &&...
        azGroupInfo{6}==config.AZIMUTH_CLUSTER
        if sum(abs(azGroupInfo{5}(:)-sort(AzObs(:)))) < ...
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
  [iPt,azimuth,groupSize,numGroups,azGroupInfo]=groupScalar(AzObs,...
      config.AZIMUTH_CLUSTER); % group azimuth data
end

OTF=zeros(size(PSF,1),size(PSF,2),numGroups);
IOTF=zeros(size(OTF));
groupMask=false(size(OTF));
for g=1:numGroups
    OTF(:,:,g)=fft2(ifftshift(imrotate(PSF,-1*azimuth(g)*180/pi,...
        'bilinear','crop')));
    IOTF(:,:,g)=fft2(ifftshift(imrotate(rot90(PSF,2),-1*azimuth(g)*180/pi,...
        'bilinear','crop')));
    mask=false(size(PSF));
    mask(iPt{g}(1:groupSize(g)))=true(groupSize(g),1);
    groupMask(:,:,g)=mask;
end

Im=ImInit;
if config.CUDA_ON
    Im=gpuArray(single(Im));
    OTF=gpuArray(single(OTF));
    IOTF=gpuArray(single(IOTF));
    ImUL=gpuArray(single(ImUL));
    ImLL=gpuArray(single(ImLL));
    Obs=gpuArray(single(Obs));
    for t=1:NUMIT
        u=gpuArray(single(zeros(size(PSF))));
        M=u;
        F=fft2(Im);
        for g=1:numGroups
            M=M+single(groupMask(:,:,g)).*real(ifft2(F.*OTF(:,:,g)));
        end
        F=fft2(Obs./M);
        for g=1:numGroups
            u=u+single(groupMask(:,:,g)).*real(ifft2(F.*IOTF(:,:,g)));
        end
        Im=u.*Im;
        % lower limit constraint
        if ~isempty(ImLL),  Im=(Im>ImLL).*(Im-ImLL)+ImLL;   end

        % upper limit constraint
        if ~isempty(ImUL),  Im=(Im<ImUL).*(Im-ImUL)+ImUL;   end
    end
    Im=gather(Im);
else
    for t=1:NUMIT
        u=zeros(size(PSF));
        M=zeros(size(PSF));
        F=fft2(Im);
        for g=1:numGroups
            M=M+groupMask(:,:,g).*real(ifft2(F.*OTF(:,:,g)));
        end
        F=fft2(Obs./M);
        for g=1:numGroups
            u=u+groupMask(:,:,g).*real(ifft2(F.*IOTF(:,:,g)));
        end
        Im=u.*Im;
        % lower limit constraint
        if ~isempty(ImLL),  Im=(Im>ImLL).*(Im-ImLL)+ImLL;   end

        % upper limit constraint
        if ~isempty(ImUL),  Im=(Im<ImUL).*(Im-ImUL)+ImUL;   end
    end
end
return
