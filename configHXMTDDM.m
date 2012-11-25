function config=configHXMTDDM
%CONFIGHXMTDDM config file builder

if gpuDeviceCount>0
    cuda_default='Yes';
else
    cuda_default='No';
end
userinput=inputdlg({'Window size, in degree:',...
    'Number of pixels along each dimension:',...
    'Pixelization method:',...
    'Cluster size (in radius of the virtual sphere) for azimuth, in arcmin: ',...
    'Turn on CUDA? No/yes: ',...
    'Acceleration method for computing scan data (FFT/convolution/parallel/none): ',...
    'Interpolation acceleration method (on-site ROTATE/PRETRI/none): ',...
    'Turn on poisson noise simulation? No/yes: ',...
    'Standard deviation of simulated gaussian noise (counts/s): ',...
    'Turn on lazy mode? Yes/no:'},'HXMT DDM configuration',1,{...
    '11.25',...
    '512',...
    'radial',...
    '120',...
    cuda_default,...
    'FFT',...
    'rotate',...
    'Yes',...
    '0',...
    'Yes'...
    });
if ~isempty(userinput)
    if ~isempty(userinput{1})
        WINDOW_SIZE=str2double(userinput{1});
    else
        WINDOW_SIZE=11.25;
    end
    WINDOW_SIZE=WINDOW_SIZE/180*pi;

    if ~isempty(userinput{2})
        SAMPLING_FREQ=round(str2double(userinput{2}));
    else
        SAMPLING_FREQ=1024;
    end
    
    PIXELIZATION=userinput{3};
    if isempty(PIXELIZATION)
      PIXELIZATION='radial';
    else
      if regexpi(PIXELIZATION,'^p')
        PIXELIZATION='parallel';
      else
        PIXELIZATION='radial';
      end
    end
    
    
    if ~isempty(userinput{4})
        AZIMUTH_CLUSTER=str2double(userinput{4});
    else
        AZIMUTH_CLUSTER=6;
    end
    AZIMUTH_CLUSTER=AZIMUTH_CLUSTER/60/180*pi;

    CUDA_ON=userinput{5};
    if isempty(CUDA_ON)
      CUDA_ON=logical(gpuDeviceCount);
    else
      switch CUDA_ON
        case {'Y','y','Yes','yes','YES'}
          CUDA_ON=true;
        otherwise
          CUDA_ON=false;
      end
    end

    SCAN_ACCELERATION=userinput{6};
    if isempty(SCAN_ACCELERATION)
      SCAN_ACCELERATION='FFT';
    end

    INTERP_METHOD=userinput{7};
    if isempty(INTERP_METHOD)
      INTERP_METHOD='rotate';
    end
    if strcmpi(INTERP_METHOD,'NONE')
      INTERP_METHOD='';
    end

    POISNOISE_ON=userinput{8};
    if isempty(POISNOISE_ON)
      POISNOISE_ON=true;
    else
      switch POISNOISE_ON
        case {'Y','y','Yes','yes','YES'}
          POISNOISE_ON=true;
        otherwise
          POISNOISE_ON=false;
      end
    end

    
    if ~isempty(userinput{9})
        GAUSSNOISE_SIGMA=abs(str2double(userinput{9}));
    else
        GAUSSNOISE_SIGMA=0;
    end

    LAZYMODE_ON=userinput{10};
    if isempty(LAZYMODE_ON)
      LAZYMODE_ON=true;
    else
      switch LAZYMODE_ON
        case {'n','N','No','no','NO'}
          LAZYMODE_ON=false;
        otherwise
          LAZYMODE_ON=true;
      end
    end
else
    errordlg('Canceled by user.')
    error('Canceled by user.')
end

config=struct('WINDOW_SIZE',WINDOW_SIZE,...
  'SAMPLING_FREQ',SAMPLING_FREQ,...
  'PIXELIZATION',PIXELIZATION,...
  'AZIMUTH_CLUSTER',AZIMUTH_CLUSTER,...
  'SCAN_ACCELERATION',SCAN_ACCELERATION,...
  'CUDA_ON',CUDA_ON,...
  'INTERP_METHOD',INTERP_METHOD,...
  'POISNOISE_ON',POISNOISE_ON,...
  'GAUSSNOISE_SIGMA',GAUSSNOISE_SIGMA,...
  'LAZYMODE_ON',LAZYMODE_ON);

disp(config)

save('config.mat','-struct','config','-mat','-v7.3')
return
