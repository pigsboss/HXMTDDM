function pixSz=getPixelSize(varargin)
%GETPIXELSIZE

global config;loadConfig
N=config.SAMPLING_FREQ;
alpha=config.WINDOW_SIZE;
METHOD=config.PIXELIZATION;

pixInfoLoaded=false;

switch length(varargin)
  case 0
  case 1
    switch class(varargin{1})
      case 'double'
        if varargin{1} > 4 % generally window_size should not be greater than 4 (in radian).
          N=varargin{1};
        else
          alpha=varargin{1};
        end
      case 'char'
        METHOD=varargin{1};
      otherwise
        error(['unsupported argument type: (',class(varargin{1}),')',varargin{1}])
    end
  case 2
    if isscalar(varargin{1}) && isscalar(varargin{2})
      alpha=varargin{1};
      N=varargin{2};
    else
      phi=varargin{1};
      theta=varargin{2};
      pixInfo=anaPixel(phi,theta,'mat=custom');
      pixSz=pixInfo.area;
      return
    end
  case 3
    alpha=varargin{1};
    N=varargin{2};
    METHOD=varargin{3};
  otherwise
    error('syntax error: too many input arguments')
end

pixInfoList=dir(['./pixel_',METHOD(1),'*']);
pixInfoFile='';
if ~isempty(pixInfoList)
  NumRows=length(pixInfoList);
  for k=1:NumRows
    row=pixInfoList(k).name;
    WINDOW_SIZE=eval(row((strfind(row,'pixel_')+7):(strfind(row,'d')-1)))*pi/180;
    SAMPLING_FREQ=eval(row((strfind(row,'d')+1):(strfind(row,'.mat')-1)));
    if abs(WINDOW_SIZE-alpha)<1e-5
      pixInfoFile=pixInfoList(k).name;
      if SAMPLING_FREQ==N
        pixInfo=load(pixInfoFile);
        pixInfoLoaded=true;
        break
      end
    end
  end
end
if ~pixInfoLoaded
  if config.LAZYMODE_ON && ~isempty(pixInfoFile)
    pixInfo=load(pixInfoFile);
    [dy,dx]=size(pixInfo.area);
    [X,Y]=meshgrid(1:dx,1:dy);
    X=(X-1)/(dx-1);
    Y=(Y-1)/(dy-1);
    [XI,YI]=meshgrid(1:N);
    XI=(XI-1)/(N-1);
    YI=(YI-1)/(N-1);
    pixSz=interp2(X,Y,pixInfo.area,XI,YI,'linear');
    pixSz=pixSz/sum(pixSz(:))*sum(pixInfo.area(:));
    return
  else
    [~,~,~,~,~,~,pixInfo]=genQuadProj([],alpha,N,METHOD);
  end
end
pixSz=pixInfo.area;
return
