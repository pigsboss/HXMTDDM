function [data,SrcNam,RA,Dec,PosErr,SrcCps,SrcCpsErr,BgrCpsa,SrcExt,SrcExtL,SrcL]=loadRASSBSC(filename)
%LOADRASSBSC Load ROSAT All Sky Survey Bright Source Catalogue.
%Reference: http://www.xray.mpe.mpg.de/rosat/survey/rass-bsc/main/rass-bsc1-fmt.html

fid=fopen(filename,'r');
cols=textscan(fid,'%s%s%f%f%f%s%s%f%f%f%f%f%f%f%f%f%f%f%f%s%f%s%s%d%s','headerlines',4);
fclose(fid);
NumCol=numel(cols);
disp([num2str(NumCol),' columns found in ',filename])
if NumCol~=25
    error('unsupported format.')
end
NumRow=numel(cols{1});
NumSrc=0;
for k=1:NumRow
    if strcmp(cols{6}{k},'_____')
        NumSrc=NumSrc+1;
    end
end
disp([num2str(NumRow),' sources found in ',filename])
disp([num2str(NumRow - NumSrc),' sources are screened.'])
data=cell(NumSrc,1);
SrcNam=cell(NumSrc,1); % Source name
RA=zeros(NumSrc,1); % R.A. (degrees) [0, 2pi)
Dec=zeros(NumSrc,1); % Dec. (degrees) [-pi/2, pi/2]
PosErr=zeros(NumSrc,1); % Positional error (arcsec)
npedm=cell(NumSrc,1); % screening flags
riv=cell(NumSrc,1); % 'new data' flags
SrcCps=zeros(NumSrc,1); % source countrates (counts/sec)
SrcCpsErr=zeros(NumSrc,1); % source countrate errors
BgrCpsa=zeros(NumSrc,1); % background countrate (counts/sec/arcmin^2)
Exposure=zeros(NumSrc,1); % exposure time (sec)
HR1=zeros(NumSrc,1); % hardness ratio 1 (B - A)/(B + A)
HR1Err=zeros(NumSrc,1);
HR2=zeros(NumSrc,1); % hardness ratio 2 (D - C)/(D + C)
HR2Err=zeros(NumSrc,1);
SrcExt=zeros(NumSrc,1); % source extent (arcsec)
SrcExtL=zeros(NumSrc,1); % likelihood of source extent
SrcL=zeros(NumSrc,1); % likelihood of source detection
ExtR=zeros(NumSrc,1); % extraction radius
k=1;
for r=1:NumRow
    if strcmp(cols{6}{r},'_____')
        SrcNam{k}=[cols{1}{r},' ',cols{2}{r}];
        RA(k)=cols{3}(r)*pi/180; % convert to rad
        if RA(k)>pi
            RA(k)=RA(k)-2*pi;
        end
        Dec(k)=cols{4}(r)*pi/180; % convert to rad
        PosErr(k)=cols{5}(r);
        npedm{k}=cols{6}{r};
        riv{k}=cols{7}{r};
        SrcCps(k)=cols{8}(r);
        SrcCpsErr(k)=cols{9}(r);
        BgrCpsa(k)=cols{10}(r)/pi/pi*(60*180)^2; % convert to counts/sec/sr 
        Exposure(k)=cols{11}(r);
        HR1(k)=cols{12}(r);
        HR1Err(k)=cols{13}(r);
        HR2(k)=cols{14}(r);
        HR2Err(k)=cols{15}(r);
        SrcExt(k)=cols{16}(r)*pi/180/3600; % convert to rad
        SrcExtL(k)=cols{17}(r);
        SrcL(k)=cols{18}(r);
        ExtR(k)=cols{19}(r);
        data{k}=struct('SrcNam',SrcNam{k},'RA',RA(k),'Dec',Dec(k),'PosErr',PosErr(k),...
            'npedm',npedm{k},'riv',riv{k},'SrcCps',SrcCps(k),'SrcCpsErr',SrcCpsErr(k),...
            'BgrCpsa',BgrCpsa(k),'Exposure',Exposure(k),'HR1',HR1(k),'HR1Err',HR1Err(k),...
            'HR2',HR2(k),'HR2Err',HR2Err(k),'SrcExt',SrcExt(k),'SrcExtL',SrcExtL(k),...
            'SrcL',SrcL(k),'ExtR',ExtR(k));
        k=k+1;
    end
end
return
