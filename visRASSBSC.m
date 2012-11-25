function [SrcMap,BgrMap,Phi,Theta]=visRASSBSC(N,FrameBoundary,filename)
%VISRASSBSC Visualize RASS BSC.

if isempty(FrameBoundary)
    FrameBoundary=[-pi,pi,-pi/2,pi/2];
end
if isscalar(N)
    NumPhi=N;
    NumTheta=N;
else
    NumPhi=N(1);
    NumTheta=N(2);
end
FrameWidth=FrameBoundary(2)-FrameBoundary(1);
FrameHeight=FrameBoundary(4)-FrameBoundary(3);
[Phi,Theta]=meshgrid(1:NumPhi,1:NumTheta);
Phi=(Phi/NumPhi)*FrameWidth+FrameBoundary(1);
Theta=((Theta-1)/(NumTheta-1))*FrameHeight+FrameBoundary(3);

switch filename
    case 'flat'
        SrcMap=zeros(size(Phi));
        BgrMap=ones(size(Phi));
        return
    case 'delta'
        SrcMap=zeros(size(Phi));
        SrcMap(round(NumPhi/2),round(NumTheta/2))=1;
        BgrMap=zeros(size(Phi));
        return
    otherwise
        [~,~,RA,Dec,~,SrcCps,~,BgrCpsa]=loadRASSBSC(filename);
end
% for periodic boundary at R.A. = -pi and R.A. = pi:
RAL=RA-2*pi;
RAR=RA+2*pi;
% for boundary at Dec. = -pi/2 and Dec. = pi/2:
DecU=pi-Dec;
DecD=-pi-Dec;
RAUL=RA-pi;
RAUR=RA+pi;
RADL=RAUL;
RADR=RAUR;
% extend map:
RAExt=[RAL;RA;RAR;RAUL;RAUR;RADL;RADR];
DecExt=[Dec;Dec;Dec;DecU;DecU;DecD;DecD];
SrcCpsExt=[SrcCps;SrcCps;SrcCps;SrcCps;SrcCps;SrcCps;SrcCps];
BgrCps=BgrCpsa*4*pi/NumPhi/NumTheta;
BgrCpsExt=[BgrCps;BgrCps;BgrCps;BgrCps;BgrCps;BgrCps;BgrCps];

NumSrc=numel(SrcCpsExt);
SrcMap=zeros(size(Phi));
for k=1:NumSrc
    if RAExt(k)>FrameBoundary(1) && RAExt(k)<=FrameBoundary(2) && ...
            DecExt(k)>=FrameBoundary(3) && DecExt(k)<=FrameBoundary(4)
        col=mod(round((RAExt(k)-FrameBoundary(1))*NumPhi/FrameWidth),NumPhi)+1;
        row=mod(round((DecExt(k)-FrameBoundary(3))*(NumTheta-1)/FrameHeight),NumTheta)+1;
        SrcMap(row,col)=SrcMap(row,col)+SrcCpsExt(k);
    end
end
Bgr=TriScatteredInterp(RAExt,DecExt,BgrCpsExt,'linear');
BgrMap=Bgr(Phi,Theta);
return
