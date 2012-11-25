function quat=getQuatDetStatus(PtPhi,PtTheta,Azimuth)
%GETQUATDETSTATUS Get quaternion of detector status from pointing and
% azimuth, given by PtPhi, PtTheta, and Roll.
% PtPhi and Roll: -pi to pi
% PtTheta: -pi/2 to pi/2
% Azimuth is optional.
%
%The quaternion indicates the rotation which turns the detector from its null
% position as well as orientation to those in its current status.

global config;loadConfig

quatY=[cos(-PtTheta/2),0,sin(-PtTheta/2),0];
quatZ=[cos(PtPhi/2),0,0,sin(PtPhi/2)];
quat=quatAMultiply(quatZ,quatY,config.CUDA_ON);
if nargin==3
    quatP=[cos(Azimuth/2),sin(Azimuth/2)*cos(PtTheta)*cos(PtPhi),...
        sin(Azimuth/2)*cos(PtTheta)*sin(PtPhi),sin(Azimuth/2)*sin(PtTheta)];
    quat=quatAMultiply(quatP,quat,config.CUDA_ON);
end
return
