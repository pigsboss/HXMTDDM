function [q,phi,theta]=genQuadProjTessella(M)
%GENQUADPROJTESSELLA Generate quadrilateral projection pixelization based
% tessellation for spherical surface.
%
% input argument:
%   M - prime meridian is covered by 4*M tessella.
%
% returns:
%   q - quaternions of all tessella
%   phi - R.A. of centers of all tessella, in rad.
%   theta - Dec. of centers of all tessella, in rad.
%
q=zeros(16*M*M,4);
phi=zeros(16*M*M,1);
theta=zeros(16*M*M,1);
NumWin=0;

% Northern polar window
NumWin=NumWin+1;
phi(NumWin)=0;
theta(NumWin)=pi/2;
q(NumWin,:)=getQuatDetStatus(phi(NumWin),theta(NumWin),0);

% Southern polar window
NumWin=NumWin+1;
phi(NumWin)=0;
theta(NumWin)=-1*pi/2;
q(NumWin,:)=getQuatDetStatus(phi(NumWin),theta(NumWin),0);

% Equatorial windows
for m=1:4*M
    NumWin=NumWin+1;
    phi(NumWin)=(m-1)*pi/(2*M);
    theta(NumWin)=0;
    q(NumWin,:)=getQuatDetStatus(phi(NumWin),theta(NumWin),0);
end
% Others
for m=1:(M-1)
    phiinc=2*atan(tan(pi/4/M)/(cos(m*pi/2/M)+sin(m*pi/2/M)*tan(pi/4/M)));
    N=ceil(2*pi/phiinc);
%   display(N)
    for n=1:N
        NumWin=NumWin+1;
        theta(NumWin)=m*pi/(2*M);
        phi(NumWin)=(n-1)*phiinc;
        q(NumWin,:)=getQuatDetStatus(phi(NumWin),theta(NumWin),0);
        NumWin=NumWin+1;
        theta(NumWin)=-1*m*pi/(2*M);
        phi(NumWin)=(n-1)*phiinc;
        q(NumWin,:)=getQuatDetStatus(phi(NumWin),theta(NumWin),0);
    end
end
q=q(1:NumWin,:);
phi=phi(1:NumWin);
theta=theta(1:NumWin);
return
