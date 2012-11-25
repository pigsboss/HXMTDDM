function [xCube,yCube,zCube,phiCube,thetaCube,rhoCube]=genQuadCube(q,Theta,N)
%GENQUADCUBE Generate quadrilateralized spherical cube. The cube has 4
%essential properties. The first is its location, represented by R.A. and
%Dec. of its center. The second is its orientation, represented by R.A. and
%Dec. of the direction of its axis of symmetry. The third is its dimension,
%represented by the length of one of its sides, in radian. The last is its
%sampling frequency on the spherical surface, represented by number of
%pixels along one of its sides. 
%An ideal quadrilateralized spherical cube is a quadrilateral area on the 
%surface of a unit sphere. Any two adjacent lateral sides are perpendicular
% to each other. All four lateral sides have the same length. All sides are
% arcs on great circles of the sphere. The area is pixelized into N x N 
%pixels. Here an approximate cube is generated.
%Define a reference cube is such a cube that its center is on (1,0,0) and
%its axis of symmetry points to (0,0,1).
%  q is the quaternion represents both the position and the orientation of
%  the cube. In fact the quaternion q represents the rotation from
%  reference cube to the current cube.
%  Theta is the length of one of its lateral sides.
%  N is the number of pixels along one of its lateral sides.


theta=1:N;
theta=((theta/N)-0.5)*Theta;
phi=theta;

xAxis=cos(theta);
yAxis=sin(theta);
zAxis=zeros(size(theta));

axis=[xAxis(:),yAxis(:),zAxis(:)]';
% axisH=[xAxis(:),yAxis(:),zAxis(:)]';
% axisV=quatARotate([cos(pi/4),sin(pi/4),0,0],axisH,0); % 90 around x-axis.
% sideEast=quatARotate([cos(Theta/4),0,0,sin(-1*Theta/4)],axisV,0);
% sideWest=quatARotate([cos(Theta/4),0,0,sin(Theta/4)],axisV,0);
% sideNorth=quatARotate([cos(Theta/4),0,sin(-1*Theta/4),0],axisH,0);
% sideSouth=quatARotate([cos(Theta/4),0,sin(Theta/4),0],axisH,0);
% xCube=[sideEast(1,:),sideWest(1,:),sideNorth(1,:),sideSouth(1,:)];
% yCube=[sideEast(2,:),sideWest(2,:),sideNorth(2,:),sideSouth(2,:)];
% zCube=[sideEast(3,:),sideWest(3,:),sideNorth(3,:),sideSouth(3,:)];

xCube=zeros(N);
yCube=zeros(N);
zCube=zeros(N);
for k=1:N
    quat=[cos(-0.5*phi(k)),0,sin(-0.5*phi(k)),0];
    axis_k=quatARotate(quat,axis,0);
    xCube(k,:)=axis_k(1,:);
    yCube(k,:)=axis_k(2,:);
    zCube(k,:)=axis_k(3,:); 
end

% rotate to position and orientation represented by q:
cube=[xCube(:),yCube(:),zCube(:)]';
cube=quatARotate(q,cube,0);
xCube=reshape(cube(1,:),[N N]);
yCube=reshape(cube(2,:),[N N]);
zCube=reshape(cube(3,:),[N N]);
phiCube=angle(xCube+yCube*1i);
thetaCube=asin(zCube);
rhoCube=sqrt(xCube.^2+yCube.^2+zCube.^2);
return
