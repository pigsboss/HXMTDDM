function [x,y,z]=ptr2xyz(phi,theta,rho)
%PTR2XYZ convert (phi, theta, rho) coordinates into (x,y,z) coordinates.
if nargin==2
  rho=ones(size(phi));
end
z=sin(theta).*rho;
x=cos(phi).*cos(theta).*rho;
y=sin(phi).*cos(theta).*rho;
return
