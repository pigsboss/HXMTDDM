function [phi,theta,rho]=xyz2ptr(x,y,z)
%XYZ2PTR convert (x, y, z) coordinates into (phi, theta, rho) coordinates.

rho=sqrt(x.^2+y.^2+z.^2);
z=z./(rho.*double(rho>eps)+double(rho<=eps));
theta=asin(z);
phi=angle(x+y*1i);
return
