function [phi_n,theta_n]=normvector(phi_0,theta_0,phi_1,theta_1)
%NORMVECTOR Find normal vector of plane defined by the origin and another
%two points on the unit spherical surface. Right-hand rule applied.
%
% phi_0 and phi_1 are the azimuthal angles of the given two points. The
% azimuthal angle goes from 0 to 2*pi.
% theta_0 and theta_1 are the altitudes, which go from -pi/2 to pi/2.
%
% phi_n is the azimuthal angle of the normal vector while theta_n is the
% altitude.
%
% n = cross(p0, p1), where n is the normal vector, i.e., the cross product
% of p0 and p1.
%
% This program is faster ( > 3x) than the build-in matlab cross program.

phi_d=phi_1-phi_0;
d0=cos(theta_0).*sin(theta_1);
d1=sin(theta_0).*cos(theta_1);
d=d0.^2 + d1.^2 - 2*d0.*d1.*cos(phi_d);
c=cos(theta_0).*cos(theta_1).*sin(phi_d);
theta_n=angle(sqrt(d)+c*1i);
x=d0.*sin(phi_0) - d1.*sin(phi_1);
y=d1.*cos(phi_1) - d0.*cos(phi_0);
phi_n=angle(x+y*1i);
return
