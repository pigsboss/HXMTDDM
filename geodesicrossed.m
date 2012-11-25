function result=geodesicrossed(phi_0,theta_0,phi_1,theta_1,...
    phi_2,theta_2,phi_3,theta_3)
%GEODESICROSSED Find whether two geodesics is crossed.
% The two geodesics are defined by two pairs of points on the unit
% spherical surface.
% phi_0, phi_1, phi_2 and phi_3 are the azimuthal angles while theta_0,
% theta_1, theta_2 and theta_3 the altitudes.

[phi_01,theta_01] = normvector(phi_0,theta_0,phi_1,theta_1);
[phi_23,theta_23] = normvector(phi_2,theta_2,phi_3,theta_3);
[phi_r,theta_r]=normvector(phi_01,theta_01,phi_23,theta_23);
if phi_r>0
    phi_l=phi_r-pi;
else
    phi_l=phi_r+pi;
end
theta_l=-1*theta_r;
if isgeodesic(phi_0,theta_0,phi_r,theta_r,phi_1,theta_1) && ...
        isgeodesic(phi_2,theta_2,phi_r,theta_r,phi_3,theta_3)
    result=true;
elseif isgeodesic(phi_0,theta_0,phi_l,theta_l,phi_1,theta_1) && ...
        isgeodesic(phi_2,theta_2,phi_l,theta_l,phi_3,theta_3)
    result=true;
else
    result=false;
end
return
