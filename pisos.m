function result=pisos(phi_test,theta_test,phi_0,theta_0,phi_1,theta_1,...
    phi_2,theta_2,phi_3,theta_3,phi_in,theta_in)
%PISOS Point-in-square-on-sphere test.
% phi_test and theta_test are azimuthal angle and altitude of a given point
% to test.
% (phi_i, theta_i), i=0,1,2,3, are coordinates of four vertices of the
% spherical square, given counter-clockwisely.
% (phi_in, theta_in) is a point which defines the inside of the square.
% If the given point lies inside the square the program returns true while
% returns false if outside.

if geodesicrossed(phi_test,theta_test,phi_in,theta_in,...
        phi_0,theta_0,phi_1,theta_1)
    result=false;
elseif geodesicrossed(phi_test,theta_test,phi_in,theta_in,...
        phi_1,theta_1,phi_2,theta_2)
    result=false;
elseif geodesicrossed(phi_test,theta_test,phi_in,theta_in,...
        phi_2,theta_2,phi_3,theta_3)
    result=false;
elseif geodesicrossed(phi_test,theta_test,phi_in,theta_in,...
        phi_3,theta_3,phi_0,theta_0)
    result=false;
else
    result=true;
end
return
