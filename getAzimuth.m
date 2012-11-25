function psi=getAzimuth(pt,az)
%GETAZIMUTH get azimuth angle from pointing vector and azimuthal vector.
% the azimuth angle is defined as the angle between the azimuthal vector
% and the meridian (from south to north) at the pointing vector.

% find the tangent vector of the local meridian:
[phi_pt,theta_pt]=xyz2ptr(pt(1,:),pt(2,:),pt(3,:));
onnorth=logical(theta_pt>0);
phi_md=phi_pt;
theta_md=pi/2-theta_pt;
phi_md(onnorth)=phi_pt(onnorth)-pi;
theta_md(~onnorth)=theta_pt(~onnorth)+pi/2;
[x_md,y_md,z_md]=ptr2xyz(phi_md,theta_md,ones(size(phi_md)));

% cross product, sin(psi):
y=sqrt(sum(cross(az,[x_md;y_md;z_md]).^2,1))./sqrt(sum(az.^2,1));

% dot product, cos(psi):
x=dot(az,[x_md;y_md;z_md])./sqrt(sum(az.^2,1));

psi=angle(x+1i*y);
return