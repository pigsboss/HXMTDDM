function result=isgeodesic(phi_0,theta_0,phi_1,theta_1,phi_2,theta_2)
%ISGEODESIC Find if arc from point 0 to point 2 through point 1 is
%geodesic.
lat0=theta_0*180/pi;
lat1=theta_1*180/pi;
lat2=theta_2*180/pi;
lon0=phi_0*180/pi;
lon1=phi_1*180/pi;
lon2=phi_2*180/pi;
d01=distance(lat0,lon0,lat1,lon1);
d12=distance(lat1,lon1,lat2,lon2);
d02=distance(lat0,lon0,lat2,lon2);
if abs(d01+d12-d02)<=1e-3
    result=true;
else
    result=false;
end
return