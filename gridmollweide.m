function gridmollweide(phi_inc_deg,theta_inc_deg,LineSpec)
%GRIDMOLLWEIDE

if isempty(LineSpec)
    LineSpec='--k';
end
phi_grid=((-180:phi_inc_deg:180)*pi/180);
theta_grid=(-90:1:90)'*pi/180;
Phi_grid=ones(size(theta_grid))*phi_grid;
Theta_grid=theta_grid*ones(size(phi_grid));
[x,y]=mollweideproj(Phi_grid,Theta_grid);
plot(x,y,LineSpec)
hold all
phi_grid=((-180:1:180)*pi/180);
for theta_grid=-90:theta_inc_deg:90
    [x,y]=mollweideproj(phi_grid,theta_grid*ones(size(phi_grid))*pi/180);
    plot(x,y,LineSpec)
    hold all
end
drawnow
return