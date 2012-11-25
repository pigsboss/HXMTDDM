function [n_sate_in,v_sate_in]=azimuthdemo(q,alpha,N)
%AZIMUTHDEMO This program shows the variance of azimuth of PSF.
[~,~,~,phi,theta]=genQuadProj(q,alpha,N);
tic;[n_sate,~,v_sate]=satedemo(3.9,[],'','');toc
[~,L]=size(v_sate);
[phi_sate,theta_sate]=xyz2ptr(n_sate(1,:),n_sate(2,:),n_sate(3,:));
tic;
phi_0=phi(1,1);
phi_1=phi(1,N);
phi_2=phi(N,N);
phi_3=phi(N,1);
theta_0=theta(1,1);
theta_1=theta(1,N);
theta_2=theta(N,N);
theta_3=theta(N,1);
phi_in=phi(round(N/2),round(N/2));
theta_in=theta(round(N/2),round(N/2));
[x_in,y_in,z_in]=ptr2xyz(phi_in,theta_in,1);
[x_0,y_0,z_0]=ptr2xyz(phi_0,theta_0,1);
[x_2,y_2,z_2]=ptr2xyz(phi_2,theta_2,1);
R=sqrt((x_2-x_0)^2+(y_2-y_0)^2+(z_2-z_0)^2)/2;
flags=double(sqrt((n_sate(1,:)-x_in).^2+(n_sate(2,:)-y_in).^2+(n_sate(3,:)-z_in).^2) <= R*1.1);
for l=1:L
    if flags(l)==1
        flags(l)=pisos(phi_sate(l),theta_sate(l),...
            phi_0,theta_0,phi_1,theta_1,phi_2,theta_2,phi_3,theta_3,...
            phi_in,theta_in);
    end
end
toc
n_sate_in=zeros(3,sum(flags));
v_sate_in=zeros(3,sum(flags));
k=1;
for l=1:L
    if flags(l)~=0
        n_sate_in(:,k)=n_sate(:,l);
        v_sate_in(:,k)=v_sate(:,l);
        k=k+1;
    end
end
figure
[v_sate_in_phi,v_sate_in_theta]=...
    xyz2ptr(v_sate_in(1,:),v_sate_in(2,:),v_sate_in(3,:));
n_ref=quatARotate([q(1),-1*q(2:4)],n_sate_in,0);
[n_sate_in_phi,n_sate_in_theta]=xyz2ptr(n_ref(1,:),n_ref(2,:),n_ref(3,:));
quiver(n_sate_in_phi,n_sate_in_theta,v_sate_in_phi,v_sate_in_theta);
% quiver3(n_sate_in(1,:),n_sate_in(2,:),n_sate_in(3,:),...
%     v_sate_in(1,:),v_sate_in(2,:),v_sate_in(3,:))
return
