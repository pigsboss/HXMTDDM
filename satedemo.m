function [n_sate,n_axis,v_sate,orbitsim]=satedemo(tau,N,aviname,epsname,matname)
%SATEDEMO Satellite demo
% This program simulates rotation and precession of a given satellite.
%
% tau: step size, in second
% N:   number of steps, simulated

alpha=43/180*pi; %inclination of the orbit of the satellite, in rad.
%beta=5.3/180*pi/86400; %precession rate, rad/sec.
beta=-5.45/180*pi/86400; %precession rate, rad/sec.
omega=2*pi/95.5/60; %rotation rate, rad/sec.

if isempty(N)
    N=ceil(2*pi/abs(beta)/tau);
end

disp(['simulation time step: ',num2str(tau),' s'])
disp(['simulation steps: ',num2str(N)])
disp(['inclination of orbit: ',num2str(alpha/pi*180),' degrees'])
disp(['precession rate: ',num2str(beta/pi*180*86400),' degrees per day'])
disp(['rotation rate: ',num2str(omega/pi*180*60),' degrees per minute'])


n_sate=zeros(3,N);
phi_sate=zeros(1,N);
theta_sate=zeros(1,N);
n_axis=zeros(3,N);
v_sate=zeros(3,N-1);
% psi_sate=zeros(1,N-1);

n_sate(:,1)=[1;0;0]; %initial position of the satellite (projected to the celestial sphere).
n_axis(:,1)=[0;-1*sin(alpha);cos(alpha)]; %set inclination of axis to alpha.
fontname='Helvetica';
fontsize=14;
if ~isempty(aviname)
    aviobj=VideoWriter(aviname);
    aviobj.FrameRate=30;
    aviobj.Quality=90;
    open(aviobj);
    avihf=figure('Position',[0,0,800,600]);
    figure(avihf)
    aviha=axes('fontname',fontname,'fontsize',fontsize);
end
[phi_sate(1),theta_sate(1),~]=xyz2ptr(-1*n_sate(1,1),n_sate(2,1),n_sate(3,1));
if ~isempty(epsname)
    moll_head=1;
    epshf=figure;
    figure(epshf)
    epsha=axes('fontname',fontname,'fontsize',fontsize);
end
hwb=waitbar(0,'Simulating orbit of satellite. Please wait...');
tic;
prog=0.1;
for n=2:N
    v_sate(:,n-1)=omega*cross(n_axis(:,n-1),n_sate(:,n-1));
    n_sate(:,n)=n_sate(:,n-1)+v_sate(:,n-1)*tau;
    n_sate(:,n)=n_sate(:,n)/sqrt(sum(n_sate(:,n).^2));
    [phi_sate(n),theta_sate(n),~]=xyz2ptr(-1*n_sate(1,n),n_sate(2,n),n_sate(3,n));
%     [~,psi_sate(n-1)]=distance(theta_sate(n-1)*180/pi,phi_sate(n-1)*180/pi,...
%         theta_sate(n)*180/pi,phi_sate(n)*180/pi);
    n_axis(:,n)=[sin(alpha)*sin(beta*(n-1)*tau);...
        -1*sin(alpha)*cos(beta*(n-1)*tau);...
        cos(alpha)];
%     n_sate(:,n)=quatARotate([cos(omega*tau/2),sin(omega*tau/2)*(n_axis(:,n-1))'],n_sate(:,n-1),0);
%     v_sate(:,n-1)=(n_sate(:,n)-n_sate(:,n-1))/tau;
%     n_axis(:,n)=quatARotate([cos(beta*tau*(n-1)/2),0,0,sin(beta*tau*(n-1)/2)],n_axis(:,1),0);
    if ~isempty(epsname)
        if abs(phi_sate(n)-phi_sate(n-1))>pi
            [x,y]=mollweideproj(phi_sate(moll_head:n-1),theta_sate(moll_head:n-1));
            plot(epsha,x,y,'k');
%             drawnow
            hold all
            moll_head=n;
        end
    end
    if ~isempty(aviname)
        plot3(aviha,n_sate(1,1:n),n_sate(2,1:n),n_sate(3,1:n)),grid on,axis square
        xlabel('x','fontname',fontname,'fontsize',fontsize)
        ylabel('y','fontname',fontname,'fontsize',fontsize)
        zlabel('z','fontname',fontname,'fontsize',fontsize)
        title('Scanning path on the celestial sphere','fontname',fontname,'fontsize',fontsize+4)
        axis([-1 1 -1 1 -1 1])
        view(135,30)
        drawnow
        currFrame=getframe(avihf);
        writeVideo(aviobj,currFrame)
    end
    if (n/N)>prog
        lasttime=toc;
        etastr=num2str((1-n/N)/(n/N)*lasttime,2);
        waitbar(n/N,hwb,['Simulating orbit of satellite. Please wait...',etastr,' seconds left.'])
        prog=prog+0.1;
    end
end
close(hwb)
if ~isempty(aviname)
    close(aviobj)
end
if ~isempty(epsname)
    phi_grid=((-180:30:180)*pi/180);
    theta_grid=(-90:1:90)'*pi/180;
    Phi_grid=ones(size(theta_grid))*phi_grid;
    Theta_grid=theta_grid*ones(size(phi_grid));
    [x,y]=mollweideproj(Phi_grid,Theta_grid);
    plot(epsha,x,y,'--k')
    hold all
    phi_grid=((-180:1:180)*pi/180);
    for theta_grid=-90:30:90
        [x,y]=mollweideproj(phi_grid,theta_grid*ones(size(phi_grid))*pi/180);
        plot(epsha,x,y,'--k')
        hold all
    end
    drawnow
    xlabel(epsha,'R.A. rad','fontname',fontname,'fontsize',fontsize);
    ylabel(epsha,'Dec. rad','fontname',fontname,'fontsize',fontsize);
    title(epsha,'Scanning path of all-sky survey, Phase I','fontname',fontname,'fontsize',fontsize+4);
    figure(epshf)
    axes(epsha)
    axis tight
    axis image
    drawnow
%     epswrite(epsname)
end
az_sate=v_sate./(ones(3,1)*sqrt(sum(v_sate.^2,1)));
orbitsim=struct('pt_sate',n_sate,...
    'v_sate',v_sate,...
    'tau',tau,...
    'n_axis',n_axis,...
    'phi_sate',phi_sate,...
    'theta_sate',theta_sate,...
    'clock',((1:N)-1)*tau,...
    'az_sate',az_sate);

if ~isempty(matname)
    save('orbitsim.mat','-struct','orbitsim','-mat','-v7.3')
end
return
