function showTessellaBorderOnMollweide(q,alpha,LineSpec)
%SHOWTESSELLABORDERONMOLLWEIDE Show border of given tessella on mollweide
%projected map.

global config;loadConfig
N=config.SAMPLING_FREQ;
if isempty(q)
    q=[1 0 0 0];
end
if isempty(alpha)
    alpha=config.WINDOW_SIZE;
end
[~,~,~,phi,theta]=genQuadProj(q,alpha,config.SAMPLING_FREQ);
[x,y]=mollweideproj(phi(1,:),theta(1,:));
plot(x,y,LineSpec);
[x,y]=mollweideproj(phi(N,:),theta(N,:));
plot(x,y,LineSpec);
[x,y]=mollweideproj(phi(:,1),theta(:,1));
plot(x,y,LineSpec);
[x,y]=mollweideproj(phi(:,N),theta(:,N));
plot(x,y,LineSpec);
return