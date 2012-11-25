function n=quatARotate(q,r,ACCELERATE)
%QUATAROTATE Quaternion accelerated rotate.
% q is 1 x 4, r is 3 x n. n has the same size as r.
% ACCELERATE: 0 -> NO Accelerating.
%             1 -> CUDA Accelerating.
[NumQuat,~]=size(q);
[~,NumVec]=size(r);
if nargin==2
  ACCELERATE=0;
end
if ACCELERATE==1
    q=gpuArray(q);
    r=gpuArray(r);
end
w=q(:,1)';
x=q(:,2)';
y=q(:,3)';
z=q(:,4)';
rx=r(1,:);
ry=r(2,:);
rz=r(3,:);
wt=-x*rx-y*ry-z*rz;
xt=w*rx+y*rz-z*ry;
yt=w*ry+z*rx-x*rz;
zt=w*rz+x*ry-y*rx;
if NumQuat>1 && NumVec==1
    n=[ (-wt.*x+xt.*w-yt.*z+zt.*y);...
        (-wt.*y+yt.*w-zt.*x+xt.*z);...
        (-wt.*z+zt.*w-xt.*y+yt.*x)];
elseif NumVec>=1 && NumQuat==1
    n=[ (-wt*x+xt*w-yt*z+zt*y);...
        (-wt*y+yt*w-zt*x+xt*z);...
        (-wt*z+zt*w-xt*y+yt*x)];
else
    error('unsupported dimensions.')
end
if ACCELERATE==1
    n=gather(n);
end
return
