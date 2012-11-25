function q3=quatAMultiply(q1,q2,ACCELERATE)
%QUATAMULTIPLY Quaternion accelerated multiply.
%q1 is 1x4, q2 is nx4. q3 has the same size as q2. or, q1 is nx4, q2 is 1x4
%and q3 has the same size as q1.
%ACCELERATE: 0 -> NO Accelerating.
%            1 -> GPU Accelerating.
if nargin==2
  ACCELERATE=0;
end
switch ACCELERATE
    case 0
        w1=q1(:,1);
        w2=q2(:,1);
        x1=q1(:,2);
        x2=q2(:,2);
        y1=q1(:,3);
        y2=q2(:,3);
        z1=q1(:,4);
        z2=q2(:,4);
        q3=[(w1*w2-x1*x2-y1*y2-z1*z2),...
            (w1*x2+w2*x1+y1*z2-y2*z1),...
            (w1*y2+w2*y1+z1*x2-z2*x1),...
            (w1*z2+w2*z1+x1*y2-x2*y1)];
    case 1
        w1=gpuArray(q1(:,1));
        w2=gpuArray(q2(:,1));
        x1=gpuArray(q1(:,2));
        x2=gpuArray(q2(:,2));
        y1=gpuArray(q1(:,3));
        y2=gpuArray(q2(:,3));
        z1=gpuArray(q1(:,4));
        z2=gpuArray(q2(:,4));
        q3=gather([(w1*w2-x1*x2-y1*y2-z1*z2),...
            (w1*x2+w2*x1+y1*z2-y2*z1),...
            (w1*y2+w2*y1+z1*x2-z2*x1),...
            (w1*z2+w2*z1+x1*y2-x2*y1)]);
    otherwise
        error('unsupported accelerating code.');
end
return