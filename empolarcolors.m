function B=empolarcolors(R,I,L,l,d,N,N1,rm,rM)
theta=linspace(0,2*pi*L/d,N);
pos=[R.*cos(theta);R.*sin(theta);theta.*(d/2/pi)];
dtheta=theta(2)-theta(1);
dl=dtheta.*[-R.*sin(theta);R.*cos(theta);ones(1,N).*(d/2/pi)];
r=linspace(rm,rM,N1);
thetap=linspace(0,2*pi,N1);
[r1,theta1]=meshgrid(r,thetap);
r1=r1(:);
theta1=theta1(:);
pos1=r1.*cos(theta1);
pos2=r1.*sin(theta1);
pos3=ones(N1*N1,1).*l;
poss=[pos1.';pos2.';pos3.'];
B=zeros(3,N1*N1);
for i=1:N
    tpos=[pos(1,i);pos(2,i);pos(3,i)]*ones(1,N1*N1);
    r=poss-tpos;
    tB=zeros(3,N1*N1);
    for j=1:N1*N1
        tB(:,j)=I.*cross(dl(:,i),r(:,j))./(norm(r(:,j))^3);
    end
    B=B+tB;
end

Bmag = sqrt(B(1,:).^2 + B(2,:).^2 + B(3,:).^2);  % 模长，1×N1^2
x_obs = pos1;  % 已经算好的 x 坐标（列向量）
y_obs = pos2;  % y 坐标

scatter(x_obs, y_obs, 15, Bmag, 'filled');  % 15为点大小
colormap(jet);
colorbar;
axis equal;
xlabel('x (m)'); ylabel('y (m)');
title(['|B|/10^{-7}T at z = ', num2str(l), ' m']);


end