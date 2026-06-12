function B=emverticalstream2(R,I,L,h1,d,N,N1,s,h2)

theta=linspace(0,2*pi*L/d,N);
pos=[R.*cos(theta);R.*sin(theta);theta.*(d/2/pi)];
dtheta=theta(2)-theta(1);
dl=dtheta.*[-R.*sin(theta);R.*cos(theta);ones(1,N).*(d/2/pi)];
x=linspace(-s,s,N1);
z=linspace(h1,h2,N1);
[pos1,pos3]=meshgrid(x,z);
pos1=pos1(:);
pos2=zeros(N1*N1,1);
pos3=pos3(:);
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

X_grid = reshape(poss(1,:), N1, N1);
Z_grid = reshape(poss(3,:), N1, N1);
Bx_grid = reshape(B(1,:), N1, N1);
By_grid = reshape(B(2,:), N1, N1);
Bz_grid = reshape(B(3,:), N1, N1);

% 计算磁场模长 |B|
Bmag = sqrt(Bx_grid.^2 + By_grid.^2 + Bz_grid.^2);

% 创建新图窗，两个子图并排
figure;

% 子图1：磁力线图 (streamslice)
subplot(1,2,1);
streamslice(X_grid, Z_grid, Bx_grid, Bz_grid, 1.5);
axis equal tight;
xlabel('x (m)');
ylabel('z (m)');
title('Magnetic field lines');
colormap(jet);   % 可选，也可单独设置

% 子图2：|B| 颜色云图（截断上限）
subplot(1,2,2);
% 设置截断阈值（例如取 Bmag 的 95% 分位数或手动给定）
threshold = prctile(Bmag(:), 95);  % 自动取 95% 分位数
% 或手动指定：threshold = 0.02;  % 单位特斯拉
Bmag_clipped = min(Bmag, threshold);

pcolor(X_grid, Z_grid, Bmag_clipped);
shading interp;
colorbar;
clim([0, threshold]);  % 强制颜色轴范围为 [0, threshold]
xlabel('x (m)');
ylabel('z (m)');
title(['|B| distribution (clipped at ', num2str(threshold, '%.2e'), ' T)']);
% 可选：整体加一个总标题
sgtitle(['Magnetic field of solenoid (R=', num2str(R), 'm, L=', num2str(L), 'm)']);

end