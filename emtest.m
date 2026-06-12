function B=emtest(R,I,L,d,N)

poss=[0;0;L/2];
theta=linspace(0,2*pi*L/d,N);
pos=[R.*cos(theta);R.*sin(theta);theta.*(d/2/pi)];
dtheta=theta(2)-theta(1);
dl=dtheta.*[-R.*sin(theta);R.*cos(theta);ones(1,N).*(d/2/pi)];
B=zeros(3,1);
for j=1:(N)
    tpos=[pos(1,j);pos(2,j);pos(3,j)];
    r=poss-tpos;
    tB=I.*cross(dl(:,j),r)./(norm(r)^3);
    B=B+tB;
end
B=B*10^(-7);

end