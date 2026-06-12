function Ball=emmidfuncd(R,I0,L0,dL,dd,d0,N0,N1,x_func)

Ball=zeros(1,N1);
for i=1:N1
    L=L0+i*dL;
    d=d0-i*dd;
    I=I0*d;

    x=x_func(L,R);
    N=floor(N0*L/d/2);
    poss=[x;0;0];
    theta1=linspace(0,pi*L/d,N);
    theta2=linspace(-pi*L/d,0,N);
    theta=[theta1,theta2];
    pos=[R.*cos(theta);R.*sin(theta);theta.*(d/2/pi)];
    dtheta=theta(2)-theta(1);
    dl=dtheta.*[-R.*sin(theta);R.*cos(theta);ones(1,2*N).*(d/2/pi)];
    B=zeros(3,1);
    for j=1:2*N
        tpos=[pos(1,j);pos(2,j);pos(3,j)];
        r=poss-tpos;
        tB=I.*cross(dl(:,j),r)./(norm(r)^3);
        B=B+tB;
    end
    Ball(1,i)=norm(B);
end

plot(L0+[1:N1]*dL,Ball);

end