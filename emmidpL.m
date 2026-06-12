function Ball=emmidpL(R,I0,L,dd,d0,N,N1,x,y)

Ball=zeros(1,N1);
for i=1:N1
    d=d0-i*dd;
    I=I0*d;
    poss=[x;y;L/2];
    N0=floor(N*L/d);
    theta=linspace(0,2*pi*L/d,N0);
    pos=[R.*cos(theta);R.*sin(theta);theta.*(d/2/pi)];
    dtheta=theta(2)-theta(1);
    dl=dtheta.*[-R.*sin(theta);R.*cos(theta);ones(1,N0).*(d/2/pi)];
    B=zeros(3,1);
    for j=1:(N0)
        tpos=[pos(1,j);pos(2,j);pos(3,j)];
        r=poss-tpos;
        tB=I.*cross(dl(:,j),r)./(norm(r)^3);
        B=B+tB;
    end
    Ball(1,i)=norm(B);
end

plot(d0-[1:N1]*dd,Ball);

end