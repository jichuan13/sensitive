function err=Lcans(lambda)
k0=1000;
n=30;M=2;ro=0.3170; 
loc=[ -0.0069,0.0095 ];
g=BoundaryValue(n,0,0);
rand('seed',0)
g=(1+0.01*(2*rand(length(g),1)-1)).*g;
dth=2*pi/n; th=[0:dth:2*pi];
for jj=1:n+1 
   for kk=0:M
       tri(jj,kk+1)=cos(kk*th(jj));
      dtri(jj,kk+1)=-kk*sin(kk*th(jj));
    end
    for kk=1:M
       tri(jj,kk+M+1)=sin(kk*th(jj));
        dtri(jj,kk+M+1)=kk*cos(kk*th(jj));
    end
end
beta=zeros(2*M+1,1); 
beta(1)=ro; 
s=tri*beta;
ds=dtri*beta;


tol=norm(g);kk=1;og1=g;
while(kk<=k0&&tol>=1e-6)
    kk
beta0=beta;
g1=og1;
og=OriginalBoundaryValue(beta0,M,n,loc);
g2=og;
og1=g2;
GG=g-og;

x0=loc(1,1);y0=loc(1,2);
r=s';
rd=ds';
y1=x0+r.*cos(th);y2=y0+r.*sin(th);
y1d=rd.*cos(th)-r.*sin(th);y2d=rd.*sin(th)+r.*cos(th);
z=cos(th)'+1i*sin(th)';x=y1+1i*y2;  
for k=1:length(th); 
    G=(log(abs(z-x(k)))-log(abs(z./abs(z)-abs(z)*x(k))))/(2*pi);
    GGG1(:,k)=G;
end
z=(1+1.e-7)*z;
for k=1:length(th); 
    G=(log(abs(z-x(k)))-log(abs(z./abs(z)-abs(z)*x(k))))/(2*pi);
    GGG2(:,k)=G; 
end
F=(GGG2-GGG1)/1.e-7;
gg=F*diag(y2d.*cos(th)-y1d.*sin(th))*tri*dth;
dg=GG'*gg*dth;
dF=-dg;
J=gg;R=-GG;
%阻尼最小二乘法
% lambda=0.002;
beta1=beta0'-(inv(J'*J+lambda*eye(length(beta0)))*(J'*R))';
%高斯牛顿法
%  beta1=beta0'-(inv(J'*J)*(J'*R))';
%最速下降法
% alpha=0.4;
%  beta1=beta0'-alpha*((J'*R))';
beta=beta1';

beta';
tol=norm(g2-g1);

s=tri*beta;
ds=dtri*beta;



kk=kk+1;
end

a=0:2*pi/n:2*pi;
x0=loc(1,1);y0=loc(1,2);
r=s';
x1=x0+r.*cos(a);y1=y0+r.*sin(a);
hold on 
plot(x1,y1,'b')
r=0.5*(0.8.*sqrt((cos(a)).^2+0.25.*(sin(a)).^2));
x4=-0.0+r.*cos(a);y4=0.0+r.*sin(a);
plot(x4,y4,'r-')
hold on
ex=x1-x4;ey=y1-y4;
err=sum(sqrt(ex.^2+ey.^2))

% figure(2)
% % plot(lam,e)
% plot(e)



