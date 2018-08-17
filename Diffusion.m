%% Initialize
%N-# of partc
nx=100;               %Number of steps in space(x)
nt=1000;               %Number of time steps 
dt=0.01;              %Width of each time step
t=0:dt:dt*nt;
length=1;
Gconst=5;
dx=length/(nx-1);         %Width of space step
x=0:dx:length;           
eps=2*dx;
D=0.01; 
temp=0;
w = zeros(nx,nt);
u = zeros(nx,nt);
V=zeros(1,nx); 

u(:,1)=0;
u(25,1)=1;

temp=0;
V(:)=1/nx;

for i=1:nx
    w(i,1)=u(i,1)*V(i);
end

%% compute distance matrix
dist=zeros(nx,nx);
for p=1:nx
    for q=1:nx
        if(abs(q-p)<nx-Gconst)
            dist(p,q)=abs(x(p) - x(q));
        else 
            dist(p,q)=length +0.01- abs(x(q) - x(p));
        end

    end
end

%% Solve PDE
dist
temp=0;
vec=[];newvec=[];
for i=1:nt-1
    for p=1:nx 
        a1=mod((p-1-Gconst),nx)+1;
        b1=mod(p-1+Gconst,nx)+1;
        %vec=  [ vec; p ,        a1         b1];
        if(a1<b1)
         for q=a1:b1
            temp=temp+(w(q,i)-w(p,i))*normpdf(dist(p,q),0,2*eps);%(w(q,i)-w(p,i))*

         end

        else

        for q=a1:nx

            temp=temp+(w(q,i)-w(p,i))*normpdf(dist(p,q),0,2*eps);%(w(q,i)-w(p,i))*

         end

         for q=1:b1

            temp=temp+(w(q,i)-w(p,i))*normpdf(dist(p,q),0,2*eps);%(w(q,i)-w(p,i))*

         end

            end

  
    %    newvec=[newvec temp];

    w(p,i+1)=w(p,i)+dt*V(p)*D/eps^2*temp;

    temp=0; 

    end

end

%% Substitute solution

temp=0;

    for i=1:nt

        for p=1:nx

        a1=mod((p-1-Gconst),nx)+1;

        b1=mod(p-1+Gconst,nx)+1;

        if(a1<b1)

         for q=a1:b1

             u(p,i)=u(p,i)+w(q,i)*normpdf(dist(p,q),0,2*eps);               

         end

        else

        for q=a1:nx

            u(p,i)=u(p,i)+w(q,i)*normpdf(dist(p,q),0,2*eps);               

         end

         for q=1:b1

      u(p,i)=u(p,i)+w(q,i)*normpdf(dist(p,q),0,2*eps);               

         end

        end

        end

    end

   yexact = zeros(nx,nt);

  
%% Plot
    for i=1:nt 
        h=plot(x,u(:,i));   
        axis([0 1 0 0.1])
        title({['1-D Diffusion with \nu =',num2str(D),' time(\itt) = ',num2str(i)]})
        drawnow; 
        refreshdata(h)
    end