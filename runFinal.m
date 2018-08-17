
N_A=2.4*10^(5);
N_P= 9.8*10^(4);
x0 = [0 N_A];
t=0:dt:dt*nt;
dt=0.01; 
nt=50;
[t,y] = ode45(@reaction_sysFINAL,t, x0);

plot(t,y(:,1),t,y(:,2))
legend(' y: A','y: P');
y(end,1)
y(end,2)


