function dydt = reaction_sysFINAL(t,x)
V_cyto= 2.5*10^(4);
S_cyto= 4.4*10^(3);
N_A=2.4*10^(5);
N_P= 9.8*10^(4);
k_offA=3.24*10^(-3);
k_offP=7.19*10^(-3);
k_onA=6.29*10^(-3)*S_cyto/V_cyto;
k_onP=7.682*10^(-2);
alpha = 2;
beta = 2;
k_AP = 0.01*(S_cyto/V_cyto)^(2 + 2*alpha);%10^(-2)
k_PA = 0.01*(S_cyto/V_cyto)^(2 + 2*beta);%10^(-5)

dydt = [-(k_offA+k_onA)*x(1)-k_AP*x(1)*x(2)^alpha+k_onA*N_A; 
    -(k_offP+k_onP)*x(2)-k_PA*x(2)*x(1)^beta+k_onP*N_P;];

end
%c1=k_onA*N_A/V_cyto 
%c2=(k_offA + (k_onA*S_cyto)/V_cyto)