
index=find(simIp.time == -0.688);
index=2432;
t=simIp.time(index);
p=Input_struct.p;
e=Input_struct.e;
Ip=simIp.signals.values(index);
IPF=simPFC.signals.values(index,:);
IIC= simIIC.signals.values(index,:);  
psi_sens=simFluxSens.signals.values(index,:);
B_sens=simMagnSens.signals.values(index,:);

[psicntrl1,psib1,r_cntrl_pnts1,z_cntrl_pnts1,status]= ...
    getCCSflux(t,p,e,Ip,IPF,IIC,psi_sens,B_sens,r_cntrl_gap_eq,z_cntrl_gap_eq);