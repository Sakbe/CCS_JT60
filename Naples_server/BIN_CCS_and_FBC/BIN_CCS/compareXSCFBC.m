close all

%%% Compare XSC and FBC
clear CntrlPFC
for i=1:10
    figure(i)
 plot(simPFC.time,simPFC.signals.values(:,i))
 hold on
 plot(simIPFc.time,simIPFc.signals.values(:,i))
 newIPFc=interp1(simIPFc.time,simIPFc.signals.values,simPFC.time);
 CntrlPFC(:,i)=simPFC.signals.values(:,i)-newIPFc(:,i);
end
%%
close all

for i=1:10
    figure(i)
plot(CntrlPFC(:,i))
grid on
end
%%
close all
for i=1:10
    figure(i)
     plot(simIPFc.time,simIPFc.signals.values(:,i))
hold on
    plot(simCurrCont.time,simCurrCont.signals.values(:,i))
    grid on
end
%%
for i=1:10
    figure(i)
 plot(simVPF.time,simVPF.signals.values(:,i))
 hold on
 plot(simVPFc.time,simVPFc.signals.values(:,i))
end
%% simFluxCntrl
for i=2:9
    figure(i)
 plot(simFluxCntrl.time,simFluxCntrl.signals.values(:,i))
 hold on
 plot(psc_jap.time,psc_jap.signals.values(:,i))
end