close all
clear all

% Equilibrium value for Ip, beta and li - JT60-SA Scenario 2 at SOF
%
Ip0 = 5.5e6;
betap0 =  0.5298;
li0 = 0.8538;

% Minor Disruption
% betap drop about 70% (according to ITER plasmadyn disturbance)
betap = [0 0; 0.1 0; 0.1+1e-3 -0.53*betap0; 0.2 -0.53*betap0];
betap(:,2) = betap0+betap(:,2);
% li drop about 12.5% (according to ITER plasmadyn disturbance)
li = [0 0; 0.1 0; 0.1+1e-3 -0.125*li0; 0.2 -0.125*li0];
li(:,2) = li0+li(:,2);

figure(1)
suptitle('Minor disruption ITER-like')
subplot(2,1,1)
plot(betap(:,1)+20,betap(:,2))
grid on
xlabel('Time [s]')
ylabel('\beta_p');
subplot(2,1,2)
plot(li(:,1)+20,li(:,2),'r')
ylabel('l_i');
grid on
xlabel('Time [s]')