% Assemble C matrix: PFC -> segments, X point field, Ip

% 1)Extract the C matrix for the desired outputs and call it C_XSC_full; 
%i.e. :


idxRxZx = signalIndexByName({'CV-RX','CV-ZX'},OutputNames, 1:length(OutputNames));
idxPsiB = signalIndexByName({'CV-FLUX'},OutputNames, 1:length(OutputNames));
idxPsiC = signalIndexByName(FluxCntrlNames, OutputNames, 1:length(OutputNames));

nPFC=10;


C_RxZx = C_contr(idxRxZx, 1:end);
C_PsiB = C_contr(idxPsiB, 1:end);
C_PsiC = C_contr(idxPsiC, 1:end);

C_PsiD = C_PsiC(2:end,:)/2/pi - repmat(C_PsiB, nOfFluxCntrl-1, 1); % remove the 1st one

C_XSC_full = [C_PsiD; C_RxZx];
C_XSC = C_XSC_full(:,1:nPFC);
nXSC = size(C_XSC,1); % n. of controlled variables


idxPFState = 1:nPFC;
C_XSC_Ip = -LinearModel.L(LinearModel.PlasmaCurrentInfo.StatePosition,idxPFState)/LinearModel.L(LinearModel.PlasmaCurrentInfo.StatePosition,LinearModel.PlasmaCurrentInfo.StatePosition);

wIp   = 1e-4;
% wPsi  = 2e0*ones(1,size(C_XSC,1));
% wPsi(end-1:end)=5e1*[1,1];
 if nOfFluxCntrl == 9
         wPsi = 2* [
        6.4616
        8.5176 
        9.7930 
        9.6720    
        8.7380    
        8.4285    
        9.1079   
        10.6375
        NaN
         NaN
        ]';
 elseif nOfFluxCntrl == 7
         wPsi = 2* [
        6.4616
        8.5176 
        9.7930 
        9.6720    
        8.7380    
        8.4285    
        NaN
         NaN
        ]';
 else
     wPsi = 2* [6.4616
            8.5176 
            9.7930 
            9.6720    
            8.7380    
            8.4285    
            9.1079   
            10.6375   
            10.2738    
            9.5358    
            9.5547    
            8.4208
            6.2867    
            4.4880    
            2.4693    
            1.4569    
            1.4295    
            1.0000    
            2.0401
            NaN
             NaN]';
 end


wPsi(end-1:end) = [20 20];

% wPsi  = 2e0*ones(1,size(C_XSC,1));
% wPsi(end-1:end)=5e1*[1,1];


C_XSC = diag([wPsi wIp])*[C_XSC; C_XSC_Ip];

% Singular Value Decomposition (svd(W * C)))
[U, S, V] = svd(C_XSC, 'econ');
SS=S;
UU=U;
VV=V;

% Compute the pinv and add weights
XSC_matrix = VV*inv(SS)*UU';
XSC_matrix = XSC_matrix(:, 1:end-1);

nDOF = 6;

XSC_matrix = VV(:,1:nDOF)*inv(SS(1:nDOF,1:nDOF))*UU(:,1:nDOF)';
XSC_matrix = XSC_matrix(:, 1:end-1);



save('XSC_SOF_flux.mat','XSC_matrix')


% eq values for test
eqPsiCntr = LinearModel.YEquil(FluxCntrlIdx)/(2*pi);
eqPsiX = equilRxZx(3);
eqRxZx = equilRxZx(1:2);

eqXSC = [eqPsiCntr(2:end)-eqPsiX; eqRxZx];

D_XSC=zeros();



%% Init PI


XSC_PID.Kp =0.1* ones(1, nPF);  
XSC_PID.Ki = 200*XSC_PID.Kp;
XSC_PID.Kd = zeros(1, nPF);
XSC_PID.Td = ones(1, nPF); % dummy value

% % % XSC PID controller
% XSC_PID.Kp = 11 * ones(1, nPF);  
% XSC_PID.Ki = 450 * ones(1, nPF);
% XSC_PID.Kd = 29*ones(1, nPF);
% XSC_PID.Td = (1/2000)*ones(1, nPF); % dummy value


XSC_PI = buildPID(XSC_PID);
return

