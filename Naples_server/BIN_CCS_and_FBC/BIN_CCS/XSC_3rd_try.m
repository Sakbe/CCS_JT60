% Assemble C matrix: PFC -> segments, X point field, Ip

% 1)Extract the C matrix for the desired outputs and call it C_XSC_full; 
%i.e. :

ind1=gapIdx(1);
%ind2=562; only 83 gaps
ind2=gapIdx(end); %83 gaps +2 gaps in strikes

nPFC=10;
idxPFState = 1:nPFC;
% % % %% All gaps
% C_XSC_gaps = LinearModel.C(ind1:1:ind2,idxPFState);

% % 30 gaps
% C_XSC_gaps = LinearModel.C([ind1:3:(ind2-2),ind2-1,ind2],idxPFState);
% 
% % %% 19 gaps
% C_XSC_gaps = LinearModel.C([ind1:5:(ind2-2),ind2-1,ind2],idxPFState);
% 
% % 20 gaps
C_XSC_gaps = LinearModel.C([round((ind1+1:4.6:(ind2-2))),ind2-1,ind2],idxPFState);


% 
% % % %Japanese selection 1
% japs_gaps=[13,22,36,52,66,76,84,85]';
% japs_gaps=japs_gaps+ind1-1;
% C_XSC_gaps = LinearModel.C(japs_gaps,idxPFState);

% % % % 
% %%%Japanese selection 2
% japs_gaps=[21,36,47,66,84,85]';
% japs_gaps=japs_gaps+ind1-1;
% C_XSC_gaps = LinearModel.C(japs_gaps,idxPFState);

C_XSC_Ip   = -LinearModel.L(LinearModel.PlasmaCurrentInfo.StatePosition,idxPFState)/LinearModel.L(LinearModel.PlasmaCurrentInfo.StatePosition,LinearModel.PlasmaCurrentInfo.StatePosition);


% Weights for the gaps and the plasma current
wIp   = 1e-5;
wGaps = 10*ones(1,size(C_XSC_gaps,1));
%wGaps([15:36])=0.10*wGaps([15:36]);
% wGaps([60:70]) = 0.05*wGaps([60:70]);

C_XSC = diag([wGaps wIp])*[C_XSC_gaps(:,idxPFState);C_XSC_Ip];

% Singular Value Decomposition (svd(W * C)))
% [U, S, V] = svd(wXSC * C_XSC * wPFC, 'econ');
[U, S, V] = svd(C_XSC, 'econ');
% Remove some dofs
% nOfDofsToBeRemoved = 0;
% nXSC=length(C_XSC(:,1))-1;
% nDOF = nXSC + 1 - nOfDofsToBeRemoved;
% SS = S(1:nDOF,1:nDOF);
% UU = U(:,1:nDOF);
% VV = V(:,1:nDOF);
SS=S;
UU=U;
VV=V;

% Compute the pinv and add weights
XSC_matrix = VV*inv(SS)*UU';
XSC_matrix = XSC_matrix(:, 1:end-1);

nDOF = 7;

XSC_matrix = VV(:,1:nDOF)*inv(SS(1:nDOF,1:nDOF))*UU(:,1:nDOF)';
XSC_matrix = XSC_matrix(:, 1:end-1);


save('XSC_SOF.mat','XSC_matrix')
