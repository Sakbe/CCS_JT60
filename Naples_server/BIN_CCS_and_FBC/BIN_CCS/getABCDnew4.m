function [A_contr, B_contr, C_contr, D_contr, varargout] = getABCD(linearModel, inputNames, outputNames, stabFlag, distFlag, cdNames)
% * [a b c d [inputList] [equilVoltages] [equilCurrents]] = getABCD(linearModel, inputNames, ...
% outputNames, [stabFlag], [distFlag], [cdNames]):
% returns the A, B, C, D matrices of the state space model.
% If distFlag is set to 1, the model accepts also Bp and Li as
% inputs.
% If the stabFlag is set to 1, the plasma is artificially
% stabilized reversing the sign of the unstable eigenvalue.
% If a cell array cdNames of current driven circuits is
% specified, the matrices of the model are modified to
% include such current driven circuits.



if nargin < 4, stabFlag = 0; end
if nargin < 5, distFlag = 0; end


%% Generate model matrices

%%%%%%% Only voltage driven circuits
if nargin < 6
  indexIn  = signalIndexByName(inputNames, linearModel.InputsInfo.Name, (1 : length(linearModel.InputsInfo.Name)));
  indexOut = signalIndexByName(outputNames, linearModel.OutputsInfo.Name, linearModel.OutputsInfo.OutputPosition);
  
  A_contr = linearModel.A;
  B_contr = linearModel.B(:,        indexIn);
  C_contr = linearModel.C(indexOut, :);
  D_contr = linearModel.D(indexOut, indexIn);
  
  % state indexes of the voltage driven circuits (needed for the disturbances matrices generation)
  indexVD  = linearModel.PoloidalCircuits.StatePosition;
  % state indexes of the eddy currents
  indexEddy = linearModel.EddyCurrentsInfo.StatePosition;
  
  % Plasma current
  if linearModel.PlasmaCurrentInfo.Kplasm
    if isempty(find(strcmp(inputNames,'Vpl'),1))
      indexEddy = [indexEddy, linearModel.PlasmaCurrentInfo.StatePosition];
    else
      indexVD = [indexVD, linearModel.PlasmaCurrentInfo.StatePosition];
    end
  end
  
  % Ip0 = linearModel.XEquil(linearModel.PlasmaCurrentInfo.StatePosition);
  
  
  %%%%%%% Current driven circuit state transformation
else
  
  vdNames = inputNames;
  if(iscell(cdNames) == 0)
    cdNames = {cdNames};
  end
  if(iscell(vdNames) == 0)
    vdNames = {vdNames};
  end
  inputNames = {vdNames{1:end}, cdNames{1:end}}; % New inputNames vector
  
  indexIn  = signalIndexByName(vdNames,linearModel.InputsInfo.Name,(1:length(linearModel.InputsInfo.Name)));
  indexOut = signalIndexByName(outputNames,linearModel.OutputsInfo.Name,linearModel.OutputsInfo.OutputPosition);
  
  % State indexes of the CD, VD, and passive circuits
  indexCD = signalIndexByName(cdNames, linearModel.PoloidalCircuits.Name, linearModel.PoloidalCircuits.StatePosition);
  indexVD = setdiff(linearModel.PoloidalCircuits.StatePosition, indexCD);
  indexEddy = linearModel.EddyCurrentsInfo.StatePosition;
  
  if linearModel.PlasmaCurrentInfo.Kplasm
    if     ~isempty(find(strcmp(cdNames,'Ipl'),1)),  indexCD = [indexCD, linearModel.PlasmaCurrentInfo.StatePosition];
    elseif ~isempty(find(strcmp(vdNames,'Vpl'),1)),  indexVD = [indexVD, linearModel.PlasmaCurrentInfo.StatePosition];
    else,                                            indexEddy = [indexEddy, linearModel.PlasmaCurrentInfo.StatePosition];
    end
  end
  
  %%% State transformation
  % 1 = VD
  % 2 = passive circuits
  % 3 = CD
  
  % Modified L matrices
  L11 = linearModel.L(indexVD,   indexVD);
  L12 = linearModel.L(indexVD,   indexEddy);
  L21 = linearModel.L(indexEddy, indexVD);
  L22 = linearModel.L(indexEddy, indexEddy);
  L1 = [L11, L12; L21, L22];
  
  L13 = linearModel.L(indexVD,   indexCD);
  L23 = linearModel.L(indexEddy, indexCD);
  L2 = [L13; L23];
  
  % Modified R matrices
  R1 = linearModel.R(indexVD,   indexVD);
  R2 = linearModel.R(indexEddy, indexEddy);
  R = [R1, zeros(size(R1, 1), size(R2, 2)); zeros(size(R2, 1), size(R1, 2)), R2];
  
  % Base Transformation matrices
  T1 = inv(L1);
  T2 = -T1 * L2;
  
  % Modified state-space model matrices
  A_contr = -R * T1;
  
  B1 = [linearModel.S(indexIn,indexIn); zeros(length(indexEddy), length(indexIn))];
  %   indexPoloidalCircuits = [linearModel.PoloidalCircuits.StatePosition linearModel.PlasmaCurrentInfo.StatePosition];
  %   indexB1 = find(ismember(indexPoloidalCircuits, indexVD)); % remove circuits that are not in the input vector from the B columns
  %   indexB1 = find(ismember(indexPoloidalCircuits(indexB1), indexIn));
  %   B1 = B1(:, indexB1);
  B2 = -R * T2;
  B_contr = [B1 B2];
  
  C1 = linearModel.C(indexOut, indexVD);
  C2 = linearModel.C(indexOut, indexEddy);
  C3 = linearModel.C(indexOut, indexCD);
  C_contr = [C1 C2] * T1;
  
  %   D1 = [linearModel.D(indexOut(1 : (end-1)), indexIn); zeros(1, length(indexIn))]; % indexOut now contains also lmszx1
  D1 = linearModel.D(indexOut(1 : (end)), indexIn);
  D2 = [C1 C2] * T2 + C3;
  D_contr = [D1 D2];
  
  
  
end

%% Disturbances

%== Old version ===========================================================
%   if distFlag
%     Edot = -linearModel.LE;
%     E = zeros(size(linearModel.LE));
%     Fdot = linearModel.Fdot;
%     F = linearModel.F;
%   else
%     Edot = zeros(size(linearModel.LE));
%     E = Edot;
%     Fdot = zeros(size(linearModel.Fdot));
%     F = Fdot;
%   end
%   B_contr = [B_contr, E([indexVD, indexEddy], :), Edot([indexVD, indexEddy], :)];
%   D3 = [F(indexOut(1 : end-1), :), Fdot(indexOut(1 : end-1), :); zeros(1, 4)]; % for lmsxz1
%   D_contr = [D_contr, D3];
%==========================================================================


if distFlag
  if ~linearModel.PlasmaCurrentInfo.Kplasm
    error('Cannot apply disturbances to a plasmaless model!');
  end
  
  if nargin >= 6 % CD coils are specified
    % Substitute E with -LE!
    B_contr = [B_contr -A_contr*linearModel.LE([indexVD, indexEddy], :)];
    D_contr = [D_contr -C_contr*linearModel.LE([indexVD, indexEddy], :) + linearModel.F(indexOut, :)];
  else
    % Use x - Ew as a new state variable; Fdot is supposed
    % to be the zero matrix
    B_contr = [B_contr A_contr*linearModel.E([indexVD, indexEddy], :)];
    D_contr = [D_contr C_contr*linearModel.E([indexVD, indexEddy], :)+linearModel.F(indexOut, :)];
  end
else
  B_contr = [B_contr zeros(size(B_contr, 1), 2)];
  D_contr = [D_contr zeros(size(C_contr, 1), 2)];
end
inputNames{end+1} = 'betapol*Ipl';
inputNames{end+1} = 'li*Ipl';

% Remove NaNs;
A_contr(isnan(A_contr)) = 0;
B_contr(isnan(B_contr)) = 0;
C_contr(isnan(C_contr)) = 0;
D_contr(isnan(D_contr)) = 0;



%% Numerically stabilize the unstable plasma

if stabFlag % stabFlag is the 4th parameter
  [V, D] = eig(A_contr);
  % Index of the unstable mode
  u_index = find(diag(D) == max(max(real(D)))); % if the eigenvalues are all negative, max(D) = 0 and u_index is empty
  % Reverse the unstable mode
  D(u_index, u_index) = -D(u_index, u_index);
  A_contr = V * D / V;
  A_contr = real(A_contr);
end


%% Assign optional outputs
if nargout >= 5
  varargout{1} = inputNames;
  if nargout >= 6
    varargout{2} = linearModel.R(indexVD,indexVD) * linearModel.XEquil(indexVD);
    if nargout ==7
      varargout{3} = linearModel.XEquil(indexCD);
    end
  end
end

end



function index = signalIndexByName(name,nameVector,indexVector)
% function index = signalIndexByName(name,nameVector,indexVector)
% If name is a cell array of strings returns an vector of indexes
% If a name is present more than once in the nameVector only the first index will be returned
index = [];
if iscell(name)
  for i = 1:length(name)
    temp = indexVector(find(strcmp(nameVector,name{i})==1));
    if ~isempty(temp)
      index(end+1) = temp(1);
    end
  end
else
  temp = indexVector(find(strcmp(nameVector,name)==1));
  index = temp(1);
end
end