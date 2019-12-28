%%%%%%%%%%%%%%% RE-write JT60-SA model
inputNames={'CS1','CS2','CS3','CS4','EF1','EF2','EF3','EF4','EF5','EF6','Vpl'};
%inputNames={'CS1','CS2','CS3','CS4','EF1','EF2','EF3','EF4','EF5','EF6','VSU','VSL','Vpl'}; 
outputNames1=OutputNames;
% outputNames1(11)={'passive1'};
% outputNames1(12)={'passive2'};
[A_contr1, B_c, C_contr1, D_c, varargout] = getABCDnew1(LinearModel, inputNames, outputNames1,[0], [1],[{'VSU','VSL'}]);
%[A_contr1, B_c, C_contr1, D_c, varargout] = getABCDnew(LinearModel, inputNames, outputNames1,[0], [1],[]);


%% Check models are okey 
xii=zeros(1,140);
equilVolts1=zeros(13,1);

inputNames={'CS1','CS2','CS3','CS4','EF1','EF2','EF3','EF4','EF5','EF6','Vpl'};
outputNames1=OutputNames;
[A_contr1, B_contr1, C_contr1, D_contr1, varargout] = getABCDnew4(LinearModel, inputNames, outputNames1,[1], [1],[{'VSU','VSL'}]);

inputNames={'CS1','CS2','CS3','CS4','EF1','EF2','EF3','EF4','EF5','EF6','VSU','VSL','Vpl'};
[A_contr2, B_contr2, C_contr2, D_contr2, varargout] = getABCDnew4(LinearModel, inputNames, outputNames1,[1], [1],[]);
%[A_contr2, B_contr2, C_contr2, D_contr2, varargout] = getABCDnew4(LinearModel, inputNames, outputNames1,0,1);

%% Check with the class
clc
object=CREATE2DModel( '../../L-models_90deg/SOF@18d66s_japanese_sensors_8cntrlFluxPnts_CL.mat','../../L-models_90deg/SOF@18d66s_japanese_sensors_8cntrlFluxPnts.mat')

inputNames={'CS1','CS2','CS3','CS4','EF1','EF2','EF3','EF4','EF5','EF6','Vpl'};
[ssModel1, Tt, stateIndex1] =object.getStateSpace( inputNames, OutputNames,[1], [1],[{'VSU','VSL'}]);

inputNames={'CS1','CS2','CS3','CS4','EF1','EF2','EF3','EF4','EF5','EF6','VSU','VSL','Vpl'};
[ssModel2, Tt, stateIndex2] =object.getStateSpace( inputNames, OutputNames,[1], [1]);

 A_contr1=ssModel1.A;
 B_contr1=ssModel1.B;
C_contr1=ssModel1.C;
 D_contr1=ssModel1.D;
 
 A_contr2=ssModel2.A;
 B_contr2=ssModel2.B;
C_contr2=ssModel2.C;
 D_contr2=ssModel2.D;
 
 %% Check Ariano stuff 3-oct
 clc
 
 