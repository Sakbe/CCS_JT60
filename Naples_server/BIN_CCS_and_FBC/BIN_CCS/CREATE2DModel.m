classdef CREATE2DModel < plasmaCircuit2DModel
    
    properties (Access = private)
        % Geometrical descriptors data 
        gapData = [];  
        
        % Flux sensors identifiers and positions contained in Input_struct
        namesSensors = [];
        rSensors = [];
        zSensors = [];
        tSensors = [];
        Kplasm = 1;

    end
    methods
        function obj = CREATE2DModel(varargin)
            
            if nargin > 2
                error('Too many input arguments');
            end
            
            if (isa(varargin{1},'char') == 1)
                if(exist(varargin{1},'file')==2 || exist([varargin{1} '.mat'],'file')==2)
                    model=load(varargin{1});
                    
                    obj.nOfInputs = model.LinearModel.NOfInputs;
                    obj.nOfOutputs = model.LinearModel.NOfOutputs;
                    obj.nOfStates = model.LinearModel.NOfStates;
                    obj.nOfDisturbances = model.LinearModel.NOfDisturbances;
                    obj.plasmalessFlag = model.LinearModel.PlasmaCurrentInfo.Kplasm;
                    obj.A = model.LinearModel.A;
                    obj.B = model.LinearModel.B;
                    obj.C = model.LinearModel.C;
                    obj.D = model.LinearModel.D;
                    obj.E = model.LinearModel.E;
                    obj.F = model.LinearModel.F;
                    obj.Fdot = model.LinearModel.Fdot;
                    obj.L = model.LinearModel.L;
                    obj.LE = model.LinearModel.LE;
                    obj.R = model.LinearModel.R;
                    obj.S = model.LinearModel.S;
                    obj.uEq = model.LinearModel.UEquil;
                    obj.wEq = model.LinearModel.WIpEquil;
                    obj.xEq = model.LinearModel.XEquil;
                    obj.yEq = model.LinearModel.YEquil;
                    obj.uNames = model.LinearModel.InputsInfo.Name;
                    obj.uPositions = 1:model.LinearModel.NOfInputs;
                    obj.xNames = [model.LinearModel.PoloidalCircuits.Name' model.LinearModel.EddyCurrentsInfo.Name' 'Ipl'];
                    obj.xPositions = [model.LinearModel.PoloidalCircuits.StatePosition model.LinearModel.EddyCurrentsInfo.StatePosition model.LinearModel.PlasmaCurrentInfo.StatePosition];%non credo funziona
                    obj.PFStatePosition = model.LinearModel.PoloidalCircuits.StatePosition;
                    obj.EddyStatePosition = model.LinearModel.EddyCurrentsInfo.StatePosition;
                    obj.IpStatePosition = model.LinearModel.PlasmaCurrentInfo.StatePosition;
                    obj.yNames = model.LinearModel.OutputsInfo.Name;
                    obj.yPositions = model.LinearModel.OutputsInfo.OutputPosition;
                    obj.wNames = model.LinearModel.DisturbancesInfo.Name;
                    obj.wPositions = model.LinearModel.DisturbancesInfo.OutputPosition;
                    obj.Kplasm = model.LinearModel.PlasmaCurrentInfo.Kplasm;
                    
                    if  obj.Kplasm == 0 % Plasmaless
                        obj.L = obj.L(1:end-1,1:end-1);
                        obj.R = obj.R(1:end-1,1:end-1);
                        obj.A = obj.A(1:end-1,1:end-1);
                        obj.B = obj.B(1:end-1,:);
                        obj.S = obj.S(1:end-1,:);
                        obj.C = obj.C(:,1:end-1);
                        obj.LE = obj.LE(1:end-1,:);
                    end
                    
                    %  Modify the model to take correctly into account the disturbances
%                     idxIp=obj.xPositions(end);
%                     idxI=setdiff(1:size(obj.L,1),idxIp);
%                     oldL=obj.L;
%                     oldC=obj.C;
%                     W0=obj.yEq(obj.wPositions)/obj.yEq(model.LinearModel.PlasmaCurrentInfo.OutputPosition);
%                     obj.L(idxI,idxIp)=obj.L(idxI,idxIp)+obj.LE(idxI,:)*W0;
%                     obj.L(idxIp,idxIp)=obj.L(idxIp,idxIp)+obj.LE(idxIp,:)*W0;
%                     obj.C(:,idxIp)=obj.C(:,idxIp)+obj.F*W0;
                    
                    
                    if nargin==2
                        if(exist(varargin{1},'file')==2 || exist([varargin{1} '.mat'],'file')==2)
                            load(varargin{2});
                            
                            obj.vesselData = struct('p', Input_struct.p, 'e', Input_struct.e);
                            obj.equilFluxData = struct('p',Input_struct.p, 't', Input_struct.t, ...
                                                       'nodeData', x_np(1:length(Input_struct.p)), ...
                                                       'boundaryFlux',y_np(strcmp(y_type,'psb_c')));

                            [obj.gapData.equilValues, obj.gapData.idx] = findGap(obj); % remove NaNs
                            obj.gapData.r = Input_struct.r_sens_gap(obj.gapData.idx);
                            obj.gapData.z = Input_struct.z_sens_gap(obj.gapData.idx);
                            obj.gapData.theta = Input_struct.theta_sens_gap_deg(obj.gapData.idx);
                            
                            obj.namesSensors = Input_struct.names_sensors; 
                            obj.rSensors = Input_struct.r_sens;
                            obj.zSensors = Input_struct.z_sens;
                            obj.tSensors = Input_struct.theta_sens;
                        end
                    end
                end
            elseif (isa(varargin{1},'CREATE2DModel') == 1)
                obj.nOfInputs = varargin{1}.nOfInputs;
                obj.nOfOutputs = varargin{1}.nOfOutputs;
                obj.nOfStates = varargin{1}.nOfStates;
                obj.nOfDisturbances = varargin{1}.nOfDisturbances;
                obj.plasmalessFlag = varargin{1}.plasmalessFlag;
                obj.A = varargin{1}.A;
                obj.B = varargin{1}.B;
                obj.C = varargin{1}.C;
                obj.D = varargin{1}.D;
                obj.E = varargin{1}.E;
                obj.F = varargin{1}.F;
                obj.Fdot = varargin{1}.Fdot;
                obj.L = varargin{1}.L;
                obj.LE = varargin{1}.LE;
                obj.R = varargin{1}.R;
                obj.S = varargin{1}.S;
                obj.uEq = varargin{1}.uEq;
                obj.wEq = varargin{1}.wEq;
                obj.xEq = varargin{1}.xEq;
                obj.yEq = varargin{1}.yEq;
                obj.uNames = varargin{1}.uNames;
                obj.uPositions = varargin{1}.uPositions;
                obj.xNames = varargin{1}.xNames;
                obj.xPositions = varargin{1}.xPositions;
                obj.PFStatePosition = varargin{1}.PFStatePosition;
                obj.EddyStatePosition = varargin{1}.EddyStatePosition;
                obj.IpStatePosition = varargin{1}.IpStatePosition;
                obj.yNames = varargin{1}.yNames;
                obj.yPositions = varargin{1}.yPositions;
                obj.wNames = varargin{1}.wNames;
                obj.wPositions = varargin{1}.wPositions;
                obj.vesselData = varargin{1}.vesselData;
                obj.equilFluxData = varargin{1}.equilFluxData;
                obj.gapData = varargin{1}.gapData;
                obj.namesSensors = varargin{1}.namesSensors;
                obj.rSensors = varargin{1}.rSensors;
                obj.zSensors = varargin{1}.zSensors;
            end
        end
                   
        % Plot the 2D plasma shape from the gaps values
        % AH = plotPlasmaShape(OBJ,GAPVALUES,XPOS,[GAPFLAG],[LINESTYLE],[INTERPMETHOD])
        % GAPVALUES    are the values of the gaps
        % XPOS         is the null point position specified as [rx, zx]
        % GAPFLAG      specifies whether to plot the gaps or no       
        % LINESTYLE    specifies the line style
        % INTERPMETHOD specifies interpolation method (linear, spline, pchip). See interp1 for details.  
        function plotPlasmaShape(obj,gapValues,xPos,varargin)
            if ~isempty(obj.gapData)
                
                obj.plotMechanicalStructure();
                hold on

                % Default parameters values
                gapFlag = 0;
                lineStyle = '-om';
                interpMethod = 'linear';
                
                if nargin > 3
                    gapFlag = varargin{1};
                    if nargin > 4
                        lineStyle = varargin{2};
                        if nargin > 5
                            interpMethod = varargin{3};
                            if nargin > 6
                                error('Too many input variables.');
                            end
                        end
                    end
                end
            end
        
            % Some manipulations
            [~,gapPlotIdx]= obj.findGap; % returns the equilibrium values of the gaps and the non-NaN indexes of the gaps
            gapValues = gapValues(gapPlotIdx);
            if isrow(gapValues)
                gapValues = gapValues';
            end
            
            % Plot gaps if requested
            outIdx = find(gapValues > 2*mean(gapValues)); 
            plotIdx = setdiff(1:length(gapValues), outIdx);
            if~isempty(outIdx)
                warning('Outliers gaps found. It is possible that some gaps do not intersect the plasma boundary.');
            end
            if(gapFlag)
                gapLength = .5 * mean(gapValues);
                rr = [obj.gapData.r, obj.gapData.r + gapLength*cosd(obj.gapData.theta)];
                zz = [obj.gapData.z, obj.gapData.z + gapLength*sind(obj.gapData.theta)];
                for i = plotIdx'
                    plot(rr(i,:), zz(i,:), '-xk')
                end
                for i = outIdx'
                    plot(rr(i,:), zz(i,:), '-x', 'Color', [.3 .3 .3])
                end
            end
                        
            % Compute boundary points coordinates from gaps   
            R = obj.gapData.r + gapValues.*cosd(obj.gapData.theta);
            Z = obj.gapData.z + gapValues.*sind(obj.gapData.theta);

            % Add X-point
            rx = xPos(1);
            zx = xPos(2);
            if ~any(isnan(xPos))
                R = [rx; R];
                Z = [zx; Z];
            end  
            
            % Find center
            rc = sum(R)/length(R);
            zc = sum(Z)/length(Z);
                       
            % Cartesian to polar coordinates
            [gapth, gaprho] = cart2pol(R-rc,Z-zc);
            [gapth,i,~] = unique(gapth);
            gaprho = gaprho(i);
            
            % Sort gaps
            if ~any(isnan(xPos))
                xth = atan2(zx-zc, rx-rc);
                gapth(gapth<xth) = gapth(gapth<xth) + 2*pi;
            end
            [~, i_sort] = sort(gapth);
            
            gapth = gapth(i_sort);
            gaprho = gaprho(i_sort);
            
            % Polar to cartesian
%             gapth = [gapth; gapth(1)+2*pi];
%             gaprho = [gaprho; gaprho(1)];
            [newR,newZ] = pol2cart(gapth,gaprho);
            
            % Assemble for plot
            rNew = newR+rc;
            zNew = newZ+zc;
            
            % Separate strike points from boundary points
            if(zx <= zc) % Lower Single Null 
                strikeGaps = find(zNew < zx);                
            else % Upper Single Null
                strikeGaps = find(zNew > zx);
            end
            boundGaps = setdiff(1:length(zNew),strikeGaps);
            
            % Close line 
            rBound = rNew(boundGaps([1:length(boundGaps), 1]));
            zBound = zNew(boundGaps([1:length(boundGaps), 1]));
            
            % Search for ill-defined gaps
            iOutliers = find(abs(zBound)>2*mean(abs(zBound)));
            iOk = setdiff(1:length(zBound), iOutliers);
            
            % Plot
            plot(rBound(iOk), zBound(iOk), lineStyle)
            plot(rBound(iOutliers), zBound(iOutliers), 'pr')
            
            hold on
            
            % Strike points
            RS = rNew(strikeGaps);
            ZS = zNew(strikeGaps);
            [ri, iMin] = min(RS);
            [rf, iMax] = max(RS);
            zi = ZS(iMin);
            zf = ZS(iMax);
            
            if (ri<rx && rf>rx)
                rStrike = [ri; rx; rf];
                zStrike = [zi; zx; zf];
            else
                rStrike = [];
                zStrike = [];
            end
            
            plot(rStrike, zStrike, lineStyle)
            hold off
        end           
       
        % Used by plot shape (make it private?)
        function [gaps, idxGap] = findGap(obj)     
            temp = find(strncmp('GAP',obj.yNames,3));
            temp = temp(1:end/2); % remove GAP*Ipl
            idxNan = find(isnan(obj.yEq(temp)));

            idxGap = setdiff(1:length(temp), idxNan);
            gaps = obj.yEq(temp(idxGap));
                        
        end

        % Returns the number of gaps
        function nGaps = getNOfGaps(obj)
            nGaps = length(obj.gapData.equilValues);
        end
        
        % Returns the gap IDs
        function gapIDs = getGapIDs(obj)
            gapIDs = obj.yNames(strncmp(obj.yNames, 'GAP', 3));
            gapIDs = gapIDs(1:end/2);
        end
        
        % Returns the number of points per control segment 
        % NOTA: I SEGMENTI POTREBBERO NON CHIAMARSI SEMPRE Fluxv!
        function nPoints = getNOfPoints(obj, varargin)
            if nargin == 1
                nPoints = length(find(strncmp(obj.yNames, 'Fluxv_01_', 9)));
            else
                nPoints = [];
                names = varargin{1};
                for i = 1 : length(names)
                    nPoints(end+1,:) = length(find(strncmp(obj.yNames, [names{i} '_'], length(names{i})+1)));
                end
            end
        end
        
        % Returns the names and the position of the magnetic probes
        function [n,r,z,t] = getSensorsInfo(obj)
            n = obj.namesSensors;
            r = obj.rSensors;
            z = obj.zSensors;
            t = obj.tSensors;
        end
        
        % Returns equilibrium flux map
        function eqMap = getEquilFluxMap(obj)
            eqMap = obj.equilFluxData;
        end
        
        % Return the number of segments, the number of points per control
        % segments and the index of the segments in the output vector of
        % the CREATE model
        % NOTA: I SEGMENTI POTREBBERO NON CHIAMARSI SEMPRE Fluxv!
        function [nSeg, nPoints, idxOut] = getSegmentInfo(obj)
            nPoints = obj.getNOfPoints;
            nSeg = length(find(strncmp(obj.yNames, 'Fluxv', 5)))/nPoints;
            idxOut = find(strncmp(obj.yNames, 'Fluxv', 5));
        end
        
        
        % Methods for getting the model matrices
        function m = getModelMatrix(obj, matrixName, varargin)
            if (nargin~=2 && nargin~= 4)
                error('Wrong number of input parameters');
            end 
            if nargin == 2
                m = getModelMatrix@plasmaCircuit2DModel(obj,matrixName);
            else
                if isnumeric(varargin{1}) % Rows
                    rowidx = varargin{1};
                elseif ischar(varargin{1})||iscell(varargin{1})
                    row = varargin{1};
                    if strcmp(row, ':')
                        rowidx = 1 : length(obj.getYNames);
                    elseif strcmp(row, 'end')
                        rowidx = length(obj.getYNames);
                    else
                        rowidx = getSignalIndex(obj,row,'y');
                    end
                end
                
                if isnumeric(varargin{2}) % Columns
                    colidx = varargin{2};
                elseif ischar(varargin{2})||iscell(varargin{2})
                    col = varargin{2};
                    if strcmp(col, ':')
                        colidx = 1 : length(obj.getXNames);
                    elseif strcmp(row, 'end')
                        colidx = length(obj.getXNames);
                    else
                        colidx = getSignalIndex(obj,col,'x');
                    end
                end
                 
                m = getModelMatrix@plasmaCircuit2DModel(obj,matrixName,rowidx,colidx);
            end
        end
        
        function m = getAMatrix(obj,varargin)
            matrix = 'A';
            if nargin == 1
                m = obj.getModelMatrix(matrix);
            elseif nargin == 3
                m = obj.getModelMatrix(matrix, varargin{1}, varargin{2});
            end
        end
                    
        function m = getBMatrix(obj,varargin)
            matrix = 'B';
            if nargin == 1
                m = obj.getModelMatrix(matrix);
            elseif nargin == 3
                m = obj.getModelMatrix(matrix, varargin{1}, varargin{2});
            end
        end
        
        function m = getCMatrix(obj,varargin)
            matrix = 'C';
            if nargin == 1
                m = obj.getModelMatrix(matrix);
            elseif nargin == 3
                m = obj.getModelMatrix(matrix, varargin{1}, varargin{2});
            end
        end
                                   
        function m = getDMatrix(obj,varargin)
            matrix = 'D';
            if nargin == 1
                m = obj.getModelMatrix(matrix);
            elseif nargin == 3
                m = obj.getModelMatrix(matrix, varargin{1}, varargin{2});
            end
        end
            
        function flag = isPlasmaless(obj)
           flag = not(obj.Kplasm); 
        end
            
        % Return a state space model for the specified IOs    
        % INPUTIDS specifies the identifiers used to select the model inputs
        %          (either indexes or names can be used)
        % OUTPUTIDS specifies the identifiers used to select the model outputs
        %          (either indexes or names can be used)
        % STABFLAG if equal to 1, the plasma is artificially
        %          stabilized by reversing the sign of the unstable eigenvalue
        % DISTFLAG If distFlag is set to 1, the model accepts also poloidal beta and
        %          internal inductance as inputs. These additional inputs are appended
        %          after the standard ones
        % CDIDS    specifies the identifiers of the current driven circuits
        function [ssModel, T1, stateIndex] = getStateSpace(obj,inputIDs,outputIDs,varargin)

             % Check on input IDs
             if (ischar(inputIDs))
                 inputIDs = {inputIDs};
                 indexIn = obj.getSignalIndex(inputIDs, 'u');
             elseif(iscell(inputIDs))
                 indexIn = obj.getSignalIndex(inputIDs, 'u');
             else
                 indexIn=inputIDs;
             end 
 
             % Check on output IDs
             if (ischar(outputIDs))
                 outputIDs = {outputIDs};
                 indexOut = obj.getSignalIndex(outputIDs, 'y');
             elseif (iscell(outputIDs))
                 indexOut = obj.getSignalIndex(outputIDs, 'y');
             else
                 indexOut = outputIDs;
             end

             % Default values of the optional outputs
             T1 = eye(size(obj.L,1));
             stateIndex = 1:size(obj.L,1);
             
             if nargin < 6
                % No optional parameters
                ssA = -inv(obj.L)*obj.R;
                ssB = (inv(obj.L))*obj.S;   
                ssB = ssB(:, indexIn);
                ssC = obj.C(indexOut, :);
                ssD = obj.D(indexOut, indexIn);
                stateNames  = obj.xNames;
                inputNames  = obj.uNames(indexIn);
                outputNames = obj.yNames(indexOut);
                if isrow(inputNames), inputNames = inputNames'; end
                

                % state indexes of the voltage driven circuits (needed for
                % the disturbances matrices generation)
                indexVD = obj.PFStatePosition;
                % state indexes of the eddy currents
                indexEddy = obj.EddyStatePosition;
                if ~obj.isPlasmaless
                    % add Ipl
                    if isempty(find(strcmp(inputIDs,'Vpl'),1))
                        indexEddy = [indexEddy, obj.IpStatePosition];
                    else
                        indexVDin = [indexVD, obj.nOfInputs];
                        indexVD = [indexVD, obj.IpStatePosition]; 
                    end
                end
                 
             % Current driven circuit state transformation
             elseif nargin == 6
                cdNames = varargin{3};
                vdNames = inputIDs;

                if(ischar(cdNames) == 1)
                    cdNames = {cdNames};
                    indexCD = obj.getSignalIndex(cdNames,'x');
                elseif (iscell(cdNames)==1)
                    indexCD = obj.getSignalIndex(cdNames,'x');                    
                else
                    indexCD = cdNames;
                end

                % state indexes of the voltage driven circuits
                indexVD = setdiff(obj.PFStatePosition, indexCD);
                indexVDin = indexVD;
                
                % state indexes of the eddy current
                indexEddy = obj.EddyStatePosition;
                
                if ~obj.isPlasmaless
                    % Add Vloop to the inputs
                    if not(isempty(find(strcmp(vdNames,'Vpl'))))
                        indexVDin = [indexVD, obj.nOfInputs];
                        indexVD = [indexVD, obj.IpStatePosition];
                    end

                    % Consider plasma as a passive circuit
                    if (isempty(find(strcmp(vdNames,'Vpl'))) == 1) && (isempty(find(strcmp(cdNames,'Ipl'))) == 1)
                        indexEddy = [indexEddy, obj.IpStatePosition];
                    end
                end
                % State transformation to make the model compute the CD
                % currents derivatives
                % 1 = VD
                % 2 = CD
                % 3 = eddy+IP
                L11 = obj.L(indexVD, indexVD);
                L13 = obj.L(indexVD, indexEddy);
                L31 = obj.L(indexEddy, indexVD);
                L33 = obj.L(indexEddy, indexEddy);
                L1 = [L11, L13; L31, L33]; 

                L12 = obj.L(indexVD, indexCD);
                L32 = obj.L(indexEddy, indexCD);
                L2 = [L12; L32]; 

                % Modified R matrices
                R1 = obj.R(indexVD, indexVD);
                R3 = obj.R(indexEddy, indexEddy);
                R = [R1, zeros(size(R1, 1), size(R3, 2)); zeros(size(R3, 1), size(R1, 2)), R3];

                % Base Transformation matrices
                T1 = inv(L1);
                T2 = -T1 * L2;
                stateIndex = [indexVD, indexEddy];

                % Modified model matrices
                ssA = -R * T1;

%                B1 = [obj.S(indexVD,indexVDin); zeros(length(indexEddy), length(indexVD))];

%                 indexPoloidalCircuits = [obj.PFStatePosition obj.IpStatePosition];
%                 indexB1 = find(ismember(indexPoloidalCircuits, indexVD)); % remove circuits that are not in the input vector from the B columns
%                 indexB1 = find(ismember(indexB1, indexIn));
%                 B1 = B1(:, indexB1);

   B1 = [eye(length(indexVD)); zeros(length(indexEddy), length(indexVD))];
  indexPoloidalCircuits = [obj.PFStatePosition obj.IpStatePosition];
  indexB1 = find(ismember(indexPoloidalCircuits, indexVD)); % remove circuits that are not in the input vector from the B columns
  indexB1 = find(ismember(indexPoloidalCircuits(indexB1), indexIn));
  B1 = B1(:, indexB1);

                B2 = -R * T2;
                ssB = [B1 B2];

                C1 = obj.C(indexOut, indexVD);
                C2 = obj.C(indexOut, indexCD);
                C3 = obj.C(indexOut, indexEddy);
                ssC = [C1 C3] * T1;

%                 D1 = [obj.D(indexOut(1 : (end-1)), indexIn); zeros(1, length(indexIn))];
                D1 = obj.D(indexOut(1 : end), indexIn);
                D2 = [C1 C3] * T2 + C2;
                ssD = [D1 D2];
                stateNames  = obj.xNames([indexVD indexEddy]);
                try
                  inputNames  = {vdNames{1:end} cdNames{1:end}}';
                catch
                  inputNames = {vdNames{1:end}; cdNames{1:end}};
                end
                outputNames = obj.yNames(indexOut);
             end
            
            % Numerically stabilize the unstable plasma
            if nargin >= 4
                if varargin{1} % stabFlag is the 4th parameter
                    ssA = obj.stabilizePlasma(ssA);
                end
            end
            
            if nargin >= 5
                distFlag = varargin{2};                
                if distFlag
                    if obj.isPlasmaless
                        error('Cannot apply disturbances to a plasmaless model!')
                    end
                    if nargin == 6 % CD coils are specified
                        % Substitute E with -LE!
                        ssB = [ssB -ssA*obj.LE([indexVD, indexEddy], :)];
                        ssD = [ssD -ssC*obj.LE([indexVD, indexEddy], :) + obj.F(indexOut, :)];
                    else
                        % Use x - Ew as a new state variable; Fdot is supposed
                        % to be the zero matrix                    
                        ssB = [ssB ssA*obj.E([indexVD, indexEddy], :)];
                        ssD = [ssD ssC*obj.E([indexVD, indexEddy], :)+obj.F(indexOut, :)];
                    end
                else
                    ssB = [ssB zeros(size(ssB, 1), 2)];
                    ssD = [ssD zeros(size(ssC, 1), 2)];
                end 
                inputNames{end+1} = 'betapol*Ipl';
                inputNames{end+1} = 'li*Ipl';
            end
            
            ssA(isnan(ssA)) = 0;
            ssB(isnan(ssB)) = 0;
            ssC(isnan(ssC)) = 0;
            ssD(isnan(ssD)) = 0;
            ssModel = ss(ssA,ssB,ssC,ssD, ...
                'StateName', stateNames, ...
                'InputName', inputNames, ...
                'OutputName', outputNames); 
        end
        
    end
end
