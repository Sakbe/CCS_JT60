classdef plasmaCircuit2DModel < handle
    % Abstract class (interface) for the plasmaCircuit2DModel object
    %
    % This interface specifies the services that should be provided by a
    % generic 2D plasma/circuit linearized model (i.e. a model that
    % describes the behaviour of the plant (plasma, active coils and passive
    % structures) around a given equilibrium.
    %
    % A number of properties have been defined. This permits to implement
    % some methods in the abstract interface.
    %
    % It will be the constructor of each class implementation that will
    % populate these properties.
    
    
    % Properties - These are all procteded (can be accessed only by the
    % derived classes, but not by the user!)
    properties (Access = protected)
        % # of manipulable inputs in the model (i.e. voltages applied to the
        % active and passive coils and to the plasma
        nOfInputs = NaN; % Is NaN whan not defined
        
        % # of non manuipulable inputs in the model (i.e. disturbaces)
        % It is assumed that both the value and the derivative of each
        % disturbance affects the model behaviour
        nOfDisturbances = NaN; % Is NaN whan not defined
        
        % # of outputs in the model
        nOfOutputs = NaN; % Is NaN whan not defined
        
        % # of states of the linear model (i.e. model order)
        nOfStates = NaN;
        
        % plasmaless flag - set equal to 1 if the model is plasmaless and 0
        % otherwise
        plasmalessFlag = NaN; % Is NaN whan not defined
        
        % Model matrices
        % A matrix
        A = NaN; % Is NaN whan not defined
        
        % B matrix
        B = NaN; % Is NaN whan not defined
        
        % C matrix
        C = NaN; % Is NaN whan not defined
        
        % D matrix
        D = NaN; % Is NaN whan not defined
        
        % E matrix
        E = NaN; % Is NaN whan not defined
        
        % F matrix
        F = NaN; % Is NaN whan not defined
        
        % Fdot matrix
        Fdot = NaN; % Is NaN whan not defined
        
        % L matrix
        L = NaN; % Is NaN whan not defined
        
        % R matrix
        R = NaN; % Is NaN whan not defined
        
        % LE matrix
        LE = NaN; % Is NaN whan not defined
        
        % S matrix (circuit polarities)
        S = NaN; % Is NaN whan not defined
        
        % Equilibrium vectors
        % Equilibrium state
        xEq = NaN; % Is NaN whan not defined
        
        % Equilibrium input
        uEq = NaN; % Is NaN whan not defined
        
        % Equilibrium disturbance
        wEq = NaN; % Is NaN whan not defined
        
        % Equilibrium output
        yEq = NaN; % Is NaN whan not defined
        
        % Input data
        % Input names
        uNames = NaN; % Is NaN whan not defined
        
        % Input positions
        uPositions = NaN; % Is NaN whan not defined
        
        % State names
        xNames = NaN; % Is NaN whan not defined
        
        % State positions
        xPositions = NaN; % Is NaN whan not defined
        
        PFStatePosition = NaN;
        EddyStatePosition = NaN;
        IpStatePosition = NaN;
        
        % Output names
        yNames = NaN; % Is NaN whan not defined
        
        % Output positions
        yPositions = NaN; % Is NaN whan not defined
        
        % Disturbance names
        wNames = NaN; % Is NaN whan not defined
        
        % Disturbance positions
        wPositions = NaN; % Is NaN whan not defined
        
        % Vessel description
        % It is assumed to that the mechanical structure is specified through
        % a pdemesh struct
        vesselData = [];
        
        % Equilibrium flux map
        % It is assumed that the equilibrium plasma configuration is
        % specified via a FEM mesh
        equilFluxData = [];
        
    end
    
    % Abstract methods (to be implemented by the derived classes)
    methods (Abstract)
        ret = plotPlasmaShape(obj,varargin);   
    end
    
    % Generic methods
    methods
        
    % Get matrix method
    % OUTMATRIX = getModelMatrix(OBJ,MATRIXNAME,[ROWS,COLS])
    % Here we can find a way to specify either the row and column indexes
    % or the names! 
    function outMatrix = getModelMatrix(obj,matrixName,varargin)
        % Selects matrix
        if nargin >= 2
            % Selects row and column indexes
            try
                switch matrixName
                    case 'A'
                        outMatrix = obj.A;
                    case 'B'
                        outMatrix = obj.B;
                    case 'C'
                        outMatrix = obj.C;
                    case 'D'
                        outMatrix = obj.D;
                    case 'E'
                        outMatrix = obj.E;
                    case 'F'
                        outMatrix = obj.F;
                    case 'Fdot'
                        outMatrix = obj.Fdot;
                    case 'L'
                        outMatrix = obj.L;
                    case 'R'
                        outMatrix = obj.R;
                    case 'LE'
                        outMatrix = obj.LE;
                    case 'S'
                        outMatrix = obj.S;
                    otherwise
                        error('The selected matrix does not exist');
                end
                if (nargin >= 3)
                    outMatrix = outMatrix(varargin{1},varargin{2});
                end
            catch
                error('The selected indexes are out of range');
            end
        else
            error('Wrong number of input parameters: getModelMatrix(matrixName,[rowIdx,columnIdx]');
        end
    end
    
    % Set matrix method
    % setModelMatrix(OBJ,MATRIXNAME,VALUE,[ROWS,COLS])
    % Here we can find a way to specify either the row and column indexes
    % or the names! 
    function setModelMatrix(obj,matrixName,value,varargin)   
        if nargin >= 4
%             value = num2str(value);
%             r = int2str(varargin{1});
%             c = int2str(varargin{2});
            r = varargin{1};
            c = varargin{2};
            for i = 1:length(r)
                for j = 1:length(c)
                    eval(['obj.' matrixName '(' int2str(r(i)) ',' int2str(c(j)) ') = ' num2str(value(i,j)) ';']);
                end
            end
        elseif nargin == 3
            for i = 1 : eval(['length(obj.' matrixName ')'])
                for j = 1 : eval(['length(obj.' matrixName ')'])
                    r = int2str(i);
                    c = int2str(j);           
                    eval(['obj.' matrixName  '(' r ',' c ') = ' num2str(value(i,j)) ';']);
                end
            end
        else
            error('Not enough input arguments');
        end
    end
        
    % Get the A matrix
    % A = getAMatrix(OBJ,[ROWS,COLS])
    function m = getAMatrix(obj,varargin)
        if (nargin == 1)
            m = obj.getModelMatrix('A');
        elseif (nargin > 1)
            m = obj.getModelMatrix('A',varargin{1},varargin{2});
        else
            error('Wrong number of input parameters: getAMatrix(rowIdx,columnIdx])');
        end
    end
        
    % Get the B matrix
    % B = getBMatrix(OBJ,[ROWS,COLS])
    function m = getBMatrix(obj,varargin)
       if (nargin == 1)
           m = obj.getModelMatrix('B');
       elseif (nargin > 1) 
           m = obj.getModelMatrix('B',varargin{1},varargin{2});
       else
           error('Wrong number of input parameters: getBMatrix(rowIdx,columnIdx])');
       end
    end
    
    % Get the C matrix
    % C = getAMatrix(OBJ,[ROWS,COLS])
    function m = getCMatrix(obj,varargin)
        if (nargin == 1)
            m = obj.getModelMatrix('C');
        elseif (nargin > 1)
            m = obj.getModelMatrix('C',varargin{1},varargin{2});
        else
            error('Wrong number of input parameters: getCMatrix(rowIdx,columnIdx])');
        end
    end
    
    % Get the D matrix
    % D = getDMatrix(OBJ,[ROWS,COLS])
    function m = getDMatrix(obj,varargin)
        if (nargin == 1)
            m = obj.getModelMatrix('D');
        elseif (nargin > 1)
            m = obj.getModelMatrix('D',varargin{1},varargin{2});
        else
            error('Wrong number of input parameters: getDMatrix(rowIdx,columnIdx])');
        end
    end
    
    % Get the E matrix
    % E = getEMatrix(OBJ,[ROWS,COLS])
    function m = getEMatrix(obj,varargin)
        if (nargin == 1)
            m = obj.getModelMatrix('E');
        elseif (nargin > 1)
            m = obj.getModelMatrix('E',varargin{1},varargin{2});
        else
            error('Wrong number of input parameters: getEMatrix(rowIdx,columnIdx])');
        end
    end
    
    % Get the F matrix
    % F = getFMatrix(OBJ,[ROWS,COLS])
    function m = getFMatrix(obj,varargin)
        if (nargin == 1)
            m = obj.getModelMatrix('F');
        elseif (nargin > 1)
            m = obj.getModelMatrix('F',varargin{1},varargin{2});
        else
            error('Wrong number of input parameters: getFMatrix(rowIdx,columnIdx])');
        end
    end
    
    % Get the Fdot matrix
    % Fdot = getFdotMatrix(OBJ,[ROWS,COLS])
    function m = getFdotMatrix(obj,varargin)
        if (nargin == 1)
            m = obj.getModelMatrix('Fdot');
        elseif (nargin > 1)
            m = obj.getModelMatrix('Fdot',varargin{1},varargin{2});
        else
            error('Wrong number of input parameters: getFdotMatrix(rowIdx,columnIdx])');
        end
    end
    
    % Get the L matrix
    % L = getLMatrix(OBJ,[ROWS,COLS])
    function m = getLMatrix(obj,varargin)
        if (nargin == 1)
            m = obj.getModelMatrix('L');
        elseif (nargin > 1)
            m = obj.getModelMatrix('L',varargin{1},varargin{2});
        else
            error('Wrong number of input parameters: getLMatrix(rowIdx,columnIdx])');
        end
    end
    
    % Get the R matrix
    % R = getRMatrix(OBJ,[ROWS,COLS])
    function m = getRMatrix(obj,varargin)
        if (nargin == 1)
            m = obj.getModelMatrix('R');
        elseif (nargin > 1)
            m = obj.getModelMatrix('R',varargin{1},varargin{2});
        else
            error('Wrong number of input parameters: getRMatrix(rowIdx,columnIdx])');
        end
    end
    
    % Set the R matrix
    % setRMatrix(OBJ,varargin)
    function setRMatrix(obj,varargin)
        if nargin == 2
            obj.setModelMatrix('R', varargin{1})
        elseif nargin == 4
            obj.setModelMatrix('R', varargin{1}, varargin{2}, varargin{3});
        else
            error('Wrong number of parameters')            
        end
        
    end
    
    % Get the LE matrix
    % LE = getLEMatrix(OBJ,[ROWS,COLS])
    function m = getLEMatrix(obj,varargin)
        if (nargin == 1)
            m = obj.getModelMatrix('LE');
        elseif (nargin > 1)
            m = obj.getModelMatrix('LE',varargin{1},varargin{2});
        else
            error('Wrong number of input parameters: getLEMatrix(rowIdx,columnIdx])');
        end
    end
    
    % Get the S matrix
    % S = getSMatrix(OBJ,[ROWS,COLS])
    function m = getSMatrix(obj)
        if (nargin == 1)
            m = obj.getModelMatrix('S');
        elseif (nargin > 1)
            m = obj.getModelMatrix('S',varargin{1},varargin{2});
        else
            error('Wrong number of input parameters: getSMatrix(rowIdx,columnIdx])');
        end
    end
    
    % Get signal names
    % yNames = getYNames(obj)
    function n = getYNames(obj,idx)
        n = obj.yNames;
        if nargin == 2, n = n(idx); end
    end
    
    % uNames = getUNames(obj)
    function n = getUNames(obj,idx)
        n = obj.uNames;
        if nargin == 2, n = n(idx); end
    end
    
    % xNames = getXNames(obj)
    function n = getXNames(obj,idx)
        n = obj.xNames;
        if nargin == 2, n = n(idx); end
    end
    
    % wNames = getWNames(obj)
    function n = getWNames(obj,idx)
        n = obj.wNames;
        if nargin == 2, n = n(idx); end
    end
    
    
    
    % Get signal indexes
    % INDEX = getSignalIndex(OBJ,NAME,SIGNALTYPE)
    % If NAME is a cell array of strings INDEX contains the vector of
    % indexes.
    %
    % SIGNALTYPE can assume the following values
    % 'x' for state signal
    % 'u' for input signal
    % 'y' for output signal
    % 'w' for disturbance signal
    %
    % If a name is present more than once in the correspondent name vector
    % then only the first index is returned
    function index = getSignalIndex(obj,name,signalType)
    index = [];
    if nargin == 3
        switch signalType
            case 'u'
                if iscell(name)
                    for i = 1:length(name)
                        temp = obj.uPositions(find(strcmp(obj.uNames,name{i})==1));
                        if ~isempty(temp)
                            index(end+1) = temp(1);
                        end
                    end
                else
                    temp = obj.uPositions(find(strcmp(obj.uNames,name)==1));
                    index = temp(1);
                end
            case 'x'
                if iscell(name)
                    for i = 1:length(name)
                        temp = obj.xPositions(find(strcmp(obj.xNames,name{i})==1));
                        if ~isempty(temp)
                            index(end+1) = temp(1);
                        end
                    end
                else
                    temp = obj.xPositions(find(strcmp(obj.xNames,name)==1));
                    index = temp(1);
                end
            case 'y'
                if iscell(name)
                    for i = 1:length(name)
                        temp = obj.yPositions(find(strcmp(obj.yNames,name{i})==1));
                        if ~isempty(temp)
                            index(end+1) = temp(1);
                        end
                    end
                else
                    temp = obj.yPositions(find(strcmp(obj.yNames,name)==1));
                    index = temp(1);
                end
            case 'w'
                if iscell(name)
                    for i = 1:length(name)
                        temp = obj.wPositions(find(strcmp(obj.wNames,name{i})==1));
                        if ~isempty(temp)
                            index(end+1) = temp(1);
                        end
                    end
                else
                    temp = obj.wPositions(find(strcmp(obj.wNames,name)==1));
                    index = temp(1);
                end
            otherwise
                error('Wrong SIGNALTYPE');
        end
    else
        error('Wrong number of input parameters: getSignalIndex(obj,name,signalType)');
    end
    end
    
    % Get input indexes
    % INDEX = getInputIndex(OBJ,NAME)
    % If NAME is a cell array of strings INDEX contains the vector of
    % indexes.
    % If a name is present more than once in the input names vector
    % then only the first index is returned
    function index = getInputIndex(obj,name)
        index = obj.getSignalIndex(name,'u');
    end
    
    % Get state indexes
    % INDEX = getStateIndex(OBJ,NAME)
    % If NAME is a cell array of strings INDEX contains the vector of
    % indexes.
    % If a name is present more than once in the state names vector
    % then only the first index is returned
    function index = getStateIndex(obj,name)
        index = obj.getSignalIndex(name,'x');
    end
    
    % Get output indexes
    % INDEX = getOutputIndex(OBJ,NAME)
    % If NAME is a cell array of strings INDEX contains the vector of
    % indexes.
    % If a name is present more than once in the output names vector
    % then only the first index is returned
    function index = getOutputIndex(obj,name)
        index = obj.getSignalIndex(name,'y');
    end
    
    % Get disturbances indexes
    % INDEX = getDisturbanceIndex(OBJ,NAME)
    % If NAME is a cell array of strings INDEX contains the vector of
    % indexes.
    % If a name is present more than once in the disturbance names vector
    % then only the first index is returned
    function index = getDisturbanceIndex(obj,name)
        index = obj.getSignalIndex(name,'w');
    end
    
    
    % Return a state space model
    % SSMODEL = getStateSpace(OBJ,INPUTIDS,OUTPUTIDS,[STABFLAG],[DISTFLAG],[CDIDS])
    % Return a state space model of the plasma/circuit model
    %
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
    function ssModel = getStateSpace(obj,inputIDs,outputIDs,varargin)
        ssA = []; % tbd
        ssB = []; % tbd
        ssC = []; % tbd
        ssD = []; % tbd
        ssModel = ss(ssA,ssB,ssC,ssD);
    end
    
    % Return equilibrium outputs
    function yEquil = getYEq(obj,varargin)
        if nargin == 1
            yEquil = obj.yEq;
        elseif nargin == 2
            yIDs = varargin{1};
            index = obj.getOutputIndex(yIDs);
            yEquil = obj.yEq(index); 
        elseif nargin > 2
            error('Too many input arguments.')
        end
    end
    
    % Plot the poloidal section of the tokamak mechanical structure
    % AH = plotMechanicalStructure(OBJ,[COLOR])
    % COLOR is the color specifier
    function ah = plotMechanicalStructure(obj,varargin)
        if ~isempty(obj.vesselData)
            ah = pdemesh(obj.vesselData.p,obj.vesselData.e,[]);
            axis equal
            xlabel('R[m]', 'Interpreter', 'latex')
            ylabel('Z[m]', 'Interpreter', 'latex')
            if nargin == 2
                set(ah,'Color',varargin{1});
            end
        else
            error('Vessel data not available');
        end
    end
    
    % Plot the 2D plasma equilibrium shape from the equilibrium mesh
    % data
    % AH = plotEquilibriumShape(OBJ,[COLOR])
    % COLOR is the color specifier
    function ah = plotEquilibriumShape(obj,varargin)
        if ~isempty(obj.equilFluxData)
            ah = pdecont(obj.equilFluxData.p,obj.equilFluxData.t,obj.equilFluxData.nodeData,obj.equilFluxData.boundaryFlux*[1 1]);
            %pdecont(Input_struct.p,Input_struct.t,x_np(1:length(Input_struct.p)),y_np(strcmp(y_type,'psb_c'))*[1 1]);
            if nargin == 2
                set(ah,'Color',varargin{1});
            end
        else
            error('Equilibrium data not available');
        end
    end
    
end

methods (Static)
    % Reverse the sign of the unstable eigenvalue of the input dynamic
    % matrix
    function stabA = stabilizePlasma(inA)
        % Spectral decomposition of the input dynamic A matrix
        [VV, DD] = eig(inA);
        % Index of the unstable mode
        u_index = find(diag(DD) == max(max(real(DD)))); % if the eigenvalues are all negative, max(D) = 0 and u_index is empty
        if(isempty(u_index))
            error('There are no positive eigenvalues in the A matrix.');
        end
        % Reverse the unstable mode
        DD(u_index, u_index) = -DD(u_index, u_index);
        stabA = real(VV * DD / VV);
    end
    
    % Generate the initial state for a VDE (expressed in cm)
    function x0 = generateVDE(inA, indexZp, vde)
        [VV,DD] = eig(inA);
        idx_unstable = find(diag(DD>1e-5));
        if(idx_unstable)
            Cz = obj.C(indexZp,:) * VV; % could use the getModelMatrix function
            z0 = 1e-2; % vde expressed in cm
            xnewi0 = z0/Cz(idx_unstable);
            xxnewi0 = zeros(obj.nOfStates,1);
            xxnewi0(idx_unstable) = xnewi0;
            x0 = VV * xxnewi0 * vde;
        else
            error('There are no positive eigenvalues in the A matrix.')
        end
    end
end
end