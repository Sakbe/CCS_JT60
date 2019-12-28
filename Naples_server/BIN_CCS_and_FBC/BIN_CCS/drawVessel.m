function ret = drawVessel(p,e,axesHandle,color)
% ret = drawVessel(P,E,[AXESHANDLE],[COLOR])
% Draw a poloidal section of the tokamak vessel 
% Returns 0 if an error occurs and 1 otherwise
%
% the P and E inputs are the parameters that can be found in the equilbrium
% file under the Input_struct structure
%
% Optional parameters:
% - COLOR      (default = 'k')
% - AXESHANDLE (default = gca)
try
    ret = 0;
    % Check input parameters
    if nargin < 2
        error('You need to define at least P and E inputs')
        return            
    end
    if nargin < 4
        color = 'k';
        if nargin < 3
            axesHandle = gca;    
        end
    end
    h = pdemesh(p,e,[]);
    h.Color = color;
    % This axis limit will work only for JT-60SA
    axis('equal')
    axis([0, 6, -4.5, 4.5])
    ret = 1;
catch
    ret = 0;
end