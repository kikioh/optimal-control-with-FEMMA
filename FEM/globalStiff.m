% assembly the global stiffness matrix and its gradient. Syntax:
%
%       [K,Kgrad,Kgrad_R] = globalStiff(fem,control)
%
% Parameters:
%
%		fem
%       	fem.T                  	- total time
%			fem.nEmems				- number of elements		
%       	fem.nNodes           	- number of nodes	
%       	fem.dofNode          	- dofs of node
%			fem.KIndex (optional)	- Index matrix for assembly K, pre-computed for 'LvN' equation
%
%		control
%       	control.offsets         - 1*numSpin cell array, resonance offset in Hz
%       	control.pwr_level      	- 1*numSpin cell array, B1 amplitude for all control channels in rad/s
%       	control.waveform      	- numChn by numP matrix, represent phase or xy waveform
%
% Output:
%
%       K                   	- global stiffness matrix
%       Kgrad               	- global stiffness matrix gradient wrt controls
%       Kgrad_R (optional)      - right-side gradient of global stiffness matrix wrt controls if Hemerite functions is used
%
% mengjia.he@kit.edu, 2025.09.04

function [K,Kgrad,Kgrad_R] = globalStiff(fem,control)

% uniform parameters
dt = fem.T/fem.nElems;              % slice duratoin
nElems = fem.nElems;                % number of elements
nNodes = fem.nNodes;                % number of nodes
dofNode = fem.dofNode;              % dofs of each node

% switch  fem.equation
    
    % case 'Helmholtz'

        % % Helmholtz filter radius
        % R = fem.R;

        % % element stiffness matrix, independent on coordinates
        % S_element = [dt/3 + R^2/(2*dt), dt/6 - R^2/(2*dt);
                     % dt/6 - R^2/(2*dt), dt/3 + R^2/(2*dt)];
        % % combine element stiffness matrix
        % S = repmat(S_element, 1, 1, nElems);

        % % assemble global stiffness matrix
        % [rows,cols] = stiffIndex(nElems,dofNode);
        % K = sparse(rows(:), cols(:), S(:), nNodes*dofNode, nNodes*dofNode);
		
		% if nargout>1, Kgrad = []; end
		% if nargout>2, Kgrad_R = []; end

    % case 'LvN'

% preallocate global stiffness matrix
if nargout == 1
	
	% element stiffness matrix
	switch fem.shapeOrder
	
		case 'linear', S = eleMatrix_linear(control);
		
		case 'quadratic', S = eleMatrix_quad(control);
		
		case 'hermite', S = eleMatrix_hermite(control);
		
	end

elseif nargout >= 2
	
	% element stiffness matrix and gradient
	switch fem.shapeOrder
	
		case 'linear', [S,Kgrad] = eleMatrix_linear(control);
		
		case 'quadratic', [S,Kgrad] = eleMatrix_quad(control); 
		
		case 'hermite', [S,Kgrad,Kgrad_R] = eleMatrix_hermite(control); 
		
	end

end

% assemble global stiffness matrix
rows = fem.KIndex.rows; cols = fem.KIndex.cols; 
K = sparse(rows(:), cols(:), S(:), dofNode*nNodes, dofNode*nNodes);

% end

end

