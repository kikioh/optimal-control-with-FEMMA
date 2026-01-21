% Solve a linear system Kx=f with the constraints that x(ifix)=rho_init. 
% Use the adjoint method to compute gradient. Syntax:
%
%   [x,grad] = femsolve(control,K,Kgrad_L,Kgrad_R)
%
% Parameters:
%
%	control				- structure for optimal control 
%   K                   - square matrix, global stifness matrix
%   Kgrad_L             - gradient of the stifness matrix to control sequence, required if grad is a return
%   Kgrad_R (optional)  - right-hand gradient of the stifness matrix for piecewise-linear waveform
%
% Outputs:
%
%   x                   - solution of K * x = f
%
%   grad                - gradient of the fedelity function respect to control sequence
%
% mengjia.he@kit.edu, 2025.09.09

% function [t3,t4,x,t5,grad] = femsolveTest(control,K,Kgrad_L,Kgrad_R)
function [x,grad] = femsolve(control,K,Kgrad_L,Kgrad_R)
% numRep = 10;
% t3 = zeros(1,numRep);
% t4 = zeros(1,numRep);
% t5 = zeros(1,numRep);

% extract parameters
rho_init = control.rho_init{1};					% initial spin state
rho_targ = control.rho_targ{1};					% target spin state
dof = size(K,1);                                % number of total dofs
givenDof = numel(rho_init);                     % number of given dofs

% index of boundary constraints
ifix = 1:givenDof;								% index of boundary constraints
zero_index = 1:givenDof*4;						% set a large number for high-order elements

% update load vector
f = sparse(dof,1);  							% all-zero load vector
f = f - K(:,ifix) * rho_init;                   % move columns of KI*dI to the rhs
f(ifix) = rho_init;                             % set rhs for dof constraint equations

% update stiffness matrix
K(zero_index,ifix) = 0;                         % zero out columns
K(ifix,zero_index) = 0;                         % zero out rows
K(ifix,ifix) = speye(givenDof);                 % put ones on the diagonal

% solve spin trajectory on the nodes

% for m = 1:numRep
	% tic;
	dK = decomposition(K, 'banded'); 
	% t3(m) = toc;
% end
% t3 = median(t3);

 
% for m = 1:numRep
	% tic;
	x = dK \ f;
	% t4(m) = toc;
% end
% t4 = median(t4);

if nargout == 2
	switch control.integrator
	
		
		case 'rectangle'	% gradient for piecewise-constant waveform
		
		 
		% for p = 1:numRep
		
			% tic;
			
			Kgrad = Kgrad_L;
			nv = size(control.waveform,1)*size(control.waveform,2);		% number of control variables
			nc = size(control.waveform,1);								% number of control channels
			eleDof = size(Kgrad,1) / nv;								% dofs of element
			nodeDof = numel(rho_targ);									% dofs of node
			nNodes = dof / nodeDof;                     				% number of nodes
			ne = (nNodes-1)/(eleDof/nodeDof-1);         				% number of elements 
			dofs_matrix = dofIndex(nodeDof,eleDof,ne,nc,'strided');     % matrix of element index
			dofs_vector = dofs_matrix(:);

			% pick solutions correspond to elements
			eleSol = x(dofs_vector);                    

			% boundary conditions
			for m = 1:nc
				rows = (m-1)*eleDof*ne+1 : (m-1)*eleDof*ne + givenDof;
				cols = (m-1)*eleDof*ne + givenDof + 1: (m-1)*eleDof*ne + eleDof;
				Kgrad(rows, cols) = 0;
			
			end

			% calculate forces correspond to elements 
			eleForce = Kgrad * eleSol;              

			% sparse global force matrix
			rows = dofs_vector;                                 % row index [eleDof*nv×1]
			cols = repelem(1:nv, eleDof)';                      % column index [eleDof*nv×1]
			globalForce = sparse(rows, cols, eleForce, dof, nv);  % assembly global forces

			% objective with respect to target state
			gx = [zeros(1,dof-nodeDof),rho_targ'];

			% solve adjoint equation
			lambda = gx / dK;

			% gradient with respect to control variables
			grad = -lambda * globalForce;			
			
			% reshape gradient
			grad = reshape(reshape(grad, ne, []).', [], 1).';
			
			% t5(p) = toc;
		% end
		% t5 = median(t5);
		
		case 'trapezium' 												% gradient for piecewise-linear waveform
		
			nv = size(control.waveform,1)*size(control.waveform,2);		% number of control variables
			nc = size(control.waveform,1);								% number of control channels
			rhoDof = numel(rho_targ);
			nodeDof = rhoDof * 2;
			nNodes = dof / nodeDof;                     				% number of nodes
			eleDof = size(Kgrad_L,1) / (nv-nc);                         % element dofs
			ne = (nNodes-1)/(eleDof/nodeDof-1);         				% number of elements 
			dofs_matrix = dofIndex(nodeDof,eleDof,ne,nc,'strided');     % matrix of element index
			dofs_vector = dofs_matrix(:);
			
			% pick solutions correspond to elements
			eleSol = x(dofs_vector); 
			
			% boundary conditions
			for m = 1:nc
				rows = (m-1)*eleDof*ne + 1 : (m-1)*eleDof*ne + givenDof;
				cols = (m-1)*eleDof*ne + givenDof + 1: (m-1)*eleDof*ne + eleDof;
				Kgrad_R(rows, cols) = 0;
			
			end
			
			% calculate forces correspond to elements 
			eleForce_L = Kgrad_L * eleSol;   
			eleForce_R = Kgrad_R * eleSol;
			
			% sparse global force matrix 
			cols_L = cols_index(ne,nc,eleDof)+1;                     			% column index for left-side force [eleDof*(nc*ne)×1]      
			cols_R = cols_index(ne,nc,eleDof);                       			% column index for right-side force [eleDof*(nc*ne)×1] 			
			force_L = sparse(dofs_vector, cols_L, eleForce_L(:), dof, nv);  	% assembly left-side global forces
			force_R = sparse(dofs_vector, cols_R, eleForce_R(:), dof, nv);  	% assembly right-side global forces
			globalForce = force_L + force_R;
			
			% objective with respect to target state
			gx = [zeros(1,dof-nodeDof), rho_targ', zeros(1,rhoDof)];

			% solve adjoint equation
			lambda = gx / dK;

			% gradient with respect to control variables
			grad = -lambda * globalForce;
			
			% reshape gradient
			grad = reshape(reshape(grad, ne+1, []).', [], 1).';
	
	end
	

end

end
