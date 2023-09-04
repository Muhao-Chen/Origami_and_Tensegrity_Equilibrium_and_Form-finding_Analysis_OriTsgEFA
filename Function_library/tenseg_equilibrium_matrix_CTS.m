function [A_1,A_1c,A_1a,A_1ac,A_2,A_2c,A_2a,A_2ac,l,l_c]=tenseg_equilibrium_matrix_CTS(N,C,S,Ia)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% This function gives the equilibrium matrix of tensegrity structures
%
% Inputs:
%   N: node matrix (3 x n array for n nodes)
%	C: members' length
%   S: clustering matrix
%   Ia: Tranform matrix to get free nodal coordinate: na=Ia'*n
%
% Outputs:
%	A_1: equilibrium matrix with no constraints, force density as variable
%	A_1c: Equilirium matrix with group, force density as variable.
%   A_1a: equilibrium matrix with constraints, force density as variable
%   A_1ac: equilibrium matrix with constraints and group, force density as variable
%	A_2: equilibrium matrix with no constraints, force as variable
%	A_2c: Equilirium matrix with group, force as variable.
%   A_2a: equilibrium matrix with constraints, force as variable
%   A_2ac: equilibrium matrix with constraints and group, force as variable
%   l: members' length vector
%   l_c: cluster members' length vector 
%%

% element length
H=N*C';                     % element's direction matrix
l=sqrt(diag(H'*H));         % elements' length
l_c=pinv(S')*l;            % elements' length in group

Cell_H=mat2cell(H,3,ones(1,size(H,2)));          % transfer matrix H into a cell: Cell_H

A_1=kron(C',eye(3))*blkdiag(Cell_H{:});     % equilibrium matrix
A_1a=Ia'*A_1;
A_1c=A_1*diag(l.^-1)*S'*diag(l_c); 
A_1ac=Ia'*A_1c;

A_2=A_1*diag(l.^-1);                           % equilibrium matrix
A_2a=Ia'*A_2;                                   % equilibrium matrix in group constraints
A_2c=A_1c*diag(l_c.^-1);
A_2ac=Ia'*A_2c;


end

