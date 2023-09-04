function M=tenseg_mass_matrix(mass,C,lumped)
% This function is written by Ma, S., Chen, M. and Skelton, R., 2022. TsgFEM: Tensegrity finite element method. Journal of Open Source Software, 7(75), p.3390. 
% The code license is Mozilla Public License, v. 2.0.
% The source code is here:
% https://github.com/Muhao-Chen/Tensegrity_Finite_Element_Method_TsgFEM/tree/main/Function_library/tenseg_mass_matrix.m
%
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
%
% This function generate mass matrix of tensegrity structure.
% 
% Inputs:
%   mass: use lumped matrix 1-yes,0-no
%   C: connectivity matrix
%   lumped: 1 for lumped matrix; 0 for consistent matrix
% Outputs:
%	M: mass matrix 
%%
switch lumped
    case 1  % lumped matrix
M=0.5*kron(diag(diag(abs(C)'*diag(mass)*abs(C))),eye(3));     
    case 0  % consistent matrix
M=1/6*kron((abs(C)'*diag(mass)*abs(C)+diag(diag(abs(C)'*diag(mass)*abs(C)))),eye(3));
end
end

