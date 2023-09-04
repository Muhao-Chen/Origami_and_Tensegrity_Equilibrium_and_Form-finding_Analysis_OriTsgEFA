function [Kt_aa,Kg_aa,Ke_aa,K_mode,k]=tenseg_stiff_CTS2(Ia,C,q,A_2ac,E_c,A_c,l0_c)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% This function calculates the tangent stiffness matrix information of
% CTS(clustered tensegrity structures). 
%
% Inputs:
%   Ia: transfer matrix of free coordinates
%   C: tensegrity members connectivity matrix
%   q: axial members force density vector 
%   A_2ac: equilibrium matrix with constraints and group, force as variable
%   E_c: cluster members' cross sectional vector
%   A_c: cluster members' Young's modulus vector
%   l0_c: cluster members' length vector
% Outputs:
%	Kt_aa: truss elements tangent matrix
%	Kg_aa: truss elements geometry matrix
%	Ke_aa: truss elements material matrix
%	K_mode: eigenvector of tangent stiffness matrix
%   k: eigenvalue of tangent stiffness matrix

Kg_aa=Ia'*kron(C'*diag(q)*C,eye(3))*Ia;
Ke_aa=A_2ac*diag(E_c.*A_c./l0_c)*A_2ac';
Kt_aa=Kg_aa+(Ke_aa+Ke_aa')/2;       % this is to 
[K_mode,D1] = eig(Kt_aa);         % eigenvalue of tangent stiffness matrix
k=diag(D1);   

end

