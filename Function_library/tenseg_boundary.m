function [Ia,Ib,a,b]=tenseg_boundary(pinned_X,pinned_Y,pinned_Z,nn)
% This function is written by Ma, S., Chen, M. and Skelton, R., 2022. TsgFEM: Tensegrity finite element method. Journal of Open Source Software, 7(75), p.3390. 
% The code license is Mozilla Public License, v. 2.0.
% The source code is here:
% https://github.com/Muhao-Chen/Tensegrity_Finite_Element_Method_TsgFEM/tree/main/Function_library/tenseg_boundary.m
%
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% This function calculate the free and pinned number of nodal coordinate,
% as well as the transfer matrix of free and pinned coordinates
% Inputs:
%   pinned_X: number of nodes pinned in X direction  
%	pinned_Y: number of nodes pinned in Y direction 
%   pinned_Z: number of nodes pinned in Z direction 
%
% Outputs:
%	Ia: transfer matrix of free coordinates
%	Ib: transfer matrix of pinned coordinates
%   a: number of free coordinates
%   b: number of pinned coordinates

%%
I=eye(3*nn);
b=sort([3*pinned_X-2;3*pinned_Y-1;3*pinned_Z]);   %index of pinned nodes
a=setdiff((1:3*nn)',b);  %index of free node direction
Ia=I(:,a);  %free node index
Ib=I(:,b);  %pinned nod index
end

