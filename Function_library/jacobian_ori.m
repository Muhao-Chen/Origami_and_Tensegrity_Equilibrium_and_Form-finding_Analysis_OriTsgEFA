function [phpn_e,phTpn,theta]=jacobian_ori(node_in_hinge,N,E_n_total)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% This function output the jacobian matrix of theta to nodal coordinate
%
% Inputs:
%   node_in_hinge: hinge related nodal coordinates
%	N: node matrix (3 x n array for n nodes)
%   E_n_total: transformation matrix of all the hinges-related nodes
%
% Outputs:
%	phpn_e: partial theta/ partial n (12 x n_h matrix)
%	phTpn:  equilibrium matrix of hinge rotation-- jacobian matrix of origami
%	theta:  hinge angle vector
%%
n_h=size(node_in_hinge,1);  %number of hinges

phpn_e=zeros(12,n_h);    % partial theta/ partial n
theta=zeros(n_h,1);      % angle of all hinges in rad
for i=1:n_h
    rij=N(:,node_in_hinge(i,1))-N(:,node_in_hinge(i,2));
    rkj=N(:,node_in_hinge(i,3))-N(:,node_in_hinge(i,2));
    rkl=N(:,node_in_hinge(i,3))-N(:,node_in_hinge(i,4));

    m_temp=cross(rij,rkj);
    n_temp=cross(rkj,rkl);

    phpx_i=norm(rkj)/norm(m_temp)^2*m_temp;
    phpx_l=-norm(rkj)/norm(n_temp)^2*n_temp;
    phpx_j=(rij'*rkj/norm(rkj)^2-1)*phpx_i-rkl'*rkj/norm(rkj)^2*phpx_l;
    phpx_k=(rkl'*rkj/norm(rkj)^2-1)*phpx_l-rij'*rkj/norm(rkj)^2*phpx_i;
    phpn_e(:,i)=[phpx_i;phpx_j;phpx_k;phpx_l];

    % calculate hinge angle
    if (m_temp'*rkl)
        eta=sign(m_temp'*rkl);
    else
        eta=1;
    end
    theta(i)=mod(eta*real(acos(m_temp'*n_temp/(norm(m_temp)*norm(n_temp)))),2*pi);
end

% jacobian matrix of the whole structure

Cell_phpn=mat2cell(phpn_e,12,ones(1,size(phpn_e,2)));          % transfer matrix H into a cell: Cell_H
phTpn=E_n_total*blkdiag(Cell_phpn{:});