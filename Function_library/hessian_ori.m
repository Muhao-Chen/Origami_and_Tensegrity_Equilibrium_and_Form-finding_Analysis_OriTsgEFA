function [ph2pn2_e,ph2pn2]=hessian_ori(node_in_hinge,N,E_n)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% This function output the hessian matrix of theta to nodal
% coordinate,  and the jocobian matrix(modified from Yangli's code)
%
% Inputs:
%   node_in_hinge: hinge related nodal coordinates
%	N: node matrix (3 x n array for n nodes)
%   E_n: transformation matrix of all the hinges-related nodes (cell form)
%
% Outputs:
%	ph2pn2_e: Hessian matrix for the hinges (every cell contains a 12 x 12 matrix)
%	ph2pn2:  Second-order partial derivative of the hinge angle for whole
%	structure (cell form)
%%
n_h=size(node_in_hinge,1);  %number of hinges

ph2pn2=cell(1,n_h);          %ph2px2 is the hessian matrix of theta to nodal coordinate
ph2pn2_e=cell(1,n_h);
for i=1:n_h
    rij=N(:,node_in_hinge(i,1))-N(:,node_in_hinge(i,2));
    rkj=N(:,node_in_hinge(i,3))-N(:,node_in_hinge(i,2));
    rkl=N(:,node_in_hinge(i,3))-N(:,node_in_hinge(i,4));

    rmj=cross(rij,rkj);
    rnk=cross(rkj,rkl);
    % ph2px2_eii= -norm(r_kj)/norm(m_temp)^4*((m_temp*cross(r_kj,m_temp)')+(m_temp*cross(r_kj,m_temp)')');


    %%
    di = norm(rkj)/(rmj'*rmj)*rmj;
    dl = -norm(rkj)/(rnk'*rnk)*rnk;
    dj = (rij'*rkj/(rkj'*rkj)-1)*di-rkl'*rkj/(rkj'*rkj)*dl;
    dk = -rij'*rkj/(rkj'*rkj)*di+(rkl'*rkj/(rkj'*rkj)-1)*dl;

    dii = -norm(rkj)/(rmj'*rmj)^2*((rmj*cross(rkj,rmj)')+(rmj*cross(rkj,rmj)')');
    % dii2=-norm(rkj)/(rmj'*rmj)^2*(skew(rkj)*(rmj'*rmj)-2*rmj*rmj'*skew(rkj));

    dtempij = -norm(rkj)/(rmj'*rmj)^2*(rmj*(cross(rij-rkj,rmj))'+(cross(rij-rkj,rmj))*rmj');
    dij = -rmj*rkj'/(rmj'*rmj*norm(rkj))+dtempij;

    dtempik = norm(rkj)/(rmj'*rmj)^2*(rmj*(cross(rij,rmj))'+(cross(rij,rmj))*rmj');
    dik = rmj*rkj'/(rmj'*rmj*norm(rkj))+dtempik;

    dli = zeros(3);

    dll = norm(rkj)/(rnk'*rnk)^2*(rnk*cross(rkj,rnk)'+(rnk*cross(rkj,rnk)')');

    dtemplk = norm(rkj)/(rnk'*rnk)^2*(rnk*(cross(rkl-rkj,rnk))'+(cross(rkl-rkj,rnk))*rnk');
    dlk = -rnk*rkj'/(rnk'*rnk*norm(rkj))+dtemplk;

    dtemplj = norm(rkj)/(rnk'*rnk)^2*(rnk*(cross(rnk,rkl))'+(rnk*(cross(rnk,rkl))')');
    dlj = rnk*rkj'/(rnk'*rnk*norm(rkj))+dtemplj;

    dT1jj = 1/(rkj'*rkj)*((-1+2*rij'*rkj/(rkj'*rkj))*rkj-rij);
    dT2jj = 1/(rkj'*rkj)*(2*rkl'*rkj/(rkj'*rkj)*rkj-rkl);
    djj = di*dT1jj'+(rij'*rkj/(rkj'*rkj)-1)*dij-(dl*dT2jj'+rkl'*rkj/(rkj'*rkj)*dlj);

    dT1jk = 1/(rkj'*rkj)*(-2*rij'*rkj/(rkj'*rkj)*rkj+rij);
    dT2jk = 1/(rkj'*rkj)*((1-2*rkl'*rkj/(rkj'*rkj))*rkj+rkl);
    djk = di*dT1jk'+(rij'*rkj/(rkj'*rkj)-1)*dik-(dl*dT2jk'+rkl'*rkj/(rkj'*rkj)*dlk);

    dT1kk = dT2jk;
    dT2kk = dT1jk;
    dkk = dl*dT1kk'+(rkl'*rkj/(rkj'*rkj)-1)*dlk-(di*dT2kk'+rij'*rkj/(rkj'*rkj)*dik);

    ph2pn2_e{i} = [ dii , dij , dik, dli';
        dij', djj , djk, dlj';
        dik' , djk' , dkk , dlk' ;
        dli , dlj, dlk, dll];

    ph2pn2{i}=  E_n{i}*ph2pn2_e{i}*E_n{i}';
end