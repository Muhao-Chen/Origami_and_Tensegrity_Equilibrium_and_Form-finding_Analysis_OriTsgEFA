function data_out=static_solver_ori_2(data)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% solve nonlinear equilibrium equations using modified Newton method
% converge to stable equilibrium, considering substep, for origami( including
% CTS, TTS), use MGDCM in the iteration
%
% Inputs:
%	data_in: data structure describing simulation task
% Outputs:
%	data_out: data structure containing simulation results


global E A l0 Ia Ib C S w ne Xb Xa dXa f_int l_int theta_0 k_h E_n_total node_in_hinge
% minimize total energy? (1: use, 0: not use) it's time consuming
use_energy=1;

%% input data
C=data.C;
ne=data.ne;
nn=data.nn;
Ia=data.Ia;
Ib=data.Ib;
S=data.S;
index_b=data.index_b;
index_s=data.index_s;
substep=data.substep;
E0=data.E;
consti_data=data.consti_data;
material=data.material;
A=data.A;
theta_0_t=data.theta_0_t;   % initial angle
k_h=data.k_h;               % stiffness of hinge
E_n=data.E_n;               % transfer matrix from matrix to structure
E_n_total=cell2mat(E_n);        % transformation matrix of the whole structure
node_in_hinge=data.node_in_hinge;       % node in triangle element in hinge

if  isfield(data,'subsubstep')
    subsubstep=data.subsubstep;
else
    subsubstep=30;          %default ssubstep
end
w_t=data.w_t;
dXb_t=data.dnb_t;
l0_t=data.l0_t;
X0=data.N(:);
data_out=data;     %initialize output data
data_out.E_out=E0*ones(1,substep);


%% calculate equilibrium
X=X0;               %initialize configuration
Xb0=Ib'*X;          %pinned node
E=E0;
Xb=Xb0;               %forced node displacement
l0=l0_t(:,1);         %forced enlongation of string
theta_0=theta_0_t(:,1);     % initial angle
% lamda=linspace(0,1,substep);    %coefficient for substep
num_slack=ne*zeros(substep+1,1);    %num of string slack
Xa0=Ia'*X;
Xa=Xa0;
cont=2;
u=1e-2;
tol = 1e-6; MaxIter = 150; 
U=zeros(3*nn,1);

 if strcmpi(data.LoadType, 'Force')             % load type: Force
    MaxIcr = data.MaxIcr;                   
    b_lambda = data.InitialLoadFactor;          
    Uhis = zeros(3*nn,MaxIcr);        
    FreeDofs = find(sum(Ia,2));
    lmd = 0; icrm = 0; MUL = [U,U];
    Fhis = zeros(MaxIcr,1);
    data_out.t_out=zeros(ne,MaxIcr);        %output member force
    data_out.theta_out=zeros(numel(theta_0),MaxIcr);      %angle of hinge
    data_out.M_out=zeros(numel(theta_0),MaxIcr);      %moment of hinge
    data_out.l_out=zeros(ne,MaxIcr);                % member length
    data_out.n_out=zeros(3*nn,MaxIcr);                % nodal coordinate
    data_out.Kt_aa_out=cell(1,MaxIcr);       %tangent stiffness of truss
    data_out.K_t_oa_out=cell(1,MaxIcr);       %tangent stiffness of whole struct.
   
    F=w_t(:,end);
    while icrm<MaxIcr && ~data.StopCriterion(U)&& lmd<=1 
        icrm = icrm+1;
        iter = 0; err = 1;
%         fprintf('icrm = %d, lambda = %6.4f\n',icrm,lmd);

        while err>tol && iter<MaxIter
            iter = iter+1;
%% equilibrium & tangent stiffness matrix

            X=[Ia';Ib']\[Xa;Xb];
            l=sqrt(sum((reshape(X,3,[])*C').^2))'; %bar length
            l_c=S*l;

            % member force (of truss)
            %         q=E.*A.*(1./l0-1./l);      %force density
            strain=(l_c-l0)./l0;        %strain of member
            [E,sigma]=stress_strain(consti_data,index_b,index_s,strain,material);
            t_c=sigma.*A;         %member force
            t=S'*t_c;
            q_c=t_c./l_c;
            q=t./l;      %reculate force density
            q_bar=diag(q);

            % angle and jacobian (of hinge)
            [phpn_e,phTpn,theta]=jacobian_ori(node_in_hinge,reshape(X,3,[]),E_n_total);       % jacobian matrix

            % moment of hinge
            M=diag(k_h)*(theta-theta_0);

            % equilibrium matrix (of truss)
            N=reshape(X,3,[]);
            H=N*C';
            Cell_H=mat2cell(H,3,ones(1,size(H,2)));          % transfer matrix H into a cell: Cell_H
            A_2a=Ia'*kron(C',eye(3))*blkdiag(Cell_H{:})*diag(l.^-1);     % equilibrium matrix
            A_2ac=A_2a*S';

            % equilibrium matrix (or truss + hinge)
            A_o_a=[A_2ac,Ia'*phTpn];
            % unbalanced force
            Ki=kron(C'*q_bar*C,eye(3));       %stiffness matrix of truss
            IF=Ki*X+phTpn*M;                   %unbalanced force
%             Fp_a=Ia'*Fp;                   %see the norm of unbalanced force
            % tangent stiffness matrix (of truss)
            Kg_aa=Ia'*Ki*Ia-A_2a*q_bar*A_2a';
            Ke_aa=A_2ac*diag(E.*A./l0)*A_2ac';
            Kt_aa=Kg_aa+(Ke_aa+Ke_aa')/2;

            % jacobian matrix (of hinge)
            [ph2pn2_e,ph2pn2]=hessian_ori(node_in_hinge,N,E_n);         % calculate hessian matrix
            G=cell2mat(ph2pn2);

            % tangenet stiffness (truss +hinge)
            K_t_oa=Kt_aa+Ia'*(phTpn*diag(k_h)*phTpn'+G*kron(M,eye(3*nn)))*Ia;


%             [IF,K] = GlobalK_fast_ver(U,Node,truss,angles);
            R = lmd*F-IF;   MRS = [F,R];
            MUL(FreeDofs,:) = K_t_oa\(Ia'*MRS);
            dUp = MUL(:,1); dUr = MUL(:,2);
            if iter==1, dUr = 0*dUr; end
            dlmd=nlsmgd(icrm,iter,dUp,dUr,b_lambda);
            dUt = dlmd*dUp+dUr;
            U = U+dUt;
            Xa=Xa0+Ia'*U;
            err = norm(Ia'*dUt);
            lmd = lmd+dlmd;
%             fprintf('    iter = %d, err = %6.4f, dlambda = %6.4f\n',iter,err,dlmd);
            if err > 1e8
                disp('Divergence!'); 
                break; 
            end
        end

        if iter>15
            b_lambda = b_lambda/2;
%             disp('Reduce constraint radius...')
            icrm = icrm-1;
            U = Uhis(:,max(icrm,1));  % restore displacement
            lmd = Fhis(max(icrm,1));   % restore load
            Xa=Xa0+Ia'*U;
        elseif iter<3
%             disp('Increase constraint radius...')
            b_lambda = b_lambda*1.5;
            Uhis(:,icrm) = U;
            Fhis(icrm) = lmd;
            data_out.t_out(:,icrm)=t_c;      %member force
            data_out.theta_out(:,icrm)=theta;      %angle of hinge
            data_out.M_out(:,icrm)=M;      %moment of hinge
            data_out.l_out(:,icrm)=l_c;      % length of bars
            data_out.n_out(:,icrm)=X;
            data_out.Kt_aa_out{:,icrm}=Kt_aa;       %tangent stiffness of truss
            data_out.K_t_oa_out{:,icrm}=K_t_oa;       %tangent stiffness of whole struct.
   
        else
            Uhis(:,icrm) = U;
            Fhis(icrm) = lmd;
            data_out.t_out(:,icrm)=t_c;      %member force
            data_out.theta_out(:,icrm)=theta;      %angle of hinge
            data_out.M_out(:,icrm)=M;      %moment of hinge
            data_out.l_out(:,icrm)=l_c;      % length of bars
            data_out.n_out(:,icrm)=X;
            data_out.Kt_aa_out{:,icrm}=Kt_aa;       %tangent stiffness of truss
            data_out.K_t_oa_out{:,icrm}=K_t_oa;       %tangent stiffness of whole struct.
   
        end
    end
elseif strcmpi(data.LoadType, 'Substep')
    substep=data.substep;
    MaxIcr =substep;                   
    b_lambda = data.InitialLoadFactor;          
    Uhis = zeros(3*nn,MaxIcr);        
    FreeDofs = find(sum(Ia,2));
    lmd = 0; icrm = 0; MUL = [U,U];
    Fhis = zeros(MaxIcr,1);
    data_out.t_out=zeros(ne,MaxIcr);        %output member force
    data_out.theta_out=zeros(numel(theta_0),MaxIcr);      %angle of hinge
    data_out.M_out=zeros(numel(theta_0),MaxIcr);      %moment of hinge
    data_out.l_out=zeros(ne,MaxIcr);                % member length
    data_out.n_out=zeros(3*nn,MaxIcr);                % nodal coordinate
    data_out.Kt_aa_out=cell(1,MaxIcr);       %tangent stiffness of truss
    data_out.K_t_oa_out=cell(1,MaxIcr);       %tangent stiffness of whole struct.
   
    F=w_t(:,end);
    while icrm<MaxIcr && ~data.StopCriterion(U)&& lmd<=1 
        icrm = icrm+1;
        iter = 0; err = 1;
        lmd=icrm/MaxIcr;
        Fhis(icrm) = lmd;
        w=w_t(:,icrm);               %external force
        Xb=Xb0+dXb_t(:,icrm);         %forced node displacement
        l0=l0_t(:,icrm);         %forced enlongation of string
        theta_0=theta_0_t(:,icrm);     % initial angle

%         fprintf('icrm = %d, lambda = %6.4f\n',icrm,lmd);

        while err>tol && iter<MaxIter
            iter = iter+1;
%% equilibrium & tangent stiffness matrix

            X=[Ia';Ib']\[Xa;Xb];
            l=sqrt(sum((reshape(X,3,[])*C').^2))'; %bar length
            l_c=S*l;

            % member force (of truss)
            %         q=E.*A.*(1./l0-1./l);      %force density
            strain=(l_c-l0)./l0;        %strain of member
            [E,sigma]=stress_strain(consti_data,index_b,index_s,strain,material);
            t_c=sigma.*A;         %member force
            t=S'*t_c;
            q_c=t_c./l_c;
            q=t./l;      %reculate force density
            q_bar=diag(q);

            % angle and jacobian (of hinge)
            [phpn_e,phTpn,theta]=jacobian_ori(node_in_hinge,reshape(X,3,[]),E_n_total);       % jacobian matrix

            % moment of hinge
            M=diag(k_h)*(theta-theta_0);

            % equilibrium matrix (of truss)
            N=reshape(X,3,[]);
            H=N*C';
            Cell_H=mat2cell(H,3,ones(1,size(H,2)));          % transfer matrix H into a cell: Cell_H
            A_2a=Ia'*kron(C',eye(3))*blkdiag(Cell_H{:})*diag(l.^-1);     % equilibrium matrix
            A_2ac=A_2a*S';

            % equilibrium matrix (or truss + hinge)
            A_o_a=[A_2ac,Ia'*phTpn];
            % unbalanced force
            Ki=kron(C'*q_bar*C,eye(3));       %stiffness matrix of truss
            IF=w-(Ki*X+phTpn*M);                   %unbalanced force
            IFa=Ia'*(w-(Ki*X+phTpn*M));            %unbalanced force
%             Fp_a=Ia'*Fp;                   %see the norm of unbalanced force
            % tangent stiffness matrix (of truss)
            Kg_aa=Ia'*Ki*Ia-A_2a*q_bar*A_2a';
            Ke_aa=A_2ac*diag(E.*A./l0)*A_2ac';
            Kt_aa=Kg_aa+(Ke_aa+Ke_aa')/2;

            % jacobian matrix (of hinge)
            [ph2pn2_e,ph2pn2]=hessian_ori(node_in_hinge,N,E_n);         % calculate hessian matrix
            G=cell2mat(ph2pn2);

            % tangenet stiffness (truss +hinge)
            K_t_oa=Kt_aa+Ia'*(phTpn*diag(k_h)*phTpn'+G*kron(M,eye(3*nn)))*Ia;


%             [IF,K] = GlobalK_fast_ver(U,Node,truss,angles);
            R = lmd*F-IF;   MRS = [F,R];
            dXa = K_t_oa\IFa;
            err = norm(dXa);
            Xa=Xa+dXa;
%             fprintf(' iter = %d, err = %6.4f\n',iter,err);
            if err > 1e8, disp('Divergence!'); break; end
        end


    data_out.n_out(:,icrm)=X;
    data_out.l_out(:,icrm)=l;
    %     data_out.q_out(:,k)=q;
    %     data_out.E_out(:,k)=E;
    data_out.t_out(:,icrm)=t;      %member force
    data_out.M_out(:,icrm)=M;      % moment of hinge
    data_out.theta_out(:,icrm)=theta;      %member force
    data_out.Kt_aa_out{:,icrm}=Kt_aa;       %tangent stiffness of truss
    data_out.K_t_oa_out{:,icrm}=K_t_oa;       %tangent stiffness of whole struct.
    % data_out.V{k}=energy_cal(data_out);
%     data_out.Fpn_out(icrm)=norm(IFa);
    end

 end
%% delete excessive part data
icrm = icrm+1;
Uhis(:,icrm:end) = [];
Fhis(icrm:end,:) = [];
data_out.t_out(:,icrm:end)=[];      %member force
data_out.theta_out(:,icrm:end)=[];      %angle of hinge
data_out.M_out(:,icrm:end)=[];      %moment of hinge
data_out.l_out(:,icrm:end)=[];     % length of bars
data_out.n_out(:,icrm:end)=[]; 
data_out.Kt_aa_out(:,icrm:end)=[];       %tangent stiffness of truss
data_out.K_t_oa_out(:,icrm:end)=[];       %tangent stiffness of whole struct.
    
% data_out.n_out=kron(X,ones(1,icrm-1))+Uhis;
data_out.Fhis=Fhis;

end


%--------------------------------------------------------------------------
function dl=nlsmgd(step,ite,dup,dur,cmp)
% Modified Generalized Displacement Control Method
global dupp1 sinal dupc1 numgsp
if ite==1
    if step==1
        sinal=sign(dot(dup,dup));
        dl=cmp;
        numgsp=dot(dup,dup);   
    else
        sinal=sinal*sign(dot(dupp1,dup));
        gsp=numgsp/dot(dup,dup);
        dl=sinal*cmp*sqrt(gsp);
    end 
    dupp1=dup;
    dupc1=dup;
else
    dl=-dot(dupc1,dur)/dot(dupc1,dup);
end
end



