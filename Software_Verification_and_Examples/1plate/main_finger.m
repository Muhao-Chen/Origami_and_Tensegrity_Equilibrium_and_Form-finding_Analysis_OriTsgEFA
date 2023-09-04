%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%An Origami finger%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%%
%EXAMPLE
clc; clear; close all;
% Global variable
[consti_data,Eb,Es,sigmab,sigmas,rho_b,rho_s]=material_lib('Steel_Q345','Steel_string');
material{1}='linear_elastic'; % index for material properties: multielastic, plastic.
material{2}=0; % index for considering slack of string (1) for yes,(0) for no (for compare with ANSYS)

% cross section design cofficient
thick=6e-4;        % thickness of hinge

% static analysis set
substep=5;              % load step
lumped=0;               % use lumped matrix 1-yes,0-no
saveimg=0;              % save image or not (1) yes (0)no
savedata=1;             % save data or not (1) yes (0)no
savevideo=1;            % make video(1) or not(0)
gravity=0;              % consider gravity 1 for yes, 0 for no
% move_ground=0;          % for earthquake, use pinned nodes motion(1) or add inertia force in free node(0)
savePath=fullfile(fileparts(mfilename('fullpath')),'data_temp'); %Save files in same folder as this code
%% %% N C of the structure
% Manually specify node positions
width=0.1;
p=3;
N=width*[0:p,0:p;zeros(1,p+1),ones(1,p+1);zeros(1,2*p+2)];       %nodal coordinate

C_in_1=[[(1:p)',(2:p+1)'];[p+1+(1:p)',p+1+(2:p+1)'];[1 p+2];[p+1 2*p+2]];   %bar in boundary
C_in_2=[(1:p)',(p+3:2*p+2)'];      %bar in rigid hinge
C_in_3=[(2:p)',(p+3:2*p+1)'];      %bars in rotational hinge
C_in=[C_in_1;C_in_2;C_in_3];
C = tenseg_ind2C(C_in,N);
[ne,nn]=size(C);        % ne:No.of element;nn:No.of node
C_b=C;C_s=[];
%% define hinge, rigid hinge
% C_in_h is the connectivity of higes, can be written in a function
C_in_h=[C_in_2;C_in_3];
n_h=size(C_in_h,1);         % number of hinge

[~,index_h]=ismember(C_in_h,C_in,'rows');   % index_h is the index number of hinge
[~,index_rh]=ismember(C_in_2,C_in,'rows');   % index_h is the index number of rigid hinge
[~,index_rh_in_h]=ismember(C_in_2,C_in_h,'rows');   % index_rh_in_h is the index of rigid hinge in all hinge

C_h=tenseg_ind2C(C_in(setdiff([1:ne]',index_rh),:),N);     % connectivity matrix of all edges
C_rh=tenseg_ind2C(C_in(index_rh,:),N); % connectivity matrix of rigid edges
% Plot the structure to make sure it looks right
tenseg_plot(N,C_b,C_s);

%% connectivity of triangle element Ca
% Ca can be written in a function
Ca=[[1:p;2:p+1;p+3:2*p+2],[1:p;p+3:2*p+2;p+2:2*p+1]];
% Ca=generate_Ca(C_in,N);
% Ca=zeros(3,1)

[~,np]=size(Ca);        % ne:No.of element;np:No.of plate

% plot the origami configuration
tenseg_plot_ori(N,[],[],C_h,C_rh,[],[],[],[] ,[],Ca);

%% transformation matrix from element to structure
E_n=cell(1,n_h);            %transformation matrix from element node to total node
node_in_hinge=zeros(n_h,4);
I=eye(3*nn);

for i=1:n_h
    node2=C_in_h(i,1);  % start node of the hinge
    node3=C_in_h(i,2);  % end node of the hinge
    for j=1:np
        if (node2==Ca(1,j)&node3==Ca(2,j))|(node2==Ca(2,j)&node3==Ca(3,j))|(node2==Ca(3,j)&node3==Ca(1,j))
            node1=setdiff(Ca(:,j),[node2;node3]);
        elseif (node2==Ca(2,j)&node3==Ca(1,j))|(node2==Ca(3,j)&node3==Ca(2,j))|(node2==Ca(1,j)&node3==Ca(3,j))
            node4=setdiff(Ca(:,j),[node2;node3]);
        end
    end
    node_in_hinge(i,:)=[node1,node2,node3,node4];
    E_n{i}=I(:,kron(node_in_hinge(i,:),3*ones(1,3))-kron(ones(1,4),[2,1,0]));
end
E_n_total=cell2mat(E_n);        % transformation matrix of the whole structure
%% Boundary constraints
pinned_X=[1,2,p+2,p+3]'; pinned_Y=[1,2,p+2,p+3]'; pinned_Z=[1,2,p+2,p+3]';
[Ia,Ib,a,b]=tenseg_boundary(pinned_X,pinned_Y,pinned_Z,nn);

%% generate group index for tensegrity torus structure
gr=[];                      %no group is used
Gp=tenseg_str_gp(gr,C);    %generate group matrix
S=Gp';                      % clustering matrix
%% equilibrium matrix

% equilibrium matrix of truss
[A_1,A_1c,A_1a,A_1ac,A_2,A_2c,A_2a,A_2ac,l,l_c]=tenseg_equilibrium_matrix_CTS(N,C,S,Ia);
% [A_1,A_1g,A_2,A_2g,l,l_gp]=tenseg_equilibrium_matrix2(N,C,Gp,Ia);

% equilibrium matrix of hinge
[phpn_e,phTpn,theta]=jacobian_ori(node_in_hinge,N,E_n_total);       % jacobian matrix
A_o=[A_2,phTpn];
A_o_a=Ia'*A_o;
% A_o_a=A_o;
%% SVD of equilibrium matrix
[U1,U2,V1,V2,S1]=tenseg_svd(A_o_a);         % equilibrium of truss with hinge
% [U1,U2,V1,V2,S1]=tenseg_svd(A_2);           % equilibrium of turss without hinge

%% self-stress design (of truss)
t=zeros(ne,1);      %member force
q=t./l;             % force density
%% self-stress design (of hinge)
M=zeros(n_h,1);
%% cross sectional design (of truss)
A_c=1e-2*ones(ne,1);
E_c=1e9*ones(ne,1);
index_b=[1:ne]';              % index of bar in compression
index_s=setdiff(1:size(S,1),index_b);	% index of strings
%% hinge section design  (of hinge)
k_h=1/12*E_c(index_h).*l(index_h)*thick^3;
k_h(index_rh_in_h)=1e2*1/12*E_c(index_rh).*l(index_rh)*thick^3;      % increase stiffness of rigid hinge
%% rest length (of truss), initial angle (of hinge)
l0_c=l;                     %rest length of truss
theta_0=pi*ones(n_h,1);     % initial angle of hinge
%% tangent stiffness matrix of bars
% [Kt_aa,Kg_aa,Ke_aa,K_mode,k]=tenseg_stiff_CTS(Ia,C,S,q,A_1a,E_c,A_c,l_c);
[Kt_aa,Kg_aa,Ke_aa,K_mode,k]=tenseg_stiff_CTS2(Ia,C,q,A_2ac,E_c,A_c,l0_c);

%% tangent stiffness matrix of hinge

%ph2px2 is the hessian matrix of theta to nodal coordinate
[ph2pn2_e,ph2pn2]=hessian_ori(node_in_hinge,N,E_n);         % calculate hessian matrix
G=cell2mat(ph2pn2);

%% tangent stiffness of the whole origami

K_t_oa=Kt_aa+Ia'*(phTpn*diag(k_h)*phTpn'+G*kron(M,eye(3*nn)))*Ia;

[K_mode,D1] = eig(K_t_oa);         % eigenvalue of tangent stiffness matrix
k=diag(D1);
% plot the mode shape of tangent stiffness matrix
num_plt=1:6;
plot_mode_ori(K_mode,k,N,Ia,[],[],C_h,C_rh,l,'tangent stiffness matrix',...
    'Order of Eigenvalue','Eigenvalue of Stiffness (N/m)',num_plt,0.2,saveimg,3,Ca);
%% mass matrix and damping matrix
rho=1;
mass=rho.*A_c.*l0_c;
M=tenseg_mass_matrix(mass,C,lumped); % generate mass matrix
% damping matrix
d=0;     %damping coefficient
D=d*2*max(sqrt(mass.*E_c.*A_c./l0_c))*eye(3*nn);    %critical damping
%% vibration mode analysis
[V_mode,D1] = eig(K_t_oa,Ia'*M*Ia);         % calculate vibration mode
w_2=diag(D1);                                    % eigen value of
omega=real(sqrt(w_2))/2/pi;                   % frequency in Hz
if 0
    plot_mode_ori(K_mode,k,N,Ia,[],[],C_h,C_rh,l,'natrual vibration',...
        'Order of Vibration Mode','Frequency (Hz)',num_plt,0.2,saveimg,3,Ca);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Nonlinear statics analysis %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% external force, forced motion of nodes, shrink of strings
% calculate external force and
ind_w=[4*3;8*3];w=-5e-3*ones(2,1);   %external force in Z
ind_dnb=[]; dnb0=[];
ind_dl0_c=[]; dl0_c=[];
ind_theta_0=[]; dtheta_0=[];        % initial angel change with time
[w_t,dnb_t,l0_ct,theta_0_t,Ia_new,Ib_new]=tenseg_load_prestress_ori(substep,ind_w,w,ind_dnb,dnb0,ind_dl0_c,dl0_c,ind_theta_0,dtheta_0,theta_0,l0_c,b,gravity,[0;9.8;0],C,mass);

%% Step1: equilibrium calculation
% input data
data.N=N; data.C=C; data.ne=ne; data.nn=nn; data.Ia=Ia_new; data.Ib=Ib_new;data.S=S;
data.E=E_c; data.A=A_c; data.index_b=index_b; data.index_s=index_s;
data.consti_data=consti_data;   data.material=material; %constitue info
data.w_t=w_t;  % external force
data.dnb_t=dnb_t;           % forced movement of pinned nodes
data.l0_t=l0_ct;            % forced change of rest length
data.theta_0_t=theta_0_t;   % forced change of initial angle
data.k_h=k_h;               % stiffness of hinge
data.E_n=E_n;               % transfer matrix from matrix to structure
data.node_in_hinge=node_in_hinge;       % node in triangle element in hinge
data.substep=substep;    % substep
data.InitialLoadFactor=0.001;
data.MaxIcr=1000;
data.LoadType='Force'; % 'Force' or 'Displacement'
data.StopCriterion=@(U)(norm(U)>0.5);

% nonlinear analysis
data_out1=static_solver_ori_2(data);

Fhis=data_out1.Fhis;          %load factor
t_t=data_out1.t_out;          %member force in every step
n_t=data_out1.n_out;          %nodal coordinate in every step
M_out=data_out1.M_out;
l_out=data_out1.l_out;
icrm=size(n_t,2);               % increment
clear ans;                      % don't save figure handle
%% save output data
if savedata==1
     save (fullfile(savePath,['Plate_all_',material{1},'.mat']));
end
% N_out=data_out1.N_out;
%% Plot final configuration
j=linspace(1e-5,1,4);
for i=1:4
    num=ceil(j(i)*size(n_t,2));
    tenseg_plot_ori(reshape(n_t(:,num),3,[]),[],[],C_h,C_rh,[],[],[30,30],[] ,[],Ca);
    %  axis off;
end
%% plot member force
tenseg_plot_result(Fhis,t_t,{},{'Load factor','Force / N'},fullfile(savePath,'plot_member_force.png'),saveimg);
%% plot member length
tenseg_plot_result(Fhis,l_out,{},{'Load factor','length / m'},fullfile(savePath,'plot_member_length.png'),saveimg);

%% plot hinge moment
tenseg_plot_result(Fhis,M_out,{'rgd1','rgd2','rgd3','sft1','sft2'},{'Load factor','Moment / N \times m'},fullfile(savePath,'plot_hinge_moment.png'),saveimg);

%% Plot nodal coordinate curve X Y
% tenseg_plot_result(Fhis,n_t([4*3,8*3],:),{'4Z','8Z'},{'Substep','Coordinate /m)'},fullfile(savePath,'plot_coordinate.png'),saveimg);
%% make video of the dynamic
name=fullfile(savePath,'origami_finger');
% name=['origami_finger'];
% tenseg_video(n_t,C_b,C_s,[],min(substep,50),name,savevideo,R3Ddata);
% tenseg_video_slack(n_t,C_b,C_s,l0_ct,index_s,[],[],[],min(substep,50),name,savevideo,material{2})
tenseg_video_ori(n_t,[],[],C_h,C_rh,Ca,[],min(icrm,50),name,savevideo,[]);