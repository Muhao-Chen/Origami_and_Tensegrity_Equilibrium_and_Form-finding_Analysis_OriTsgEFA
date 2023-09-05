%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%An simple miura pattern%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%%
%EXAMPLE
clc; clear; close all;
% Global variable
[consti_data,Eb,Es,sigmab,sigmas,rho_b,rho_s]=material_lib('Plastic','Plastic');
material{1}='linear_elastic'; % index for material properties: multielastic, plastic.
material{2}=0; % index for considering slack of string (1) for yes,(0) for no (for compare with ANSYS)

% cross section design cofficient
thick=6e-4;             % thickness of hinge

% static analysis set
substep=40;             % load step
lumped=0;               % use lumped matrix 1-yes,0-no
saveimg=0;              % save image or not (1) yes (0)no
savedata=0;             % save data or not (1) yes (0)no
savevideo=0;            % make video(1) or not(0)
gravity=0;              % consider gravity 1 for yes, 0 for no
% move_ground=0;          % for earthquake, use pinned nodes motion(1) or add inertia force in free node(0)
savePath=fullfile(fileparts(mfilename('fullpath')),'data_temp'); %Save files in same folder as this code
%% %% N C of the structure
% Manually specify node positions
% N=[0 1 0;1 1 0;0 0 0;1 0 0;0.7 0.3 0]';
T=[cosd(45) -sind(45) 0;sind(45) cosd(45) 0; 0 0 1];
N=T*[-1 1 0;sqrt(2)*[cosd(60), sind(60) 0;-sind(60) -cosd(60) 0];1 -1 0;0 0 0]';

% Manually specify connectivity indices.
C_b_in_v=[1 2;2 4;4 3;3 1];% vertical horizontal bar
C_b_in_d=[1 5;2 5;4 5;3 5];% diagonal bar
C_in = [C_b_in_v;C_b_in_d];% This is indicating the bar connection
% Convert the above matrices into full connectivity matrices.
C = tenseg_ind2C(C_in,N);
[ne,nn]=size(C);        % ne:No.of element;nn:No.of node
C_b=C;C_s=[];
% Plot the structure to make sure it looks right
tenseg_plot(N,C_b,[]);

%% define hinge, rigid hinge
% C_in_h is the connectivity of higes, can be written in a function
C_in_h=[C_b_in_d];
n_h=size(C_in_h,1);         % number of hinge

[~,index_h]=ismember(C_in_h,C_in,'rows');   % index_h is the index number of hinge
index_rh=[];index_rh_in_h=[];
% [~,index_rh]=ismember(C_in_2,C_in,'rows');   % index_h is the index number of rigid hinge
% [~,index_rh_in_h]=ismember(C_in_2,C_in_h,'rows');   % index_h is the index of rigid hinge in all hinge

C_h=tenseg_ind2C(C_in(setdiff([1:ne]',index_rh),:),N);     % connectivity matrix of all edges
C_rh=tenseg_ind2C(C_in(index_rh,:),N); % connectivity matrix of rigid edges
% Plot the structure to make sure it looks right
tenseg_plot(N,C_b,C_s);

%% connectivity of triangle element Ca
% Ca can be written in a function
Ca=[5 2 1;5 4 2;5 3 4;5 1 3]';% Ca=generate_Ca(C_in,N);
% Ca=zeros(3,1)

[~,np]=size(Ca);        % ne:No.of element;np:No.of plate

% plot the origami configuration
tenseg_plot_ori(N,[],[],C_h,C_rh,[],[],[],[] ,[],Ca);
axis off

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
    %  E_ntt{i}=E_n{i}(1:3:end,1:3:end);
end

E_n_total=cell2mat(E_n);        % transformation matrix of the whole structure
%% Boundary constraints
pinned_X=[1 2 5]'; pinned_Y=[1 2 5]'; pinned_Z=[1 2 5]';
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
A_c=1e-4*ones(ne,1);
E_c=1e6*ones(ne,1);
index_b=[1:ne]';              % index of bar in compression
index_s=setdiff(1:size(S,1),index_b);	% index of strings
%% hinge section design  (of hinge)
% k_h=1/12*E_c(index_h).*l(index_h)*thick^3;
% k_h(index_rh_in_h)=1e2*1/12*E_c(index_rh).*l(index_rh)*thick^3;      % increase stiffness of rigid hinge
k_h=1/12*E_c(index_h).*l(index_h)*thick^3;
%% rest length (of truss), initial angle (of hinge)
l0_c=l;                     %rest length of truss
theta_0=theta;     % initial angle of hinge (give different value)
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
[k, ind] = sort(k);
K_mode = K_mode(:, ind);

% plot the mode shape of tangent stiffness matrix
num_plt=1:6;
plot_mode_ori(round(K_mode,12),k,N,Ia,[],[],C_h,C_rh,l,'tangent stiffness matrix',...
    'Order','Eigenvalue (N/m)',num_plt,0.15,saveimg,[0,30],Ca);

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
ind_w=[];w=[];   %external force in Z
ind_dnb=[]; dnb0=[];
ind_dl0_c=[]; dl0_c=[];
ind_theta_0=[1:4]'; dtheta_0=pi*0.99*[1 -1 -1 -1]';        % initial angel change with time
[w_t,dnb_t,l0_ct,theta_0_t,Ia_new,Ib_new]=tenseg_load_prestress_ori(substep,ind_w,w,ind_dnb,dnb0,ind_dl0_c,dl0_c,ind_theta_0,dtheta_0,theta_0,l0_c,b,gravity,[0;9.8;0],C,mass);

% add initial defect by first mode shape
N_d=N+0.0*max(l)*reshape(Ia*K_mode(:,1),3,[]);
%% Step1: equilibrium calculation
% input data
data.N=N_d; data.C=C; data.ne=ne; data.nn=nn; data.Ia=Ia_new; data.Ib=Ib_new;data.S=S;
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
data.LoadType='Substep'; % 'Force' or 'Displacement' or 'Substep'
data.StopCriterion=@(U)(norm(U)>20);


% nonlinear analysis
data_out1=static_solver_ori_2(data);

Fhis=data_out1.Fhis;          %load factor
t_t=data_out1.t_out;          %member force in every step
n_t=data_out1.n_out;          %nodal coordinate in every step
M_out=data_out1.M_out;
l_out=data_out1.l_out;
Kt_aa_out= data_out1.Kt_aa_out;       %tangent stiffness of truss
K_t_oa_out= data_out1.K_t_oa_out;       %tangent stiffness of whole struct.
theta_t=data_out1.theta_out;

icrm=size(n_t,2);               % increment
% N_out=data_out1.N_out;
clear ans;                      % don't save figure handle
%% save output data
if savedata==1
    save (fullfile(savePath,['Mirua_all_',material{1},'.mat']));
end
%% Large deformation analysis
%% plot member force
tenseg_plot_result(Fhis,t_t,{},{'Load factor','Force / N'},fullfile(savePath,'plot_member_force.png'),saveimg);

%% plot angle, rest angle, hinge moment

%plot moment
tenseg_plot_result2(Fhis,M_out,{'$M_1$','$M_2$','$M_3$','$M_4$'},{'Load factor','Moment / N \times m'},...
    fullfile(savePath,'plot_hinge_moment.png'),0,{'-or','-.vg','--xb',':^k'});
legend('NumColumns',2,'Interpreter','latex');
fig=gcf;
fig.Position(3:4)=[800,350];   %change fig size
if saveimg==1
    saveas(gcf,fullfile(savePath,'plot_hinge_moment.png'));
end

% plot angle, rest angle
tenseg_plot_result2(Fhis,180/pi*[theta_0_t;theta_t],{'$\theta_{0,1}$','$\theta_{0,2}$','$\theta_{0,3}$','$\theta_{0,4}$','$\theta_{1}$','$\theta_{2}$','$\theta_{3}$','$\theta_{4}$'}...
    ,{'Load factor','Folding Angle ({\circ})'},fullfile(savePath,'plot_hinge_angle.png'),0,{'-r','-.g','--b',':k','-or','-.vg','--xb',':^k'});
legend('NumColumns',2,'Interpreter','latex');
fig=gcf;
fig.Position(3:4)=[800,350];
if saveimg==1
    saveas(gcf,fullfile(savePath,'plot_hinge_angle.png'));
end
%% plot stiffness in large deformation in XYZ direction

F_dir=zeros(3*nn,3);
F_dir([4*3-2:4*3],:)=eye(3);   % force with direction X Y Z

compliance_dir=zeros(3,substep);    % compliance with direction X Y Z
stiff_dir=zeros(3,substep);         % stiff with direction X Y Z

for i=1:substep
    K_t_oa=K_t_oa_out{i};
    disp_dir=K_t_oa\(Ia'*F_dir);
    compliance_dir(:,i)=diag(disp_dir'*K_t_oa'*disp_dir);
    stiff_dir(:,i)=1./compliance_dir(:,i);
end

%plot stiffness
figure
semilogy(Fhis,stiff_dir(1,:),'-r',...
    Fhis,stiff_dir(2,:),'-.g',...
    Fhis,stiff_dir(3,:),'--b','linewidth',2); %semilogy
set(gca,'fontsize',18);
xlabel('Load factor','fontsize',18,'Interpreter','tex');
ylabel('Stiffness (N/m)','fontsize',18);
lgd =legend('X','Y','Z');
legend('NumColumns',3);
title(lgd,'Direction')
grid on;
fig=gcf;
fig.Position(3:4)=[800,350];   %change fig size
if saveimg==1
    saveas(gcf,fullfile(savePath,'plot_stiffness.png'));
end

%% Plot final configuration
num_t=6;
j=linspace(1e-5,1,num_t);

for i=1:num_t
    hf=figure;
    num=ceil(j(i)*size(n_t,2));
    %     tenseg_plot(N,C_b,C_s,hf);
    tenseg_plot_ori(reshape(n_t(:,num),3,[]),[],[],C_h,C_rh,hf,[],[45,30],[] ,[],Ca);
    tenseg_plot_ori_dash(N,C_b,C_s,C_h,C_rh,hf,[],[45,30],[],[],Ca);

    xlim([min(n_t(1:3:end,:),[],'all'),max(n_t(1:3:end,:),[],'all')]);
    ylim([min(n_t(2:3:end,:),[],'all'),max(n_t(2:3:end,:),[],'all')]);
    zlim([min(n_t(3:3:end,:),[],'all'),max(n_t(3:3:end,:),[],'all')]);
    axis off;
end
%% make video of the dynamic
name=fullfile(savePath,'A_Mirua');
% name=['miura_ori_1'];
% tenseg_video(n_t,C_b,C_s,[],min(substep,50),name,savevideo,R3Ddata);
% tenseg_video_slack(n_t,C_b,C_s,l0_ct,index_s,[],[],[],min(substep,50),name,savevideo,material{2})
tenseg_video_ori(n_t,[],[],C_h,C_rh,Ca,[],min(icrm,50),name,savevideo,[])
%% Inifinitesimal deformation analysis
percent=0.85; num=round(percent*substep);

Kt_aa=Kt_aa_out{num};       %tangent stiffness of truss
K_t_oa=K_t_oa_out{num};     %tangent stiffness of whole struct.

num_tot=9e3;
dir_angle=[linspace(0,5/6*pi,0.1*num_tot),linspace(5/6*pi,7/6*pi,0.4*num_tot),...
    linspace(7/6*pi,11/6*pi,0.1*num_tot),linspace(11/6*pi,12/6*pi,0.4*num_tot)]';            %direction
% dir_angle=linspace(0,2*pi,num_tot)';            %direction
F_dir=zeros(3*nn,numel(dir_angle));
F_dir([4*3-1,4*3],:)=[cos(dir_angle)';sin(dir_angle)'];   % force with direction
F_dir=F_dir/diag(sqrt(diag(F_dir'*F_dir))); %normalized
%displacement
disp_dir=K_t_oa\(Ia'*F_dir);
% disp_dir=Ia'*F_dir;
% complaint coefficient
complaint_dir_ori=diag(disp_dir'*K_t_oa'*disp_dir);     %for ori
complaint_dir_truss=diag(disp_dir'*Kt_aa'*disp_dir);     %for truss
complaint_dir_hinge=diag(disp_dir'*(K_t_oa-Kt_aa)'*disp_dir);     %for hinge
% stiffness coefficient
stiff_dir_ori=1./complaint_dir_ori;
stiff_dir_truss=stiff_dir_ori.*complaint_dir_truss./complaint_dir_ori;
stiff_dir_hinge=stiff_dir_ori.*complaint_dir_hinge./complaint_dir_ori;
% plot stiffness in Bode
figure
semilogy(dir_angle*180/pi,stiff_dir_ori,'r-',...
    dir_angle*180/pi,stiff_dir_truss,'k--',...
    dir_angle*180/pi,stiff_dir_hinge,'c-.','linewidth',1.5); %semilogy
set(gca,'fontsize',18);
xlabel('Angle \alpha ({\circ})','fontsize',18,'Interpreter','tex');
ylabel('Stiffness (N/m)','fontsize',18,'Interpreter','latex');
legend('Origami','Bars','Hinges')
xlim([0 360]);
grid on;
% plot percentage
figure
plot(dir_angle*180/pi,stiff_dir_ori./stiff_dir_ori,'r-',...
    dir_angle*180/pi,stiff_dir_truss./stiff_dir_ori,'k--',...
    dir_angle*180/pi,stiff_dir_hinge./stiff_dir_ori,'c-.','linewidth',1.5); %semilogy
set(gca,'fontsize',18);
xlabel('Angle \alpha ({\circ})','fontsize',18,'Interpreter','tex');
ylabel('Percentage of Stiffness','fontsize',18,'Interpreter','latex');
legend('Origami','Bars','Hinges')
xlim([0 360]);
grid on;



