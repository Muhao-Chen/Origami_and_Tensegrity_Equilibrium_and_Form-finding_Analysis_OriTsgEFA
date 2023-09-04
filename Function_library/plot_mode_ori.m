function plot_mode_ori(Mode,value,N,Ia,C_b,C_s,C_h,C_rh,l,title,xlb,ylb,num_plt,ampli,saveimg,view_vec,Ca)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% This function plot the mode shape(eigenvalue) of origami-tensegrity structure
% with eigenvalue
%
%Inputs:
%   Mode: eigenvcector of the whole structure tangent stiffness matrix
%   value: eigenvalue of the whole structure tangent stiffness matrix
%   N: node matrix (3 x n array for n nodes)
%   Ia: transfer matrix of free coordinates
%   C_b (optional): bar connectivity matrix (beta x n array for beta bars)
%   C_s (optional): string connectivity matrix (alpha x n array for alpha strings)
%   C_h: hinge connectivity matrix
%   C_rh: rigid hinge connectivity matrix
%   l: members length vector
%   title: the title saved image
%   xlb: tangent stiffness image xlabel
%   ylb: tangent stiffness image ylabel
%   num_plt: number of modes to be ploted
%   ampli: coefficient of amplification for deformed shape from initial shape
%   saveimg: save image or not (1) yes (0)no
%   view_vec (optional): plotting view direction (see view())
%   Ca: connectivity matrix of triangle elements
%   Example: plot_mode_ori(K_mode,k,N,Ia,[],[],C_h,C_rh,l,'tangent stiffness matrix',...
%    'Order of Eigenvalue','Eigenvalue of Stiffness (N/m)',num_plt,0.2,saveimg,3,Ca);
%%
switch nargin
    case 13+2
        view_vec=2;
end

figure
% plot(1:numel(value),value,'k-o','linewidth',1.5); % original plot
semilogy(1:numel(value),value,'k-o','linewidth',1.5); %semilogy
set(gca,'fontsize',18);
xlabel(xlb,'fontsize',18,'Interpreter','latex');
ylabel(ylb,'fontsize',18,'Interpreter','latex');
grid on;
if saveimg==1
    saveas(gcf,[title,'.png']);
end

for i=1:numel(num_plt)
    f1=figure;
    title2=({['Mode ',num2str(num_plt(i)),'  \lambda_{',num2str(num_plt(i)),'} = ',num2str(value(num_plt(i)),'%.1e')]});
    % title2=[];
    %plot buckling mode
    tenseg_plot_ori(N+ampli*max(l)*reshape(Ia*Mode(:,num_plt(i))/norm(Mode(:,num_plt(i))),3,[]),C_b,C_s,C_h,C_rh,f1,[],[],title2,[],Ca);
    tenseg_plot_ori_dash(N,C_b,C_s,C_h,C_rh,f1,[],[],title2,[],Ca);

    axis off;
    view(view_vec);
    if saveimg==1
        saveas(gcf,[title,' of ',num2str(num_plt(i)),'.png']);
    end
end

end

