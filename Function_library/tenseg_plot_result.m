function tenseg_plot_result(out_tspan,data,legend1,label,name,saveimg)
% This function is written by Ma, S., Chen, M. and Skelton, R., 2022. TsgFEM: Tensegrity finite element method. Journal of Open Source Software, 7(75), p.3390. 
% The code license is Mozilla Public License, v. 2.0.
% The source code is here:
% https://github.com/Muhao-Chen/Tensegrity_Finite_Element_Method_TsgFEM/tree/main/Function_library/tenseg_plot_result.m
%
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% This function plot the result of static or dynamic analysis.
%
% Inputs:
%   data: data to be ploted, can be coordinate, velocity, force, each row
%   is a term changing with time
%	legend: lengend(1),(2)...
%   lable:[xlabel,ylabel];
% Outputs:
%%	plot the results
figure
for i=1:size(data,1)
%      n='rb';
%      plot(out_tspan,data(i,:),n(1,i),'linewidth',1);hold on
       plot(out_tspan,data(i,:),'-.o','linewidth',2);hold on
end
% plot(out_tspan,data,'linewidth',2);
set(gca,'fontsize',18,'linewidth',1.15);
legend(legend1,'location','best','fontsize',15);
ylabel(label(2),'fontsize',18);
xlabel(label(1),'fontsize',18);
if saveimg==1
    saveas(gcf,name);
end

grid on;


end

