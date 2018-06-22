%% FIGURE 3 PARAMETER SPACE BASED ON AREA OF EFFICIENT DRAINAGE
%% calculate daily efficient drainage area for a single run
% 8 June 2017
clear all; 

%% 2009 
% Distributed 22212 2009
% load('/fuxi_over5/nevis_22212_zzall_nevis.mat');
% [channel_area, channel_area_percentage] = nevis_doychannelarea2009365(zzall_nevis);

load('/fuxi_over5/nevis_22222_zzall_nevis.mat');
    [channel_area2, channel_area_percentage2] = nevis_doychannelarea2009365(zzall_nevis);

% % Discrete 22222 2009
% Files2 = dir(fullfile('fuxi_over5','nevis_22222_gammasqrVp.mat')); % load files
% Files22 = dir(fullfile('fuxi_over5','nevis_22222_gammasqrVp_Q.mat')); % load files
% for i=1:length(DOI_2009)
%     [channel_area2(i), channel_area_percentage2(i)] = ...
%         nevis_doychannelarea2009(Files2.name,Files22.name,J2009_2(i));
% end

%%
figure(1)
plot(1:1:365, channel_area_percentage2, 'r');
