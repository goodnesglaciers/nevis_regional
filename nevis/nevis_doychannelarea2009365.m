function [channel_area, channel_area_percentage] = nevis_doychannelarea2009365(zzall_nevis)
%% Find the percentage in the terrasarx plot that is channelized (N>0)
% 7 June 2017 LAS

% OUTPUT
% channel_area = square km area that is channelized (N>0) on given DOY

%% get into TSX grid
load('/nevis/terravel.mat')
load('nevis_170207b.mat');
load('nevis_122122/0185.mat');
%get rid of points outside domain
nx(gg.nout) = NaN; ex(gg.eout) = NaN; fx(gg.fout) = NaN; cx(gg.cout) = NaN;

%% get bed topo and ice sheet surface data
xx = (ps.x/10^3)*gg.nx; xxx = reshape(xx, size(xx,1)*size(xx,2), 1);
yy = (ps.x/10^3)*gg.ny; yyy = reshape(yy, size(yy,1)*size(yy,2), 1);

% interpolate nevis to terrasarX grid
terraxx = (terravel.winter.vel.*(terravel.polar.xx_moulin./1000))./(terravel.winter.vel);
terrayy = (terravel.winter.vel.*(terravel.polar.yy_moulin./1000))./(terravel.winter.vel);

for i=1:1:365
    disp(['First interp ',num2str(i)]);
    F = scatteredInterpolant(xxx,yyy,reshape(zzall_nevis(i).N,size(xx,1)*size(xx,2), 1)); 
    N_ter(:,:,i) = F(terraxx, terrayy);
    
    F = scatteredInterpolant(xxx,yyy,reshape(zzall_nevis(i).Q,size(xx,1)*size(xx,2), 1)); 
    Q_ter(:,:,i) = F(terraxx, terrayy);
end

%% rotate into a size that works for everyone...MUST BE SQUARE, No NaNs
% shift couterclockwise wise in degrees
theta = -10.4*(pi/180); 
R = [cos(theta) -sin(theta); sin(theta) cos(theta)]; % rotation matrix
v = [reshape(terraxx, size(terraxx,1)*size(terraxx,2), 1)';...
    reshape(terrayy, size(terraxx,1)*size(terraxx,2), 1)'];
so = (R*v)'; 
SO = so; SO(any(isnan(SO),2),:) = []; % remove NaNs
SO = sortrows(so,1);  % sort based on x value
% create new terrasarx X and Y meshgrid from the rotated positions at 1 km resolution
[txx, tyy] = meshgrid(min(SO(:,1))+1:0.1:max(SO(:,1)), min(SO(:,2)):0.1:max(SO(:,2))); 

for i=1:1:365
    disp(['Second interp ',num2str(i)]);
    all = horzcat(so, reshape(N_ter(:,:,i), size(terraxx,1)*size(terraxx,2), 1)); % rotated
    all(any(isnan(all),2),:) = []; % remove NaNs
    F = scatteredInterpolant(all(:,1),all(:,2),all(:,3));   
    % rotated N
    N_rot(:,:,i) = F(txx, tyy);

    all = horzcat(so, reshape(Q_ter(:,:,i), size(terraxx,1)*size(terraxx,2), 1)); % rotated
    all(any(isnan(all),2),:) = []; % remove NaNs
    FQ = scatteredInterpolant(all(:,1),all(:,2),all(:,3));   
    % rotated N
    Q_rot(:,:,i) = FQ(txx, tyy);
end

%% Find area that is channelized. N>0.
    for i=1:1:365
    disp(['Channelize condition ',num2str(i)]);
    % N condition
    B = 0; N_channel = gt(N_rot(:,:,i),B);     % logical 1 if greater than 0
    % q_critical condition: q>q_critical = 0.001 m2/s
    q_critical = 0.0001; Q_channel = gt(Q_rot(:,:,i),q_critical); % logical 1 if q>q_critical
    % meets both conditions
    NQ = double(N_channel) + double(Q_channel);
    channel_area(i) = (0.9*0.9)*length(find(NQ==2)); % km^2 area
    channel_area_percentage(i) = 100*channel_area(i)/(0.9*0.9*552*309); 
    end

end