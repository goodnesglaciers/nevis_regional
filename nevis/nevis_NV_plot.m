function nevis_NV_plot(fn,tis,fnA);
% plot timeseries of GPS velocity and effective pressure at the Flowline
% stations for 2011-2013 years
% 17 August 2016
% Inputs
%   fn            filename
%   tis           time steps for plotting
%   fnA           filename for GPS velocities mat file

%% load initial timestep
if exist([fn,'/0000']), load([fn,'/0000']); else load(fn); end
load([fn,'/0365']);
str = sprintf('%s.mat',[fn]); load(str)
str_save = sprintf('%s_transect',[oo.root,oo.fn]);
if isfield(oo,'reversey'), reversey = oo.reversey; else reversey = 0; end

%% load Joughin 2013 The Cryosphere TerraSAR-X footprint and Flowline Stations
radius=6378137.0; eccen=0.08181919; lat_true=70; lon_posy=-45; % projection parameters
foot = [68.884795, -50.182106; 68.952545, -49.409098;
        %68.393966, -48.867272; 68.282561, -49.831050];
        68.469744, -48.970285; 68.399076, -49.709055];
moulin = [68.723585, -49.536195]; % moulin is (0,0)    
[moulin_x,moulin_y] = polarstereo_fwd(moulin(1),moulin(2),radius,eccen,lat_true,lon_posy);
[foot_x,foot_y] = polarstereo_fwd(foot(:,1),foot(:,2),radius,eccen,lat_true,lon_posy);
footprint(:,1) = foot_x-moulin_x; footprint(:,2) = foot_y - moulin_y;

FLXX = [68.7735, 310.0752-360;
        68.7559, 310.2587-360;
        68.7373, 310.4469-360;
        68.7253, 310.5542-360;
        68.7158, 310.8561-360;
        68.7192, 311.1320-360];
[flxx_x,flxx_y] = polarstereo_fwd(FLXX(:,1),FLXX(:,2),radius,eccen,lat_true,lon_posy);
flxx_plot(:,1) = flxx_x-moulin_x; flxx_plot(:,2) = flxx_y - moulin_y;    

%% load velocity file
str = sprintf('%s.mat',[fnA]); vels = load(str); vels = vels.(fnA);
  
% %% extract initial variables
% if isfield(vv,'nbdy'), gg = nevis_label(gg,vv.nbdy); aa.phi = aa.phi_a(gg.nbdy); end % re-mask
% oo.evaluate_variables = 1; oo.evaluate_residual = 1; oo.evaluate_jacobian = 0; [vv2] = nevis_backbone(inf,vv,vv,aa,pp,gg,oo); % expand solution variables
% vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo); % calculate node discharge
% %nevis_unpack(aa,gg);

% %get rid of points outside domain
% nx(gg.nout) = NaN;
% ex(gg.eout) = NaN;
% fx(gg.fout) = NaN;
% cx(gg.cout) = NaN;
    
for i_t=1:1:365;
%% load timestep
%str = sprintf('%s.mat',[fn]); load(str)
load(['nevis_160813a','/',int2four(i_t)]);
      
% %% extract new variables for location specific E
if isfield(vv,'nbdy'), gg = nevis_label(gg,vv.nbdy); aa.phi = aa.phi_a(gg.nbdy); end % remask
aa = nevis_inputs(vv.t,aa,pp,gg,oo);
oo.evaluate_variables = 1; oo.evaluate_residual = 1; oo.evaluate_jacobian = 0; [vv2] = nevis_backbone(inf,vv,vv,aa,pp,gg,oo); % expand solution variables
vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo); % calculate node discharge
nevis_unpack(aa,vv2);
%clear tt vv vv2

% interpolate E to FLXX locations
E_grid = reshape(E, gg.nI, gg.nJ); 
xx = (ps.x/10^3)*reshape(gg.nx,gg.nI*gg.nJ,1); 
yy = (ps.x/10^3)*reshape(gg.ny,gg.nI*gg.nJ,1); 
F = scatteredInterpolant(xx,yy,E);
E_FLXX(:,i_t) = F(flxx_plot(:,1)./10^3, flxx_plot(:,2)./10^3);

end

for i_t=1:1:365;
%% load timestep
%str = sprintf('%s.mat',[fn]); load(str)
load([fn,'/',int2four(i_t)]);
      
% %% extract new variables for location specific E
if isfield(vv,'nbdy'), gg = nevis_label(gg,vv.nbdy); aa.phi = aa.phi_a(gg.nbdy); end % remask
aa = nevis_inputs(vv.t,aa,pp,gg,oo);
oo.evaluate_variables = 1; oo.evaluate_residual = 1; oo.evaluate_jacobian = 0; [vv2] = nevis_backbone(inf,vv,vv,aa,pp,gg,oo); % expand solution variables
vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo); % calculate node discharge
nevis_unpack(aa,vv2);
%clear tt vv vv2

% save N at every point
N_all = (ps.phi/10^6)*reshape(phi_0-phi,gg.nIJ,1);
N_time(:,i_t) = N_all;

end

load([fn,'/0365']);
    t = (ps.t/(24*60*60))*[tt.t];
    Q_in = ps.Q*[tt.Q_in];
    Q_out = ps.Q*[tt.Q_out];
    Q_outQ = ps.Q*[tt.Q_outQ];
    Q_outq = ps.Q*[tt.Q_outq];
    m = (ps.m*ps.x^2)*[tt.m];
    %E = (ps.m*ps.x^2)*[tt.E];

    phi = (ps.phi/10^6)*[tt.phi];
    N = (ps.phi/10^6)*[tt.N];
    pw = (ps.phi/10^6)*([tt.phi]-aa.phi_a(1));
    hs = ps.x^2*ps.h*[tt.hs];               % total cavity sheet volume scaled
    he = ps.x^2*ps.h*[tt.he];               % total hel volume scaled
    S = ps.x*ps.S*[tt.S];                   % total channel vilume scaled
    A = ps.x^2*sum(gg.Dx.*gg.Dy);           % area of the domain
    total_water = S + hs + he;              % total volume of water in domain (hcav+hel+channel)
   
    if isfield(tt,'pts_phi')    
    pts_phi = (ps.phi/10^6)*[tt.pts_phi];
    pts_hs = ps.hs*[tt.pts_hs];
    pts_he = ps.he*[tt.pts_he];
    pts_N = (ps.phi/10^6)*[tt.pts_N];
    %pts_N = (ps.phi/10^6)*(aa.phi_0(oo.pts_ni)*[tt.t].^0 - [tt.pts_phi]);
    pts_pw = (ps.phi/10^6)*([tt.pts_phi]-aa.phi_a(oo.pts_ni)*[tt.t].^0);
    pts_prat = ([tt.pts_phi]-aa.phi_a(oo.pts_ni)*[tt.t].^0)./(aa.phi_0(oo.pts_ni)*[tt.t].^0-aa.phi_a(oo.pts_ni)*[tt.t].^0);
    pts_S = ps.x*ps.S*([tt.pts_Sx]+[tt.pts_Sy]+[tt.pts_Sr]+[tt.pts_Ss]);
    pts_A = ps.x^2*4*(gg.Dx(1)*gg.Dy(1));
    lc = pd.l_c;
    ks = pd.k_s;
    end
    
%% get daily average for N, V, and E
% E daily average is E_FLXX from above
% N daily average
day = [tt.t].*10; 
for i=2:1:365;
    [ind,ing] = find(day>= (i-1) & day < i);
    pts_N_daymean(:,i) = mean(pts_N(:,ing),2);
end

% V daily average
for i=2:1:365
    [ind,ing] = find(vels(61,:)>= (i-1) & vels(61,:) < i);
    pts_V_daymean(1,i) = nanmean(vels(1,ing));
    pts_V_daymean(2,i) = nanmean(vels(2,ing));
    pts_V_daymean(3,i) = nanmean(vels(3,ing));
    pts_V_daymean(4,i) = nanmean(vels(4,ing));
    pts_V_daymean(5,i) = nanmean(vels(5,ing));
    pts_V_daymean(6,i) = nanmean(vels(6,ing));
end

%% set up figure
fig2=figure(2); clf; set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',[5 2 15 20]);
set(0,'DefaultAxesFontSize',10);
set(0,'DefaultTextFontSize',10);
axes_size = [0.1 0.1 0.8 0.8];
axe1 = axes('Position',[ 0.1 0.81 0.8 0.14 ],'Box','on','NextPlot','add'); 
axe2 = axes('Position',[ 0.1 0.66 0.8 0.14 ],'Box','on','NextPlot','add'); 
axe3 = axes('Position',[ 0.1 0.51 0.8 0.14 ],'Box','on','NextPlot','add');
axe4 = axes('Position',[ 0.1 0.36 0.8 0.14 ],'Box','on','NextPlot','add');
axe5 = axes('Position',[ 0.1 0.21 0.8 0.14 ],'Box','on','NextPlot','add');
axe6 = axes('Position',[ 0.1 0.06 0.8 0.14 ],'Box','on','NextPlot','add');

axes(axe1)
vel_test = 100.*ones(33798,1);
[AX,H1,H2] = plotyy([tt.t].*10, -1.*pts_N(1,:), vels(61,:), vels(1,:),'plot');
set(H1,'linestyle','-','marker','none'); 
set(H2,'linestyle','none','marker','.','MarkerSize',3);
set(AX(1),'XTickLabel',[],'ylim',[-2.5 2.5],'xlim',[0 365],'ytick',[-2:1:2]); 
set(AX(2),'XTickLabel',[],'ylim',[50 300],'xlim',[0 365],'ytick',[50:50:250]);
hold on; plot([1:365],(E_FLXX(1,:)./1.25)-2.5,'k');
text(5,2,'FL01','FontSize',12);
title([fn,'   ',fnA],'FontSize',14)
text(0,3.75,['K_s = ',num2str(ks,'%2.1e')]);
text(0,3,['lambda_c = ',num2str(pd.l_c),' m']);

axes(axe2)
vel_test = 100.*ones(33798,1);
[AX,H1,H2] = plotyy([tt.t].*10, -1.*pts_N(2,:), vels(61,:), vels(2,:),'plot');
set(H1,'linestyle','-','marker','none'); 
set(H2,'linestyle','none','marker','.','MarkerSize',3);
set(AX(1),'XTickLabel',[],'ylim',[-2.5 2.5],'xlim',[0 365],'ytick',[-2:1:2]); 
set(AX(2),'XTickLabel',[],'ylim',[50 300],'xlim',[0 365],'ytick',[50:50:250]);
hold on; plot([1:365],(E_FLXX(2,:)./1.25)-2.5,'k');
text(5,2,'FL02','FontSize',12);

axes(axe3)
vel_test = 100.*ones(33798,1);
[AX,H1,H2] = plotyy([tt.t].*10, -1.*pts_N(3,:), vels(61,:), vels(3,:),'plot');
set(H1,'linestyle','-','marker','none'); 
set(H2,'linestyle','none','marker','.','MarkerSize',3);
set(AX(1),'XTickLabel',[],'ylim',[-2.5 2.5],'xlim',[0 365],'ytick',[-2:1:2]); 
set(AX(2),'XTickLabel',[],'ylim',[50 300],'xlim',[0 365],'ytick',[50:50:250]);
set(get(AX(1),'Ylabel'),'String','-1*N [ MPa ]')
set(get(AX(2),'Ylabel'),'String','Velocity [ m yr^{-1} ]')
hold on; plot([1:365],(E_FLXX(3,:)./1.25)-2.5,'k');
text(5,2,'FL03','FontSize',12);

axes(axe4)
vel_test = 100.*ones(33798,1);
[AX,H1,H2] = plotyy([tt.t].*10, -1.*pts_N(5,:), vels(61,:), vels(4,:),'plot');
set(H1,'linestyle','-','marker','none'); 
set(H2,'linestyle','none','marker','.','MarkerSize',3);
set(AX(1),'XTickLabel',[],'ylim',[-2.5 2.5],'xlim',[0 365],'ytick',[-2:1:2]); 
set(AX(2),'XTickLabel',[],'ylim',[50 300],'xlim',[0 365],'ytick',[50:50:250]);
hold on; plot([1:365],(E_FLXX(4,:)./1.25)-2.5,'k');
text(5,2,'FL04','FontSize',12);

axes(axe5)
vel_test = 100.*ones(33798,1);
[AX,H1,H2] = plotyy([tt.t].*10, -1.*pts_N(6,:), vels(61,:), vels(5,:),'plot');
set(H1,'linestyle','-','marker','none'); 
set(H2,'linestyle','none','marker','.','MarkerSize',3);
set(AX(1),'XTickLabel',[],'ylim',[-2.5 2.5],'xlim',[0 365],'ytick',[-2:1:2]); 
set(AX(2),'XTickLabel',[],'ylim',[50 300],'xlim',[0 365],'ytick',[50:50:250]);
%plot([tt.t].*10, (E./(1e3))-2.5,'k')
hold on; plot([1:365],(E_FLXX(5,:)./1.25)-2.5,'k');
text(5,2,'FL05','FontSize',12);

axes(axe6)
vel_test = 100.*ones(33798,1);
[AX,H1,H2] = plotyy([tt.t].*10, -1.*pts_N(7,:), vels(61,:), vels(6,:),'plot');
set(H1,'linestyle','-','marker','none'); 
set(H2,'linestyle','none','marker','.','MarkerSize',3);
set(AX(1),'XTickLabel',[0:50:350],'ylim',[-2.5 2.5],'xlim',[0 365],'ytick',[-2:1:2]); 
set(AX(2),'XTickLabel',[],'ylim',[50 300],'xlim',[0 365],'ytick',[50:50:250]);
text(5,2,'FL06','FontSize',12);
%plot([tt.t].*10, (E./(1e3))-2.5,'k')
hold on; plot([1:365],(E_FLXX(6,:)./1.25)-2.5,'k');
%%
print(gcf,'-dpng','-r500',[fn,'_NV_',fnA]);


%% set up second figure: plot N versus V with color of point the DOY
fig3=figure(3); clf; set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',[5 2 15 20]);
set(0,'DefaultAxesFontSize',10);
set(0,'DefaultTextFontSize',10);
axes_size = [0.1 0.1 0.8 0.8];
axe1 = axes('Position',[ 0.1 0.66 0.35 0.2 ],'Box','on','NextPlot','add'); 
axe2 = axes('Position',[ 0.55 0.66 0.35 0.2 ],'Box','on','NextPlot','add'); 
axe3 = axes('Position',[ 0.1 0.36 0.35 0.2 ],'Box','on','NextPlot','add');
axe4 = axes('Position',[ 0.55 0.36 0.35 0.2 ],'Box','on','NextPlot','add');
axe5 = axes('Position',[ 0.1 0.06 0.35 0.2 ],'Box','on','NextPlot','add');
axe6 = axes('Position',[ 0.55 0.06 0.35 0.2 ],'Box','on','NextPlot','add');

SS = 15; %size of scatter
days = 2:1:365;

axes(axe1)
scatter(pts_N_daymean(1,2:end), pts_V_daymean(1,2:end), SS, days);
colormap(parula); title('FL01');
ylabel('V [ m/yr ]'); 
xlim([-1 2.5]); ylim([50 350]);

axes(axe2)
scatter(pts_N_daymean(2,2:end), pts_V_daymean(2,2:end), SS, days);
colormap(parula); title('FL02');
xlim([-1 2.5]); ylim([50 350]);

axes(axe3)
scatter(pts_N_daymean(3,2:end), pts_V_daymean(3,2:end), SS, days);
colormap(parula); title('FL03');
ylabel('V [ m/yr ]');
xlim([-1 2.5]); ylim([50 350]);

axes(axe4)
scatter(pts_N_daymean(5,2:end), pts_V_daymean(4,2:end), SS, days);
colormap(parula); title('FL04');
xlim([-1 2.5]); ylim([50 350]);

axes(axe5)
scatter(pts_N_daymean(6,2:end), pts_V_daymean(5,2:end), SS, days);
colormap(parula); title('FL05');
ylabel('V [ m/yr ]'); xlabel('N [ MPa ]');
xlim([-1 2.5]); ylim([50 350]);

axes(axe6)
scatter(pts_N_daymean(7,2:end), pts_V_daymean(6,2:end), SS, days);
colormap(parula); title('FL06'); xlabel('N [ MPa ]');
xlim([-1 2.5]); ylim([50 350]);

% color bar
   cx = colorbar('peer',gca,'horizontal','position',[0.4 0.90 0.2 0.02],'XAxislocation','top');
   colormap(parula(365));
   text(-4.5,1400,'t [ d ]','FontSize',12,'Interpreter','tex');
   text(-500,180,([oo.root,oo.fn]),'FontSize',12,'Interpreter','tex');
print(gcf,'-dpng','-r500',[fn,'_NV_scatter_',fnA]);   

%% set up third figure: all N through time
fig4 = figure(4); clf; set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',[5 2 20 15]);
set(0,'DefaultAxesFontSize',10);
set(0,'DefaultTextFontSize',10);

axe1 = axes('Position',[ 0.1 0.1 0.8 0.8 ],'Box','on','NextPlot','add'); 

days2 = 1:1:365;
for i=1:1:365
    plot(days2(i).*ones(length(N_time(:,i)),1),N_time(:,i),'.'); hold on;
end
xlabel(' t [ d ]'); ylabel('N [ MPa ]');
ylim([-5 4]); xlim([0 365]);

print(gcf,'-dpng','-r500',[fn,'_Ntime_',fnA]);

%title('N across entire domain ',[oo.root,oo.fn]);



end
