function nevis_Nwindow_mean_2010(fn,tis)
%% nevis surface inputs (E) 8 panel plot to go with Joughin 2013 velocities
%% LAS August 2016
DOI_2009 = [163, 174, 185, 196, 207, 218, 229, 240]; % Date of Terra-SARX windows 2009
DOI_2010 = [150, 161, 172, 183, 194, 211, 233, 249]; % Date of Terra-SARX windows 2010

% DOI = TIS
span = 5; % 11-day interval

%% load initial timestep
if exist([fn,'/0000']), load([fn,'/0000']); else load(fn); end
load([fn,'/0001']);
if isfield(oo,'reversey'), reversey = oo.reversey; else reversey = 0; end

%% contour levels
% z_conts = -1000:50:5000;
z_conts = -1000:10:5000;
b_conts = z_conts;
s_conts = z_conts;
p_conts = z_conts/10;

%% axes limits
axx = (ps.x/10^3)*[gg.xl-0.01 gg.xr gg.yb-0.01 gg.yt];

%% extract initial variables
if isfield(vv,'nbdy'), gg = nevis_label(gg,vv.nbdy); aa.phi = aa.phi_a(gg.nbdy); end % re-mask
oo.evaluate_variables = 1; oo.evaluate_residual = 0; oo.evaluate_jacobian = 0; [vv2] = nevis_backbone(inf,vv,vv,aa,pp,gg,oo); % expand solution variables
vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo); % calculate node discharge
nevis_unpack(aa,gg);

%get rid of points outside domain
nx(gg.nout) = NaN;
ex(gg.eout) = NaN;
fx(gg.fout) = NaN;
cx(gg.cout) = NaN;
    
%boundary curve
x_out = gg.nx(gg.n1); y_out = gg.ny(gg.n1);
if ~isempty(x_out), tmp = nevis_orderboundary(x_out,y_out); x_out = x_out(tmp); y_out = y_out(tmp); end % reorder to follow boundary

%% set up figure
figure(2); clf; set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',1.5.*[1 1 20 8]);
set(0,'DefaultAxesFontSize',10);
set(0,'DefaultTextFontSize',10);
axes_size = [0.1 0.1 0.4 0.4];

axe1(1) = axes('Position',[0.08 0.45 0.2 0.4],'Box','on','NextPlot','add'); 
axe1(2) = axes('Position',[0.30 0.45 0.2 0.4],'Box','on','NextPlot','add'); 
axe1(3) = axes('Position',[0.52 0.45 0.2 0.4],'Box','on','NextPlot','add'); 
axe1(4) = axes('Position',[0.74 0.45 0.2 0.4],'Box','on','NextPlot','add'); 
axe1(5) = axes('Position',[0.08 0.05 0.2 0.4],'Box','on','NextPlot','add'); 
axe1(6) = axes('Position',[0.30 0.05 0.2 0.4],'Box','on','NextPlot','add'); 
axe1(7) = axes('Position',[0.52 0.05 0.2 0.4],'Box','on','NextPlot','add'); 
axe1(8) = axes('Position',[0.74 0.05 0.2 0.4],'Box','on','NextPlot','add');

%adjust axes positions to account for aspect ratio of figure
aspect = (axx(4)-axx(3))/(axx(2)-axx(1)); % aspect of axes limits
pos = axes_size; tmp = get(gcf,'position'); sca = tmp(4)/tmp(3);
if aspect>(pos(4)/pos(3))*sca, 
    pos = pos + [pos(4)*sca*(1-1/aspect)/2 0 -pos(3)+pos(4)*sca*(1/aspect) 0]; % constrained by y
else
    pos = pos + [0 pos(3)/sca*(1-aspect)/2 0 -pos(4)+pos(3)/sca*aspect]; % constrained by x
end 
axes_size = pos;

%% load Joughin 2013 The Cryosphere TerraSAR-X footprint
load('nevis/terravel.mat')
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

%% determine wintertime average of 30 days values of N
for i_s=10:40
    load([fn,'/',int2four(i_s)]); 
    %% extract new variables
    if isfield(vv,'nbdy'), gg = nevis_label(gg,vv.nbdy); aa.phi = aa.phi_a(gg.nbdy); end % remask
    aa = nevis_inputs(vv.t,aa,pp,gg,oo);
    oo.evaluate_variables = 1; oo.evaluate_residual = 0; oo.evaluate_jacobian = 0; [vv2] = nevis_backbone(inf,vv,vv,aa,pp,gg,oo); % expand solution variables
    vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo); % calculate node discharge
    nevis_unpack(aa,vv2);
    % save variable N
    N(:,i_s) = (ps.phi/10^6)*(phi_0-phi); 
    clear tt vv vv2
end
Nwinter = nanmean(N,2);

%% plot effective pressure on central date of the 11-day interval
for i_t = 1:length(tis)

disp(['Frame ',num2str(i_t),' / ',num2str(length(tis)),' ...']);    

%% load timesteps +/- 5 days from tis
span = (-5:1:5);
for i_s=1:length(span)    
    load([fn,'/',int2four(tis(i_t)-span(i_s))]); 
    %% extract new variables
    if isfield(vv,'nbdy'), gg = nevis_label(gg,vv.nbdy); aa.phi = aa.phi_a(gg.nbdy); end % remask
    aa = nevis_inputs(vv.t,aa,pp,gg,oo);
    oo.evaluate_variables = 1; oo.evaluate_residual = 0; oo.evaluate_jacobian = 0; [vv2] = nevis_backbone(inf,vv,vv,aa,pp,gg,oo); % expand solution variables
    vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo); % calculate node discharge
    nevis_unpack(aa,vv2);
    % save plotting variable phi
    phi0_i_s(:,i_s) = phi_0; 
    phi_i_s(:,i_s) = phi;
    clear tt vv vv2
end
% average plotting variable over 11 days
phi0mean = nanmean(phi0_i_s,2);
phimean = nanmean(phi_i_s,2);

%% plot variable averages
axes(axe1(i_t))   
    % effective pressure
    xx = (ps.x/10^3)*nx; yy = (ps.x/10^3)*ny;  
    zz = (ps.phi/10^6)*reshape(phi0mean-phimean,gg.nI,gg.nJ); 
    cax = [-1 1]; % [-3 3]; 
    
    % save zz variable at time i_t for covariance stats with terrasarX
    % velocities later down
    zz_nevis(1,i_t).N = zz;
    zz_nevis(1,i_t).t = tis(i_t);
    
%     zz = max(cax(1),min(zz,cax(2))); % cap so as to avoid going off color axis     
%     hand = pcolor(xx,yy,zz); set(hand,'linestyle','none'); % shading interp;
%     caxis(cax); 
%     axis image;
%     axis(axx);
%     
%     load cmapbluered; colormap(cmap(end:-1:1,:));
%     
% %     % add moulins
% %     if ~isfield(pp,'ni_m'), pp.ni_m = []; end
% %     mscale = 100;  %amount to scale moulin input by to make size of moulin dot
% %     x = (ps.x/10^3)*nx(pp.ni_m);
% %     y = (ps.x/10^3)*ny(pp.ni_m);
% %     hold on;
% %     for i_m = 1:length(pp.ni_m),
% %         if E(pp.ni_m(i_m))>0,
% %             plot(x(i_m),y(i_m),'ko','Markersize',4+E(pp.ni_m(i_m))/mscale,'MarkerFaceColor',1*[1 1 1]); % mark moulins   
% %         else
% %             plot(x(i_m),y(i_m),'ko','Markersize',4,'MarkerFaceColor',0.8*[1 1 1]); % mark moulins  
% %         end
% %     end
%     
%     % add outlets
%     plot((ps.x/10^3)*nx(gg.nbdy),(ps.x/10^3)*ny(gg.nbdy),'ko','markersize',4,'markerfacecolor','y');
%     
%     % add boundary
%     plot((ps.x/10^3)*x_out,(ps.x/10^3)*y_out,'k','linewidth',1.5); 
%     
%     % add Joughin 2013 footprint
%     patch(footprint(:,1)./1000, footprint(:,2)./1000,'k','FaceColor','none'); 
%     
%     % add flowline
%     scatter(flxx_plot(:,1)./10^3, flxx_plot(:,2)./10^3, 30, 'k^','filled');
%     
%     %% formatting   
%     if reversey, set(gca,'YDir','reverse'); end
% %     text(-0.08,0.5,'y [ km ]','HorizontalAlignment','center','units','normalized','rotation',90);
% %     text(0.5,-0.08/aspect,'x [ km ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized');
%     
%     text(-80,33,['t = ',num2str(tis(i_t)),' +/- 5.5']);
%     set(gcf,'Color','w'); set(gca,'yticklabel',[]);
%     drawnow; shg;
end

% % color bar
%    cx = colorbar('peer',gca,'horizontal','position',[0.4 0.9 0.2 0.02],'XAxislocation','top');
%    load cmapbluered; colormap(cmap(end:-1:1,:)); 
%    text(-410,175,'N [ MPa ]','FontSize',12,'Interpreter','tex');
%    text(-500,180,([oo.root,oo.fn]),'FontSize',12,'Interpreter','tex');
%    
%    print(gcf,'-dpng','-r500',[fn,'_Nwindow_mean']);

%% save zz   
   zz_nevis(1,9).N = reshape(Nwinter, size(xx,1), size(xx,2));
   zz_nevis(1,9).t = 365;
   save(num2str([fn,'_zz_nevis.mat']),'zz_nevis','-v7.3')
   
% RUN COHERENCE SCRIPT FROM HERE
year = 2010;
[COH2, covar] = nevis_coherence_2010(fn, tis, zz_nevis); 
%[COH2, covar] = nevis_coherence_V2010(fn, tis, zz_nevis); 
   
%% PLOT TERRASARX VELOCITIES   
%    
% % load TerraSARX velocities from tis days 2010
% DOI_2010 = [150, 161, 172, 183, 194, 211, 233, 249]; % Date of Terra-SARX windows 2010
% 
% figure(3); clf; set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',1.5.*[1 1 20 8]);
% set(0,'DefaultAxesFontSize',10);
% set(0,'DefaultTextFontSize',10);
% axes_size = [0.1 0.1 0.4 0.4];
% 
% axe1(9) = axes('Position',[0.08 0.45 0.2 0.4],'Box','on','NextPlot','add'); 
% axe1(10) = axes('Position',[0.30 0.45 0.2 0.4],'Box','on','NextPlot','add'); 
% axe1(11) = axes('Position',[0.52 0.45 0.2 0.4],'Box','on','NextPlot','add'); 
% axe1(12) = axes('Position',[0.74 0.45 0.2 0.4],'Box','on','NextPlot','add'); 
% axe1(13) = axes('Position',[0.08 0.05 0.2 0.4],'Box','on','NextPlot','add'); 
% axe1(14) = axes('Position',[0.30 0.05 0.2 0.4],'Box','on','NextPlot','add'); 
% axe1(15) = axes('Position',[0.52 0.05 0.2 0.4],'Box','on','NextPlot','add'); 
% axe1(16) = axes('Position',[0.74 0.05 0.2 0.4],'Box','on','NextPlot','add');
% 
% axes(axe1(9))
%     hand = pcolor(terravel.polar.xx_moulin./10^3,terravel.polar.yy_moulin./10^3,...
%         ((terravel.x16317.vel-terravel.winter.vel)./terravel.winter.vel).*100); 
%     set(hand,'linestyle','none'); % shading interp;
%     cax = [-25 125]; caxis(cax); axis image; axis(axx);   
%     % add boundary
%     plot((ps.x/10^3)*x_out,(ps.x/10^3)*y_out,'k','linewidth',1.5);     
%    % add flowline
%     scatter(flxx_plot(:,1)./10^3, flxx_plot(:,2)./10^3, 30, 'k^','filled');
%     % formatting   
%     if reversey, set(gca,'YDir','reverse'); end
%     text(-80,33,['t = 150 +/- 5.5']);
%     set(gcf,'Color','w'); set(gca,'yticklabel',[]);
% 
% axes(axe1(10))
%     hand = pcolor(terravel.polar.xx_moulin./10^3,terravel.polar.yy_moulin./10^3,...
%           ((terravel.x16484.vel-terravel.winter.vel)./terravel.winter.vel).*100); 
%     set(hand,'linestyle','none'); % shading interp;
%     cax = [-25 125]; caxis(cax); axis image; axis(axx);  
%     % add boundary
%     plot((ps.x/10^3)*x_out,(ps.x/10^3)*y_out,'k','linewidth',1.5);     
%     % add flowline
%     scatter(flxx_plot(:,1)./10^3, flxx_plot(:,2)./10^3, 30, 'k^','filled');
%     % formatting   
%     if reversey, set(gca,'YDir','reverse'); end
%     text(-80,33,['t = 161 +/- 5.5']);
%     set(gcf,'Color','w'); set(gca,'yticklabel',[]);
% 
% axes(axe1(11))
%     hand = pcolor(terravel.polar.xx_moulin./10^3,terravel.polar.yy_moulin./10^3,...
%     ((terravel.x16651.vel-terravel.winter.vel)./terravel.winter.vel).*100);        
%     set(hand,'linestyle','none'); % shading interp;
%     cax = [-25 125]; caxis(cax); axis image; axis(axx);   
%     % add boundary
%     plot((ps.x/10^3)*x_out,(ps.x/10^3)*y_out,'k','linewidth',1.5);     
%     % add flowline
%     scatter(flxx_plot(:,1)./10^3, flxx_plot(:,2)./10^3, 30, 'k^','filled');
%     % formatting   
%     if reversey, set(gca,'YDir','reverse'); end
%     text(-80,33,['t = 172 +/- 5.5']);
%     set(gcf,'Color','w'); set(gca,'yticklabel',[]);    
% 
% axes(axe1(12))
%     hand = pcolor(terravel.polar.xx_moulin./10^3,terravel.polar.yy_moulin./10^3,...
%    ((terravel.x16818.vel-terravel.winter.vel)./terravel.winter.vel).*100);        
%     set(hand,'linestyle','none'); % shading interp;
%     cax = [-25 125]; caxis(cax); axis image; axis(axx);   
%     % add boundary
%     plot((ps.x/10^3)*x_out,(ps.x/10^3)*y_out,'k','linewidth',1.5);     
%     % add flowline
%     scatter(flxx_plot(:,1)./10^3, flxx_plot(:,2)./10^3, 30, 'k^','filled');
%     % formatting   
%     if reversey, set(gca,'YDir','reverse'); end
%     text(-80,33,['t = 183 +/- 5.5']);
%     set(gcf,'Color','w'); set(gca,'yticklabel',[]);    
% 
% axes(axe1(13))
%     hand = pcolor(terravel.polar.xx_moulin./10^3,terravel.polar.yy_moulin./10^3,...
%     ((terravel.x16985.vel-terravel.winter.vel)./terravel.winter.vel).*100);        
%     set(hand,'linestyle','none'); % shading interp;
%     cax = [-25 125]; caxis(cax); axis image; axis(axx);  
%     % add boundary
%     plot((ps.x/10^3)*x_out,(ps.x/10^3)*y_out,'k','linewidth',1.5);     
%     % add flowline
%     scatter(flxx_plot(:,1)./10^3, flxx_plot(:,2)./10^3, 30, 'k^','filled');
%     % formatting   
%     if reversey, set(gca,'YDir','reverse'); end
%     text(-80,33,['t = 194 +/- 5.5']);
%     set(gcf,'Color','w'); set(gca,'yticklabel',[]);  
% 
% axes(axe1(14))
%     hand = pcolor(terravel.polar.xx_moulin./10^3,terravel.polar.yy_moulin./10^3,...
%     ((terravel.x17152.vel-terravel.winter.vel)./terravel.winter.vel).*100);      
%     set(hand,'linestyle','none'); % shading interp;
%     cax = [-25 125]; caxis(cax); axis image; axis(axx); 
%     % add boundary
%     plot((ps.x/10^3)*x_out,(ps.x/10^3)*y_out,'k','linewidth',1.5);     
%     % add flowline
%     scatter(flxx_plot(:,1)./10^3, flxx_plot(:,2)./10^3, 30, 'k^','filled');
%     % formatting   
%     if reversey, set(gca,'YDir','reverse'); end
%     text(-80,33,['t = 211 +/- 5.5']);
%     set(gcf,'Color','w'); set(gca,'yticklabel',[]);      
% 
% axes(axe1(15))
%     hand = pcolor(terravel.polar.xx_moulin./10^3,terravel.polar.yy_moulin./10^3,...
%     ((terravel.x17486.vel-terravel.winter.vel)./terravel.winter.vel).*100);        
%     set(hand,'linestyle','none'); % shading interp;
%     cax = [-25 125]; caxis(cax); axis image; axis(axx);  
%     % add boundary
%     plot((ps.x/10^3)*x_out,(ps.x/10^3)*y_out,'k','linewidth',1.5);     
%     % add flowline
%     scatter(flxx_plot(:,1)./10^3, flxx_plot(:,2)./10^3, 30, 'k^','filled');
%     % formatting   
%     if reversey, set(gca,'YDir','reverse'); end
%     text(-80,33,['t = 233 +/- 5.5']);
%     set(gcf,'Color','w'); set(gca,'yticklabel',[]);  
%     
% axes(axe1(16))
%     hand = pcolor(terravel.polar.xx_moulin./10^3,terravel.polar.yy_moulin./10^3,...
%     ((terravel.x17820.vel-terravel.winter.vel)./terravel.winter.vel).*100);      
%     set(hand,'linestyle','none'); % shading interp;
%     cax = [-25 125]; caxis(cax); axis image; axis(axx);  
%     % add boundary
%     plot((ps.x/10^3)*x_out,(ps.x/10^3)*y_out,'k','linewidth',1.5);     
%     % add flowline
%     scatter(flxx_plot(:,1)./10^3, flxx_plot(:,2)./10^3, 30, 'k^','filled');
%     % formatting   
%     if reversey, set(gca,'YDir','reverse'); end
%     text(-80,33,['t = 249 +/- 5.5']);
%     set(gcf,'Color','w'); set(gca,'yticklabel',[]);  
%     
% % color bar
%    cx = colorbar('peer',gca,'horizontal','position',[0.4 0.9 0.2 0.02],'XAxislocation','top');
%    text(-410,175,'Velocity Difference [ % ]','FontSize',12,'Interpreter','tex');
%    text(-500,180,([oo.root,oo.fn]),'FontSize',12,'Interpreter','tex');
%    
%    print(gcf,'-dpng','-r500',[fn,'_Nwindow_vel']);
   
