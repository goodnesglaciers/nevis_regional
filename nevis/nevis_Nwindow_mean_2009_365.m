function nevis_Nwindow_mean_2009_365(fn)
%% nevis N and Q for entire year save
%% LAS June 2017

tis = [1:1:365];

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


%% load timesteps over the year
for i_s=1:1:365    
    disp(['DOY ',num2str(i_s),' / ',num2str(length(tis)),' ...']); 
    load([fn,'/',int2four(tis(i_s))]); 
    %% extract new variables
    if isfield(vv,'nbdy'), gg = nevis_label(gg,vv.nbdy); aa.phi = aa.phi_a(gg.nbdy); end % remask
    aa = nevis_inputs(vv.t,aa,pp,gg,oo);
    oo.evaluate_variables = 1; oo.evaluate_residual = 0; oo.evaluate_jacobian = 0; 
    [vv2] = nevis_backbone(inf,vv,vv,aa,pp,gg,oo); % expand solution variables
    vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo); % calculate node discharge
    nevis_unpack(aa,vv2);
    clear tt vv vv2

    % effective pressure and Q 
    zz = (ps.phi/10^6)*reshape(phi_0-phi,gg.nI,gg.nJ); 
    zzQ = ps.qs*reshape(qs+qe+qQ, gg.nI,gg.nJ);
    
    zzall_nevis(1,i_s).N = zz;
    zzall_nevis(1,i_s).Q = zzQ;
    zzall_nevis(1,i_s).t = tis(i_s);

end
%% save zz   
   save(num2str([fn,'_zzall_nevis.mat']),'zzall_nevis','-v7.3')
end

% RUN COHERENCE SCRIPTS FROM HERE
%year = 2009;
%[COH2, covar] = nevis_coherence_2009(fn, tis, zz_nevis); 
%[COH2, covar] = nevis_coherence_V2009(fn, tis, zz_nevis); 

