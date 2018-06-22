function [ps,pp,sc] = nevis_nondimension(pd,ps,oo)
% Nondimensionalise parameters in pd with scales in ps
% Inputs 
%   pd struct of dimensional parameters
%   ps struct of scales to use for non-dimensionalization [optional]
%   oo struct of options [optional]
% Outputs
%   ps struct of scales
%   pp struct of scaled parameters
%   sc scales for plot [optional; for compatability with older versions]
%
% 20 July 2014 : taken from hydro_nondimension
         
if nargin<2 || isempty(ps), ps = struct; end 
if nargin<3, oo = struct; end

%% Default scales
if ~isfield(ps,'u_b'), ps.u_b = pd.u_b; end
if ~isfield(ps,'x'), ps.x = 10*10^3; end
if ~isfield(ps,'z'), ps.z = 1*10^3; end
if ~isfield(ps,'phi'), ps.phi = pd.rho_i*pd.g*ps.z; end
if ~isfield(ps,'l_c'), ps.l_c = pd.l_c; end
if ~isfield(ps,'m'), ps.m = 10*10^(-3)/pd.td; end
if ~isfield(ps,'h'), ps.h = 0.1; end
if ~isfield(ps,'sigma'), ps.sigma = 1e-5; end

if ~isfield(ps,'Psi'), ps.Psi = ps.phi/ps.x; end
if ~isfield(ps,'hs'), ps.hs = ps.h; end
if ~isfield(ps,'hv'), ps.hv = ps.h; end
if ~isfield(ps,'he'), ps.he = ps.h; end
if ~isfield(ps,'qs'), ps.qs = ps.m*ps.x; end
if ~isfield(ps,'qe'), ps.qe = ps.qs; end
if ~isfield(ps,'Q'), ps.Q = ps.qs*100; end
if ~isfield(ps,'S'), ps.S = (ps.Q/pd.k_c/ps.Psi^pd.beta_c).^(1/pd.alpha_c); end
if ~isfield(ps,'t'), ps.t = ps.h*ps.x/ps.qs; end
if ~isfield(ps,'Xi'), ps.Xi = ps.Q*ps.Psi; end

%% Dimensionless parameters [ many of these can be chosen to be 1 by suitable choice of scales ]
pp.c1 = ps.hs/ps.h;
pp.c2 = ps.he/ps.h;
pp.c3 = ps.hv/ps.h;
pp.c4 = ps.qs*ps.t/ps.h/ps.x;
pp.c5 = ps.qe*ps.t/ps.h/ps.x;
pp.c6 = ps.m*ps.t/ps.h;
pp.c7 = ps.m*ps.t/ps.h;
pp.c8 = ps.S/ps.h/ps.x;
pp.c9 = ps.Q*ps.t/ps.h/ps.x^2;
pp.c10 = ps.m*ps.t/ps.h;
pp.c11 = ps.Xi*ps.t/ps.h/ps.x/pd.rho_w/pd.L;
pp.c12 = 1;
pp.c13 = ps.Xi*ps.t/ps.S/pd.rho_i/pd.L;
pp.c14 = pd.K_c*ps.t*ps.phi^pd.n_Glen;
pp.c15 = ps.hs/ps.h;
pp.c16 = ps.u_b*pd.h_r/pd.l_r*ps.t/ps.h;
pp.c17 = ps.hs/max(pd.h_r,eps);
pp.c18 = pd.K_s*ps.t*ps.phi^pd.n_Glen;
pp.c19 = ps.Q*ps.phi/ps.Xi/ps.x;
pp.c20 = ps.qs*ps.l_c*ps.phi/ps.Xi/ps.x;
pp.c21 = pd.k_s*ps.hs^pd.alpha_s*ps.Psi^pd.beta_s/ps.qs;
pp.c22 = pd.k_e*ps.he^pd.alpha_e*ps.Psi^pd.beta_e/ps.qe;
pp.c23 = pd.k_c*ps.S^pd.alpha_c*ps.Psi^pd.beta_c/ps.Q;
pp.c24 = pd.c_e_power/ps.he;
pp.c25 = ps.sigma/pd.rho_w/pd.g*ps.phi/ps.hv;
pp.c26 = pd.rho_w/pd.rho_i*ps.m*ps.t/ps.h;
pp.c27 = pd.rho_w/pd.rho_i;
pp.c28 = pd.c_e_log/ps.he; 
pp.c29 = pd.c_e_reg1/ps.he; 
pp.c30 = pd.N_reg1/ps.phi; 
pp.c31 = pd.c_e_reg2*ps.phi/ps.he; 
pp.c32 = pd.N_reg2/ps.phi; 
pp.c33 = pd.V_m_reg/max(pd.A_m,eps)*pd.rho_w*pd.g/ps.phi; 
pp.c34 = pd.p_w_reg/ps.phi; 
pp.c35 = ps.u_b*pd.h_rc*ps.t/ps.S;
pp.c36 = ps.S/max(pd.S_rc,eps);
pp.c37 = pd.c_e_reg3/ps.he; 
pp.c38 = pd.N_reg3/ps.phi; 
pp.c39 = pd.sigma_log/ps.he;
pp.c40 = pd.N_sigma/ps.phi;
pp.c41 = pd.rho_w*pd.gamma_cc*pd.c;

pp.n_Glen = pd.n_Glen;
pp.alpha_s = pd.alpha_s;
pp.beta_s = pd.beta_s;
pp.alpha_c = pd.alpha_c;
pp.beta_c = pd.beta_c;
pp.alpha_e = pd.alpha_e;
pp.beta_e = pd.beta_e;
pp.gamma_e = pd.gamma_e;
pp.r = pd.rho_w/pd.rho_i;
pp.sigma = pd.sigma/ps.sigma;
pp.A_m = pd.A_m/ps.x^2/ps.sigma;  
pp.phi_s = pd.phi_s/ps.phi;
pp.td = pd.td/ps.t;
pp.ty = pd.ty/ps.t;

pp.Psi_reg = pd.Psi_reg/ps.Psi; 
pp.N_reg = pd.N_reg/ps.phi;
%pp.N_reg1 = pd.N_reg1/ps.phi;
pp.p_w_reg = pd.p_w_reg/ps.phi;
pp.p_a_reg = pd.p_a_reg/ps.phi;
pp.V_m_reg = pd.V_m_reg/ps.x^2/ps.sigma*pd.rho_w*pd.g/ps.phi;

pp.E_lapse = pd.E_lapse/ps.m*ps.z;
pp.E_amp = pd.E_amp; 
pp.melt = pd.melt/ps.m;
pp.meltinterior = pd.meltinterior/ps.m;
pp.u_b = pd.u_b/ps.u_b;
pp.l_c = pd.l_c/ps.l_c;

%% Scalings for plots 
sc.x = ps.x/1000;   % km
sc.z = ps.z;        % m
sc.S = ps.S;        % m
sc.h = ps.h*100;    % cm
sc.phi = ps.phi/10^6; % MPa
sc.m = ps.m*10^3*pd.td;   % mm / d
sc.qs = ps.qs;      % m^2 / s
sc.qe = ps.qe;      % m^2 / s
sc.t = ps.t/pd.td;  % d
sc.Q = ps.Q;        % m^3 / s

end

