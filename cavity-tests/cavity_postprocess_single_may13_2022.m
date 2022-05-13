

addpath ~/git/inverse-obstacle-scattering2d/src;

dir_data = '~/ceph/rla-monograph-tests/cavity-data/';
dir_sol = '~/ceph/rla-monograph-tests/cavity-sol/';
dir_diary = '~/ceph/rla-monograph-tests/cavity-diary/';


k0 = 1;
dkinv = 4;
dk = 1.0/dkinv;
%nk = 117;

khmax = 50;
nk = (khmax-1)*dkinv+1;

inc_type = 3;


noise_type = 0;
noise_lvl = 0.02;

% Boundary condition parameters
bc = [];
bc.type = 'Dirichlet';
bc.invtype = 'o';

ifcons = 1;
a = [0.2 0.3];
b = [2 3];
optimtype = ["gn" "sd" "min(gn,sd)" "sd-gn" "sd-min(gn,sd)"];
filtertype = ["gauss-conv" "step-length"];
eps_curv = [0.01 0.1];
ncurvmin = [0 20];
[ff,ee,dd,cc,bb,aa] = ndgrid(ncurvmin,eps_curv,optimtype,filtertype,b,a);
aa = aa(:);
bb = bb(:);
cc = cc(:);
dd = dd(:);
ee = ee(:);
ff = ff(:);

ncases = length(aa);

icases = [1:20 41:60 81:100 121:140];
icase_start = 44;
icase_end = 56;

nn = length(icases);
res_opt_all = zeros(nn,1);
res_min_all = zeros(nn,1);

aause = aa(icases);
bbuse = bb(icases);
ccuse = cc(icases);
dduse = dd(icases);
eeuse = ee(icases);
ffuse = ff(icases);

iii = 68;
icase = icases(iii);

% define geometry type
% a is a measure of the width
% b is a measure of the closing angle

disp(icase);
binv = bb(icase);
a = aa(icase);
b = pi/binv;

% optimization parameters

optim_opts = [];
opts = [];
opts.verbose=true;
optim_opts.optim_type = convertStringsToChars(dd(icase));
optim_opts.filter_type = convertStringsToChars(cc(icase));
opts.store_src_info = true;
optim_opts.maxit = 100;
opts.use_lscaled_modes = true;
optim_opts.eps_upd = 1e-3;

optim_opts.eps_curv = ee(icase);

optim_opts.n_curv_min = ff(icase);
optim_opts.sd_iter = 30;


% Data and solution directories


fname = [dir_data 'cavity_ik' num2str(k0) '_nk' int2str(nk) '_dk' ...
 num2str(dk) '_a' num2str(a) '_binv' num2str(binv)  '_inctype' ...
 int2str(inc_type) ...
 '_noise' int2str(noise_type) 'noise_lvl' num2str(noise_lvl) ... 
 '_data_' bc.type '.mat'];

fname_sol = [dir_sol 'cavity_ik' num2str(k0) '_nk' int2str(nk) '_dk' ...
 num2str(dk) '_a' num2str(a) '_binv' num2str(binv) '_inctype' ...
 int2str(inc_type) ...
 '_noise' int2str(noise_type) 'noise_lvl' num2str(noise_lvl) ... 
 '_data_' bc.type '_optimtype_' optim_opts.optim_type '_filtertype_' ...
 optim_opts.filter_type '_ifcons' int2str(ifcons) '_ncurvmin' ...
 int2str(optim_opts.n_curv_min) '_epscurv' num2str(optim_opts.eps_curv) ... 
 '_lscaled.mat'];

fname_diary = [dir_diary 'cavity_ik' num2str(k0) '_nk' int2str(nk) '_dk' ...
 num2str(dk) '_a' num2str(a) '_binv' num2str(binv) '_inctype' ...
 int2str(inc_type) ...
 '_noise' int2str(noise_type) 'noise_lvl' num2str(noise_lvl) ... 
 '_data_' bc.type '_optimtype_' optim_opts.optim_type '_filtertype_' ...
 optim_opts.filter_type '_ifcons' int2str(ifcons) '_ncurvmin' ...
 int2str(optim_opts.n_curv_min) '_epscurv' num2str(optim_opts.eps_curv) ...
 '_lscaled.mat'];

load(fname)
try
  A = load(fname_sol);
catch
  return
end
rla.post_process(A.inv_data_all,fname)
%clear A inv_tmp res_opt




