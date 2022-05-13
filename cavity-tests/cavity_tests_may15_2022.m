% This script generates the data for a star shaped domain, for a fixed
% number of sensors and incident directions where data is available for all
% sensors at each incident direction



run ~/git/inverse-obstacle-scattering2d/startup.m

warning('off');
ifgenerate = 0; % flag for generating the required data

maxNumCompThreads(6);

% define k0 (starting frequncy, dk, spacing in frequency and 
% number of frequencies (nk)

k0 = 1;
dkinv = 4;
dk = 1.0/dkinv;
%nk = 117;

khmax = 50;
nk = (khmax-1)*dkinv+1;


% setting incident waves and receptors
% inc_type = 1, 10 inc direction with 200 receptors
% inc_type = 2, 50 inc directions with 200 receptors
% inc_type = 3, 100 inc directions with 200 receptors
% inc_type = 4  2*k incident directions with 8*k receptors


% choice of noise levels
% noise_type = 0; no noise
% noise_type = 1; additive
% noise_type = 2; multiplicative

noise_type = 0;
noise_lvl = 0.02;

% Boundary condition parameters
bc = [];
bc.type = 'Dirichlet';
bc.invtype = 'o';

dir_data = '~/ceph/rla-monograph-tests/cavity-data/';
dir_sol = '~/ceph/rla-monograph-tests/cavity-sol/';
dir_diary = '~/ceph/rla-monograph-tests/cavity-diary/';

ifcons = 1;
a = [0.1 0.2 0.3];
b = [2 3 6 12];
optimtype = ["sd" "sd-gn"];
filtertype = ["gauss-conv"];
inc_type = [3 4 5];
eps_curv = [0.1];
ncurvmin = [0 20];
[gg,ff,ee,dd,cc,bb,aa] = ndgrid(inc_type,ncurvmin,eps_curv,optimtype,filtertype,b,a);
aa = aa(:);
bb = bb(:);
cc = cc(:);
dd = dd(:);
ee = ee(:);
ff = ff(:);
gg = gg(:);

ncases = length(aa);

icases = [3 6 15 18 27 30 39 42];

for iii=1:length(icases)
    icase = icases(iii); 

    % define geometry type
    % a is a measure of the width
    % b is a measure of the closing angle


    binv = bb(icase);
    a = aa(icase);
    b = pi/binv;
    inc_type = gg(icase);

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
 
     fname_sol2 = [dir_sol 'cavity_residue_ik' num2str(k0) '_nk' int2str(nk) '_dk' ...
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

    diary(fname_diary);
    fprintf('a=%d,   binv=%d\n',a,binv);
    disp(opts);
    disp(optim_opts);
 

    S = load(fname);
    u_meas = S.u_meas;


    % start inverse problem
    try
       tic, [inv_data_all,src_info_out] = rla.rla_inverse_solver(u_meas,bc,...
                          optim_opts,opts); toc;

    catch
       diary off
       fprintf('error in icase %d\n',icase);
       clear u_meas S fname_diary fname fname_sol
       clear optim_opts opts b a
       continue
    end
    diary off
    save(fname_sol,'inv_data_all','src_info_out','-v7.3');       
    inv_tmp = cell2mat(inv_data_all);
    res_opt = vertcat(inv_tmp.res_opt);
    
    save(fname_sol2,'res_opt');
    
    fprintf('\n\n\n\n\n\n\n\n');
    fprintf('done with icase: %d\n',icase);
    fprintf('\n\n\n\n\n\n\n\n');
    clear u_meas S inv_data_all src_info_out fname_diary fname fname_sol
    clear optim_opts opts b a inv_tmp inc_type res_opt fname_sol2
end
%exit                      



