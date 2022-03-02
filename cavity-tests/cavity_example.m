% This script generates the data for a star shaped domain, for a fixed
% number of sensors and incident directions where data is available for all
% sensors at each incident direction



run ~/git/inverse-obstacle-scattering2d/startup.m

warning('off');
ifgenerate = 1; % flag for generating the required data

% define geometry type
% a is a measure of the width
% b is a measure of the closing angle

maxNumCompThreads(4);


binv = 3;
a = 0.3;
b = pi/binv;

% define k0 (starting frequncy, dk, spacing in frequency and 
% number of frequencies (nk)

k0 = 1;
dkinv = 4;
dk = 1.0/dkinv;
%nk = 117;

khmax = 5;
nk = (khmax-1)*dkinv+1;


% setting incident waves and receptors
% inc_type = 1, 10 inc direction with 200 receptors
% inc_type = 2, 50 inc directions with 200 receptors
% inc_type = 3, 100 inc directions with 200 receptors
% inc_type = 4  2*k incident directions with 8*k receptors

inc_type = 1;


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


% optimization parameters

optim_opts = [];
opts = [];
opts.verbose=true;
optim_opts.optim_type = 'min(gn,sd)';
optim_opts.filter_type = 'gauss-conv';
opts.store_src_info = true;
optim_opts.eps_curv = 1e-2;
optim_opts.maxit = 200;
opts.use_lscaled_modes = true;
optim_opts.n_curv_min = 30;


% Data and solution directories
dir_data = '../cavity-data/';
dir_sol = '../cavity-sol/';
dir_diary = '../cavity-diary/';


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
     optim_opts.filter_type '_lscaled.mat'];

fname_diary = [dir_diary 'cavity_ik' num2str(k0) '_nk' int2str(nk) '_dk' ...
     num2str(dk) '_a' num2str(a) '_binv' num2str(binv) '_inctype' ...
     int2str(inc_type) ...
     '_noise' int2str(noise_type) 'noise_lvl' num2str(noise_lvl) ... 
     '_data_' bc.type '_optimtype_' optim_opts.optim_type '_filtertype_' ...
     optim_opts.filter_type '_lscaled.mat'];

diary(fname_diary);
 





src0 = [0.01;-1.1];
opts.test_analytic = true;
opts.src_in = src0;
opts.verbose=true;

% Set of frequencies (k_{i})
kh = 1:dk:(1+(nk-1)*dk);



if (ifgenerate)
     n  = 300;

    src_info = geometries.larrycup(a,b,n);
    L = src_info.L;

    save(fname,'src_info');
    u_meas = cell(nk,1);

    nppw = 20;

    for ik=1:nk
        if(inc_type == 1)
            n_dir = 10;
            n_tgt = 200;
        end
        if(inc_type == 2)
            n_dir = 50;
            n_tgt = 200;
        end
        if(inc_type == 3)
            n_dir = 100;
            n_tgt = 200;
        end
        if(inc_type == 4)
            n_dir = ceil(2*kh(ik));
            n_tgt = ceil(8*kh(ik));
        end

        % set target locations
        %receptors (r_{\ell})
        r_tgt = 10;
        t_tgt = 0:2*pi/n_tgt:2*pi-2*pi/n_tgt;

        % Incident directions (d_{j})
        t_dir = 0:2*pi/n_dir:2*pi-2*pi/n_dir;

        [t_tgt_grid,t_dir_grid] = meshgrid(t_tgt,t_dir);
        t_tgt_grid = t_tgt_grid(:);
        t_dir_grid = t_dir_grid(:);
        xtgt = r_tgt*cos(t_tgt_grid);
        ytgt = r_tgt*sin(t_tgt_grid);
        tgt   = [ xtgt'; ytgt'];


        sensor_info = [];
        sensor_info.tgt = tgt;
        sensor_info.t_dir = t_dir_grid;


        n = ceil(nppw*L*abs(kh(ik))/2/pi);
        n = max(n,300);
        src_info = geometries.larrycup(a,b,n);

       [mats,erra] = rla.get_fw_mats(kh(ik),src_info,bc,sensor_info,opts);
       fields = rla.compute_fields(kh(ik),src_info,mats,sensor_info,bc,opts);

       u_meas0 = [];
       u_meas0.kh = kh(ik);
       u_meas0.uscat_tgt = fields.uscat_tgt;
       u_meas0.tgt = sensor_info.tgt;
       u_meas0.t_dir = sensor_info.t_dir;
       u_meas0.err_est = erra;
       u_meas{ik} = u_meas0;
       fprintf('kh = %d    err= %d\n',kh(ik),u_meas0.err_est);
    end


    save(fname,'u_meas','-append');
else
    S = load(fname);
    u_meas = S.u_meas;
end



% start inverse problem
[inv_data_all,src_info_out] = rla.rla_inverse_solver(u_meas,bc,...
                          optim_opts,opts);

diary off
save(fname_sol,'inv_data_all','src_info_out','-v7.3');                      
                      



