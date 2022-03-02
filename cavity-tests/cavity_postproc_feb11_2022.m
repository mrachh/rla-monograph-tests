fname1 = '~/git/rla-monograph-tests/cavity-data/cavity_ik1_nk117_dk0.25_a0.3_binv2_inctype1_noise0noise_lvl0.02_data_Dirichlet.mat'; 
fname2 = '~/git/rla-monograph-tests/cavity-data/cavity_ik1_nk117_dk0.25_a0.3_binv3_inctype3_noise0noise_lvl0.02_data_Dirichlet.mat'; 
fname3 = '~/git/rla-monograph-tests/cavity-data/cavity_ik1_nk41_dk0.125_a0.3_binv3_inctype3_noise0noise_lvl0.02_data_Dirichlet.mat';
% gn
A1 = load('~/git/rla-monograph-tests/cavity-sol/cavity_ik1_nk117_dk0.25_a0.3_binv2_inctype3_noise0noise_lvl0.02_data_Dirichlet_optimtype_gn_filtertype_gauss-conv.mat');
ff1 = fname1;
% min(gn-sd)
A2 = load('~/git/rla-monograph-tests/cavity-sol/cavity_ik1_nk117_dk0.25_a0.3_binv2_inctype3_noise0noise_lvl0.02_data_Dirichlet_optimtype_min(gn,sd)_filtertype_gauss-conv.mat');
ff2 = fname1;
% sd
A3 = load('~/git/rla-monograph-tests/cavity-sol/cavity_ik1_nk117_dk0.25_a0.3_binv2_inctype3_noise0noise_lvl0.02_data_Dirichlet_optimtype_sd_filtertype_gauss-conv.mat');
ff3 = fname1;

% min(gn,sd)
A4 = load('~/git/rla-monograph-tests/cavity-sol/cavity_ik1_nk117_dk0.25_a0.3_binv3_inctype3_noise0noise_lvl0.02_data_Dirichlet_optimtype_min(gn,sd)_filtertype_gauss-conv.mat');
ff4 = fname2;

% min(gn,sd) lscaled
A5 = load('~/git/rla-monograph-tests/cavity-sol/cavity_ik1_nk117_dk0.25_a0.3_binv3_inctype3_noise0noise_lvl0.02_data_Dirichlet_optimtype_min(gn,sd)_filtertype_gauss-conv_lscaled.mat');
ff5 = fname2;

% min(gn,sd) finer freq scaling
A6 = load('~/git/rla-monograph-tests/cavity-sol/cavity_ik1_nk41_dk0.125_a0.3_binv3_inctype3_noise0noise_lvl0.02_data_Dirichlet_optimtype_min(gn,sd)_filtertype_gauss-conv_lscaled.mat');
ff6 = fname3;