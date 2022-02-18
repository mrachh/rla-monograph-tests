fname1 = '~/ceph/rla-monograph-tests/cavity-data/cavity_ik1_nk117_dk0.25_a0.3_binv2_inctype1_noise0noise_lvl0.02_data_Dirichlet.mat'; 
fname2 = '~/ceph/rla-monograph-tests/cavity-data/cavity_ik1_nk117_dk0.25_a0.3_binv3_inctype1_noise0noise_lvl0.02_data_Dirichlet.mat'; 
% gn
A1 = load('~/ceph/rla-monograph-tests/cavity-sol/cavity_ik1_nk117_dk0.25_a0.3_binv2_inctype3_noise0noise_lvl0.02_data_Dirichlet_optimtype_gn_filtertype_gauss-conv.mat');
ff1 = fname1;
% min(gn-sd)
A2 = load('~/ceph/rla-monograph-tests/cavity-sol/cavity_ik1_nk117_dk0.25_a0.3_binv2_inctype3_noise0noise_lvl0.02_data_Dirichlet_optimtype_min(gn,sd)_filtertype_gauss-conv.mat');
ff2 = fname1;
% sd
A3 = load('~/ceph/rla-monograph-tests/cavity-sol/cavity_ik1_nk117_dk0.25_a0.3_binv2_inctype3_noise0noise_lvl0.02_data_Dirichlet_optimtype_sd_filtertype_gauss-conv.mat');
ff3 = fname1;

% min(gn,sd)
A4 = load('~/ceph/rla-monograph-tests/cavity-sol/cavity_ik1_nk117_dk0.25_a0.3_binv3_inctype3_noise0noise_lvl0.02_data_Dirichlet_optimtype_min(gn,sd)_filtertype_gauss-conv.mat');
ff4 = fname2;


