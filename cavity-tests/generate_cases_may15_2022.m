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



