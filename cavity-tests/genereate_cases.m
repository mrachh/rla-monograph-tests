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



