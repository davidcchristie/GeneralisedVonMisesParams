function a1_b1_a2_b2 = params2moms(A_kap1_mu1_kap2_mu2, trunc)

% From a p x 5 (or p x 4) matrix of GvM parameters (A, kappa1, mu1, kappa2,
% mu2)  or (kappa1, mu1, kappa2, % mu2) returns a p x 4 matrix of moments 
% (a1, b1, a2, b2).   This is done using modified Bessel series: the
% truncation term can be set using trunc.

if nargin < 2 
    trunc = 40;
end
    
    addpath([pwd, filesep, 'functions']);
    
    OffsetParams = Params_2_Offset_Pars(A_kap1_mu1_kap2_mu2);
    
    nlines = size(OffsetParams,1);
    t0_m1_m2_n2 = zeros(nlines, 4);
    for l = 1:nlines
        t0_m1_m2_n2(l,2:end) = MomsAndJacobianFlexi(OffsetParams(l,:), trunc);
        locsT0 = MomsAndJacobianFlexi(OffsetParams(l,:),trunc, true);
        mu1 = A_kap1_mu1_kap2_mu2(l,end-2);
        t0_m1_m2_n2(l,1) =  mu1-locsT0(1,1);
    end
    
    a1_b1_a2_b2 = Offset_Moms_2_Moms(t0_m1_m2_n2);
end