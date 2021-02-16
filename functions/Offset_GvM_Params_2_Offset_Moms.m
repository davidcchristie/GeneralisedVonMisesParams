function OffsetMoms = Offset_GvM_Params_2_Offset_Moms(K1K2Psi, trunc)
%     Takes a vector or array of vectors containing kappa1, kappa2 and Psi (the two
%     concentration parameters) and the relative location parameter mu2-mu1
%     and calculates the offset moments
            
    if nargin<2
        trunc = 40;
    end
    
    OffsetMoms = MomsAndJacobianFlexi(K1K2Psi, trunc);
    
        
end