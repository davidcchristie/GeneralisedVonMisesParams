function K1K2Psi = Params_2_Offset_Pars(K1mu1K2mu2)
    K1K2Psi = [K1mu1K2mu2(:,end-3), K1mu1K2mu2(:,end-1), wrapToPi(K1mu1K2mu2(:,end)-K1mu1K2mu2(:,end-2))];
    
end