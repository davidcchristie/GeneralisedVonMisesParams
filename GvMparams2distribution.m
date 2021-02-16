function [spec2D, DirBins] = GvMparams2distribution(A_K1mu1_K2mu2, DirBins)

if nargin<2
    DirBins = linspace(-pi,pi);
end

DirBins = DirBins(:)';

% ParamSet: [1 A, 2 kappa1, 3 mu1, 4 kappa2, 5 mu2]
% D(theta) = exp(kappa1*cos(theta-theta0-mu1)+kappa2*cos(2*theta-theta0-mu2)/D0
ndirbins = length(DirBins);
nspeclines = size(A_K1mu1_K2mu2, 1);

MU1bar = repmat(A_K1mu1_K2mu2(:,3), 1, ndirbins);
MU2bar = repmat(A_K1mu1_K2mu2(:,5), 1, ndirbins); 
KAP1 = repmat(A_K1mu1_K2mu2(:,2), 1, ndirbins);
KAP2 = repmat(A_K1mu1_K2mu2(:,4), 1, ndirbins);
THETA = wrapToPi(repmat(DirBins, nspeclines, 1));
A =  repmat(A_K1mu1_K2mu2(:,1), 1, ndirbins);

spec2D =exp(KAP1.*cos(wrapToPi(THETA-MU1bar))+KAP2.*cos(2*wrapToPi(THETA-MU2bar)))./A;

end
