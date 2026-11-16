function [feat, names] = buildFeatureStack_LogRatio(Aclean, cfg)
% buildFeatureStack_LogRatio  Construct mode-specific Phase 5 feature stack.

B  = Aclean(:,:,cfg.BANDS.BLUE);
G1 = Aclean(:,:,cfg.BANDS.GREEN1);
G2 = Aclean(:,:,cfg.BANDS.GREEN2);
R  = Aclean(:,:,cfg.BANDS.RED);

B_over_G1 = B ./ max(G1, cfg.LR.EPS);
G1_over_G2 = G1 ./ max(G2, cfg.LR.EPS);
B_over_R = B ./ max(R, cfg.LR.EPS);
DII_BR = log(max(B, cfg.LR.EPS)) - cfg.LR.K_BR * log(max(R, cfg.LR.EPS));
DII_GR = log(max(G1, cfg.LR.EPS)) - cfg.LR.K_GR * log(max(R, cfg.LR.EPS));
VisMean = (B + G1 + R) ./ 3;
RGBstd = std(cat(3, B, G1, R), 0, 3);
Whiteness = 1 - RGBstd ./ max(VisMean, cfg.LR.EPS);
localKernel = ones(3,3,'single') ./ 9;
muLocal = conv2(single(DII_GR), localKernel, 'same');
mu2Local = conv2(single(DII_GR).^2, localKernel, 'same');
LocalVar_DII_GR = max(mu2Local - muLocal.^2, 0);

feat = cat(3, B_over_G1, G1_over_G2, B_over_R, DII_BR, DII_GR);
names = {'B_over_G1','G1_over_G2','B_over_R','DII_BR','DII_GR'};

switch upper(string(cfg.FEATURE_MODE))
    case "A3"
        % Production default: ratios + DII only.

    case "A5B"
        % Baseline operational stack: ratios + DII only.

    case "A4"
        feat = cat(3, feat, VisMean, Whiteness, LocalVar_DII_GR);
        names = [names, {'VisMean','Whiteness','LocalVar_DII_GR'}];

    case "C1"
        blue = Aclean(:,:,cfg.BANDS.BLUE);
        nir  = Aclean(:,:,cfg.BANDS.NIR);
        [blue_contrast, blue_entropy] = computeOldStyleGLCM(blue, cfg.LR.GLCM_WINDOW);
        [nir_contrast,  nir_entropy]  = computeOldStyleGLCM(nir,  cfg.LR.GLCM_WINDOW);
        feat = cat(3, feat, blue_contrast, blue_entropy, nir_contrast, nir_entropy);
        names = [names, {'Blue_GLCM_Contrast','Blue_GLCM_Entropy', ...
                         'NIR_GLCM_Contrast','NIR_GLCM_Entropy'}];

    otherwise
        error('Unsupported FEATURE_MODE: %s', cfg.FEATURE_MODE);
end
end
