function [contrast, entropy_] = computeOldStyleGLCM(img, win)
% computeOldStyleGLCM  Match the legacy texture surrogate from Preprocess_DII
img_norm = mat2gray(img);
img_u8   = im2uint8(img_norm);
se       = true(win, win);
contrast = single(stdfilt(img_norm, se).^2);
entropy_ = single(entropyfilt(img_u8, se));
end
