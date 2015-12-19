function [XDEN] = func_denoise_dw2d(X)
% FUNC_DENOISE_DW2D Saved Denoising Process.
%   X: matrix of data
%   -----------------
%   XDEN: matrix of denoised data
%   cfsDEN: decomposition vector (see WAVEDEC2)
%   dimCFS: corresponding bookkeeping matrix

%  Auto-generated by Wavelet Toolbox on 16-Dec-2015 17:45:45

% Analysis parameters.
%---------------------

% Denoising parameters.
%-----------------------
% meth = 'penalhi';
% scal_OR_alfa = ;
sorh = 'h';    % Specified soft or hard thresholding

[thr,sorh,keepapp] = ddencmp('den','wv',X);

% Denoise using WDENCMP.
%----------------------

XDEN = wdencmp('gbl',X,'sym4',2,thr,sorh,keepapp);
%XDEN = wdencmp('gbl',X,'db4',2,thr,sorh,keepapp);
