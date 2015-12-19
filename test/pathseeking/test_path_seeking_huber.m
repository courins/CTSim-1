load 'temp.mat';

% Load simulation parameters and datas
[phan, map] = loadMaterialsDensityPhantom(p);

[ geom ] = loadProjectionGeometryCT( p );

spectrum = loadSpectraCT(p, geom, 2e5);

%Compute ground trut
[imgGtAtt, imgGtHu ] = computeGroundTruth(phan, spectrum);
%figure; imdisp( imgGtAtt, [0.2 0.24] );


%% compute sinogram

[sinoRaw, tubeCurrentProfile] = simulateCTRawData( phan, geom, spectrum, sinosDir, 1, 1, 2 );

sinoAtt = processCTRawData( sinoRaw, spectrum, tubeCurrentProfile);

weights = computeWeightsPwls( sinoRaw, 0, spectrum.electronicNoise );

%% parameters
img_fbp_sharp = reconFBP( sinoAtt, geom, 'ram-lak');
% imdisp( img_fbp_sharp, [0.18 0.22]   );

img_fbp_soft = reconFBP( sinoAtt, geom, 'hamming');
% imdisp( img_fbp_soft, [0.18 0.22]    );

%%
nitn = 50;
numos = 20;
noFrames = 40;

delta_list = [2e-4 1e-3 2e-3 3e-3 4e-3 5e-3];

for i = 1 : length( delta_list )
    
    delta = delta_list(i);
    beta1 = 5 / delta;
    beta2 = 500 / delta;
    
    [img_pwls1] = reconPwlsLALMOs14( sinoAtt, weights, geom, beta1, 'huber', nitn, delta, numos );
    %figure; imdisp( img_pwls1,  [0.19 0.23]);
    
    [img_pwls2] = reconPwlsLALMOs14( sinoAtt, weights, geom, beta2, 'huber', nitn, delta, numos );
    %figure; imdisp( img_pwls2,  [0.19 0.23]);
    
    % [img_pwls0] = reconPwlsLALMOs14( sinoAtt, weights, geom, 2 * sqrt( beta1 * beta2 ), 'huber', nitn, delta, numos );
    %figure; imdisp( img_pwls0,  [0.19 0.23]);
    
    [ img_psadmma, betas_psadmm ] = reconPwlsPathSeekingADMM( sinoAtt,weights, geom, 'huber', delta, ...
        numos / 2, img_pwls1, img_pwls2, beta1, beta2, noFrames, 4 );
    %figure; imdisp( img_psadmma,  [0.19 0.23]);
    
    save(['liver_high_' num2str(delta) '.mat'], 'img_pwls1', 'img_pwls2', 'img_psadmma', 'betas_psadmm' );
    
end




