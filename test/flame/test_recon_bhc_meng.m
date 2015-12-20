% load reconstruction parameters
load 'temp.mat';

geom = loadProjectionGeometryCT( p );

spectrum = loadSpectraCT(p, geom, 2e6);

Dir = 'D:\MATLAB\CTData\Dec_01_2015_Study\3ppi_sic_60KV_50mA_lowflow\';

%% load air scan data

dataPathAir = [Dir 'air_04' '\'];

process_seq_file( Dir, 'air_04' );

sinoAttAir = loadTableTopData( dataPathAir, geom, 0, [701 800], [101 500] );

% load burn scan data

dataPath = [Dir 'burn_05' '\'];

process_seq_file( Dir, 'burn_05' );

sinoAtt = loadTableTopData( dataPath, geom, 0, [701 800], [101 500] );


%% first pass reconstruction

sinoAttAirPoly = beamHardeningMaterialCorrection(sinoAttAir, spectrum, 'Quartz', 3 );

imgAir = reconFBP( sinoAttAirPoly, geom, 'hamming' );

%imgAir = reconFBP( sinoAttAir, geom, 'hamming' );

%clear sinoAttAirPoly;

%return;


%% second pass beam hardening correction

mapTube = single( imgAir > 0.6 );

sinoTube = forwardProjectMex( mapTube, geom ) ;

sinoTube = imfilter3( sinoTube, fspecial('gaussian', [5 5], 1 ) );

%%

sinoDiff = sinoAtt - sinoAttAir;

sinoAttAirBHC = beamHardeningMaterialCorrectionBurner(sinoDiff, sinoAttAir, sinoTube, spectrum);

imgAir = reconFBP( sinoAttAirBHC, geom, 'hamming' );

% 
% % final reconstruction
% 
% imgAir = reconFBP( sinoAttAirBHC, geom, 'hamming' );
% 
% figure(21); 
% if geom.reconSize(3) < 40
%     imdisp( imgAir, [0 0.5]   );
% else
%     imdisp( squeeze(imgAir(end / 2, :, : ))', [0 0.8]   );
% end


%% first pass reconstruction

 %sinoAttPoly = beamHardeningMaterialCorrection(sinoAtt, spectrum, 'Quartz', 3 );

imgKr = reconFBP( sinoAtt, geom, 'hamming' );

%clear sinoAttPoly;

%% second pass beam hardening correction

% mapTube = single( imgKr > 0.6 );
% 
% sinoTube = forwardProjectMex( mapTube, geom ) ;
% 
% sinoTube = imfilter3( sinoTube, fspecial('gaussian', [5 5], 1 ) );
% 
% sinoAttBHC = beamHardeningMaterialCorrectionBurner(sinoAtt, sinoTube, spectrum);
% 
% % final reconstruction
% 
% imgKr = reconFBP( sinoAttBHC, geom, 'hamming' );
% 
% figure(22);
% if geom.reconSize(3) < 40
%     imdisp( imgKr, [0 0.5] );
% else
%      imdisp( squeeze(imgKr(end / 2, :, : ))', [0 0.8]   );
% end

%%
figure(23)
imagesc( squeeze(imgKr(end / 2, :, : ) - imgAir(end / 2, :, : ))', [0 0.1]   );


return;

%%
imgKrReg = imgKr;

[optimizer, metric] = imregconfig('monomodal');
for i = 1 : size( imgAir, 3 )
    if mod(i, 50) == 0
        fprintf('(%i/%i)... ', i, size(imgAir, 3 ) );
    end
    slice = imregister(imgKr(:,:,i), imgAir(:, :, i), 'rigid', optimizer, metric);
    imgKrReg( :,:,i) = slice;
end

figure(23); imdisp( imgKrReg(:,end/2,:) - imgAir(:,end/2,:)  , [-0.1 0.1] );