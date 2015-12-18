% load reconstruction parameters
load 'temp.mat';

geom = loadProjectionGeometryCT( p );

% geom.detSize(2) = 256;
% geom.reconSize(3) = 256;
% geom.reconOffset(3) = 0;

spectrum = loadSpectraCT(p, geom, 2e6);

dir = 'D:\MATLAB\CTData\Dec_01_2015_Study\3ppi_sic_60KV_50mA_lowflow\';

%% load air scan data

dataPathAir = [dir 'air_05' '\'];
process_seq_file( dir, 'air_05' );

sinoAttAir = loadTableTopData( dataPathAir, geom );

%% load burn scan data

dataPath = [dir 'burn_04' '\'];
process_seq_file( dir, 'burn_04' );

sinoAtt = loadTableTopData( dataPath, geom );

%% first pass reconstruction

sinoAttAirPoly = beamHardeningMaterialCorrection(sinoAttAir, spectrum, 'Quartz', 3 );

imgAir = reconFBP( sinoAttAirPoly, geom, 'hamming' );

clear sinoAttAirPoly;

%% second pass beam hardening correction

mapTube = single( imgAir > 0.6 );

sinoTube = forwardProjectMex( mapTube, geom ) ;

sinoTube = imfilter3( sinoTube, fspecial('gaussian', [5 5], 1 ) );

sinoAttAirBHC = beamHardeningMaterialCorrectionBurner(sinoAttAir, sinoTube, spectrum);

% final reconstruction

imgAir = reconFBP( sinoAttAirBHC, geom, 'hamming' );

figure(21); 
if geom.reconSize(3) < 40
    imdisp( imgAir, [0 0.5]   );
else
    imdisp( squeeze(imgAir(end / 2, :, : ))', [0 0.8]   );
end


%% first pass reconstruction

sinoAttPoly = beamHardeningMaterialCorrection(sinoAtt, spectrum, 'Quartz', 3 );

imgKr = reconFBP( sinoAttPoly, geom, 'hamming' );

clear sinoAttPoly;

%% second pass beam hardening correction

mapTube = single( imgKr > 0.6 );

sinoTube = forwardProjectMex( mapTube, geom ) ;

sinoTube = imfilter3( sinoTube, fspecial('gaussian', [5 5], 1 ) );

sinoAttBHC = beamHardeningMaterialCorrectionBurner(sinoAtt, sinoTube, spectrum);

% final reconstruction

imgKr = reconFBP( sinoAttBHC, geom, 'hamming' );

figure(22);
if geom.reconSize(3) < 40
    imdisp( imgKr, [0 0.5] );
else
     imdisp( squeeze(imgKr(end / 2, :, : ))', [0 0.8]   );
end

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