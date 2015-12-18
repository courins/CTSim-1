% load reconstruction parameters
load 'temp.mat';

geom = loadProjectionGeometryCT( p );

% geom.detSize(2) = 256;
% geom.reconSize(3) = 256;
% geom.reconOffset(3) = 0;

spectrum = loadSpectraCT(p, geom, 2e6);

my_dir = 'D:\TabletopScannerData\Dec_14_2015_Study\Air\';

%% load air scan data

% dataPathAir = [my_dir 'air_16' '\'];
% 
% process_seq_file( my_dir, 'air_16' );
% 
% sinoAttAir = loadTableTopData( dataPathAir, geom );


validPixelsX = geom.detSize(1);
validPixelsY = geom.detSize(2);
noViews = geom.noViews;
sinoAttSum = zeros(validPixelsY, validPixelsX, noViews, 'single');

numOfData = 15;
counter = 1;
files = dir('D:\TabletopScannerData\Dec_14_2015_Study\Air\*.seq');
for file = files'
    seq = file.name;
    fileLabel = seq(1:end-4)


    dataPath = [my_dir fileLabel '\'];
    process_seq_file( my_dir, fileLabel );
    
    sinoAttAir = loadTableTopData( dataPath, geom );

    sinoAttSum =  sinoAttSum + sinoAttAir; 
    
    if (counter >= numOfData)
        break;
    end
    
    counter = counter + 1;
    
end

sinoAttAir = sinoAttSum./counter;

%% load burn scan data
% valid detector size defined by geometry paramters
validPixelsX = geom.detSize(1);
validPixelsY = geom.detSize(2);
noViews = geom.noViews;
sinoAttSum = zeros(validPixelsY, validPixelsX, noViews, 'single');

% files = dir('D:\TabletopScannerData\Dec_14_2015_Study\*.seq');
% for file = files'
%     seq = file.name;
%     fileLabel = seq(1:end-4)
% 
% 
%     dataPath = [my_dir fileLabel '\'];
%     process_seq_file( my_dir, fileLabel );
% 
% end

% dataPath = [my_dir 'burn_16' '\'];
% 
% process_seq_file( my_dir, 'burn_16' );
% 
% sinoAtt = loadTableTopData( dataPath, geom );

% We can load multiple datasets and average them

my_dir = 'D:\TabletopScannerData\Dec_14_2015_Study\Burn\';
counter = 1;
files = dir('D:\TabletopScannerData\Dec_14_2015_Study\Burn\*.seq');
for file = files'
    seq = file.name;
    fileLabel = seq(1:end-4)


    dataPath = [my_dir fileLabel '\'];
    process_seq_file( my_dir, fileLabel );
    
    sinoAtt = loadTableTopData( dataPath, geom );

    sinoAttSum =  sinoAttSum + sinoAtt; 
    
    if (counter >= numOfData)
        break;
    end
    
    counter = counter + 1;
    
end

sinoAtt = sinoAttSum./counter;

% Test: Denoising using wavelets

% for iv = 1:noViews
%     noissi2d = sinoAtt(:, :, iv); 
% 
%     [XDEN,~,~] = func_denoise_dw2d_v2(noissi2d);
% 
%     if iv == 1
%         figure;
%         imagesc(XDEN-noissi2d);
%         colormap jet;
%         
%         fid = fopen('Before.bin','w');
%         fwrite(fid, noissi2d, 'single');
%         
%         fid = fopen('After.bin','w');
%         fwrite(fid, XDEN, 'single');
%         
%     end
%     sinoAtt(:, :, iv) = XDEN;
% end


%% first pass reconstruction

%sinoAttAirPoly = beamHardeningMaterialCorrection(sinoAttAir, spectrum, 'Quartz', 3 );

%imgAir = reconFBP( sinoAttAirPoly, geom, 'hamming' );

imgAir = reconFBP( sinoAttAir, geom, 'hamming' );

%clear sinoAttAirPoly;


%return;


%% second pass beam hardening correction

% mapTube = single( imgAir > 0.6 );
% 
% sinoTube = forwardProjectMex( mapTube, geom ) ;
% 
% sinoTube = imfilter3( sinoTube, fspecial('gaussian', [5 5], 1 ) );
% 
% sinoAttAirBHC = beamHardeningMaterialCorrectionBurner(sinoAttAir, sinoTube, spectrum);
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