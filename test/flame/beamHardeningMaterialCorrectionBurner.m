function sinoOut = beamHardeningMaterialCorrectionBurner(sinoIn, sinoAirScan, pathLengthQuartz, spectrum)
% A special two-pass dual material beam hardening correction algorithm
% designed for the porous median bruner projection
% The materials of the burner are Quartz and SiC, the correction is mean
% get the correct scaling for the Kr density after the subtraction of the
% reconstruction. 
% Inputs:
%   SinoIn - original uncorrected sinogram
%   pathLengthQuartz - extracted path length of the quartz tube from the
%       prior reconstruction
%   spectrum - spectrum information
% Output:
%   sinoOut - beam hardening corrected sinogram
%
%
% Copyright (c) 2015.12 by Meng Wu, Stanford University.


% default material assumed for beam hardening correction is water

if spectrum.energyIntegrating
    photonsPerEnergyBin = spectrum.DQE * spectrum.photonsPerEnergyBin .* spectrum.energyBinLabels * spectrum.detectorGain ;
else
    photonsPerEnergyBin = spectrum.DQE * spectrum.photonsPerEnergyBin * spectrum.detectorGain  ;
end

if spectrum.useBowtie
    photonsPerEnergyBin = computeResultingPhotons( photonsPerEnergyBin, ...
        spectrum.energyBinLabels, spectrum.bowtieMaterial, 10 );
end


%% get material energy-dependent attenuation coefficients
muQuartz = materialAttenuation( spectrum.energyBinLabels, 'Quartz' ) ;

muSiC = materialAttenuation( spectrum.energyBinLabels, 'SiC' );

muKr = XrayMu( 'Kr', spectrum.energyBinLabels );

muKrRef = XrayMu( 'Kr', spectrum.energyAverage + 5 );

%% estimate gain correction with spectrum info

% create a meshgrid of path lengths
quartzLengths = linspace( 0.6, 1.0, 64);
sicLengths = linspace( 0, 0.7, 64);

[quartzLengthTable,  sicLengthTable] = meshgrid(quartzLengths, sicLengths );

quartzLengthTable = quartzLengthTable(:);
sicLengthTable = sicLengthTable(:);

n = length( quartzLengthTable );

% estimated beam harden values
lineIntegralTable = zeros( n, 1);
gainCorrection = zeros( n, 1);

for i = 1 : n
    
    transmittedSpectrum = photonsPerEnergyBin .* exp( - ( muQuartz * quartzLengthTable(i) + muSiC * sicLengthTable(i) ) );
    
    lineIntegralTable(i) = log( sum(photonsPerEnergyBin) / sum( transmittedSpectrum ) );
    
    muKrEffective = log( sum( transmittedSpectrum ) / sum( transmittedSpectrum .* exp( - 0.001 * muKr ) ) ) / 0.001;
    
    gainCorrection(i) = muKrRef / muKrEffective;
    
end


mean( lineIntegralTable(:))
%% calculate beam hardening gain correction polynomial

A = [ones(n, 1) lineIntegralTable lineIntegralTable.^2 quartzLengthTable quartzLengthTable.*lineIntegralTable ];

polyCoeffs = A \ gainCorrection;

% check MSE 
sqrt( mean( (gainCorrection -  A * polyCoeffs ).^2 ) );

%  figure(11)
%  imagesc( reshape( gainCorrection -  A * polyCoeffs , [64 64] ) ); colorbar;


%% Final gain correction

sinoGainBHC = ones( size( sinoAirScan ), 'single' );
sinoGainBHC = sinoGainBHC * polyCoeffs(1) + sinoAirScan * polyCoeffs(2) + sinoAirScan.^2 * polyCoeffs(2) + pathLengthQuartz * polyCoeffs(3) + pathLengthQuartz .*sinoAirScan * polyCoeffs(4) ;

sinoGainBHC( sinoGainBHC > 3 ) = 3; 

sinoOut = sinoGainBHC .* sinoIn;

%imdisp( sinoGainBHC(:,:,end/2), [0 5])

end