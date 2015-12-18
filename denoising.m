
sinoAtt_CV = zeros(768, 512, 625, 'single');

for iv = 1:625
     fprintf('Filter in image slice %2.2f \n', iv);
     
     noissi2d = (sinoAtt_aver(:, :, iv) - sinoAttAir_aver(:, :, iv)); 
     de = gauss_denoise_cvt1(noissi2d(1:512, :));
     noissi2d(1:512, :) = de;
     sinoAtt_CV(:, :, iv) = noissi2d;
end
    
 