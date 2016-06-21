clear;
% scalar = [1.66,2.21,3.08];
% compression_ratio = [12.0150,15.0304,20.0290];
% snr = [0.0104,0.0113,0.0132];


scalar = [0.89,1.23,1.82];
compression_ratio = [11.9862,15.0934,19.9950];
snr = [0.0029,0.0036,0.0035];
plot(scalar,compression_ratio);
xlabel('Scalar');
ylabel('Compression_Ratio');
title('Compression_Ratio and Scalar Graph');

plot(scalar,snr);
xlabel('Scalar');
ylabel('SNR');
title('SNR and Scalar Graph');
