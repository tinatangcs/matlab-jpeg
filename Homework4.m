clear;
% [I,map] = imread('River.gif');
[I,map] = imread('Lena.gif');

G = ind2gray(I,map);
G = double(G);
scalar = 10;
[header,data] = MyJPEG(G,scalar);
rebuildImage = MyJPEGDecoder(header,data,size(G),scalar);
SNR = snr(G,rebuildImage)
MSE = mse(G,rebuildImage)
compression_len = length(data.DCbitstream) + length(data.ACbitstream);
original_len = size(G,1)*size(G,2)*8;
compression_ratio = original_len/compression_len
imagesc(rebuildImage);
colormap(gray);