function [header,data] = MyJPEG(image,scalar)
quanMatrix = [
    16	11	10	16	24	40	51	61
    12	12	14	19	26	58	60	55
    14	13	16	24	40	57	69	56
    14	17	22	29	51	87	80	62
    18	22	37	56	68	109	103	77
    24	35	55	64	81	104	113	92
    49	64	78	87	103	121	120	101
    72	92	95	98	112	100	103	99];
blocksize = [8,8];
image = double(image);
sub128 = image - 128;
applydct = blkproc(sub128, blocksize, @dct2);
applyquan = blkproc(applydct,blocksize,@quantization,scalar*quanMatrix);

%     blkproc
%     ggg = blkproc(applydct,blocksize,@bbb);
DCtermMat = blkproc(applyquan,blocksize,@getDCterm);
DCterm = reshape(DCtermMat',[1],[]);
ACtermMat = blkproc(applyquan,blocksize,@getACterm);
for i = 1:size(ACtermMat,1)
    for j = 1:size(ACtermMat,2)/63
        ACterm{i,j} = ACtermMat(i,63*(j-1)+1:j*63);
    end
end

% First Code DC Residuals
DCresiduals(1) = DCterm(1,1);
for i = 2:size(DCterm,2)
    DCresiduals(i) = DCterm(1,i) - DCterm(1,i-1);
end
dc_symbols = [0:11];
DCquantizer = [0,1,2^1,2^2,2^3,2^4,2^5,2^6,2^7,2^8,2^9,2^10,2^11];
[dc_k,dc_t,dc_s,dc_times] = quantizeDC(DCresiduals,DCquantizer);
dc_prob = dc_times/size(DCresiduals,2);
[DCdict,avglen] = huffmandict(dc_symbols,dc_prob);
for i = 1:size(DCdict,1)
    DCdict{i,2} = num2str(DCdict{i,2});
    DCdict{i,2} = DCdict{i,2}(~isspace(DCdict{i,2}));
end
for i = 1:size(dc_k,2)
    dc_h{1,i} = DCdict{dc_k(1,i)+1,2};
    
end
for i = 1:size(dc_t,2)
    temp = dec2bin(dc_t(1,i));
    l = length(temp);
    k = dc_k(1,i);
    if k == 0 || k == 1
        dc_m{1,i} = '-';
    else
        morebit = k - 1 - l;
        for j = 1:morebit
            temp = ['0',temp];
        end
        dc_m{1,i} = temp;
    end
end

%Then Code AC terms
ACquantizer = [1,2^1,2^2,2^3,2^4,2^5,2^6,2^7,2^8,2^9,2^10];

for i = 1:size(ACterm,1)
    for j = 1:size(ACterm,2)
        for l = length(ACterm{i,j}):-1:0
            if l == 0
                ACterm{i,j}(1) = 9999;
                break;
            elseif ACterm{i,j}(l) == 0
                ACterm{i,j}(l) = [];
                continue;
            else
                ACterm{i,j}(l+1) = 9999;
                break;
            end
        end
        zero_runs = 0;
        index = 1;
        for l = 1:length(ACterm{i,j})%calculate zero runs tuple
            if ACterm{i,j}(l) == 9999
                AC_Pairs{i,j}{index} = 9999;
            elseif ACterm{i,j}(l) == 0
                zero_runs = zero_runs + 1;
            else
                AC_Pairs{i,j}{index} = [zero_runs,ACterm{i,j}(l)];
                zero_runs = 0;
                index = index + 1;
            end
        end
    end
end

AC_Good_Pair = reshape(AC_Pairs',1,[]);
AC_NUM = getACNum(AC_Good_Pair);
[ac_d,ac_k,ac_t,ac_s] = quantizeAC(AC_Good_Pair,ACquantizer);
% build huffman tree for (d,k) tuple
ac_symbol = {[9999,9999]};
ac_times = [0];
d_k_tuple = [];
for i = 1:length(ac_d)
    for j = 1:length(ac_d{i})
        flag = 0;
        d_k_tuple = [ac_d{i}(j),ac_k{i}(j)];
        for l = 1:length(ac_symbol)
            if ac_symbol{l} == d_k_tuple
                ac_times(l) = ac_times(l) + 1;
                flag = 1;
                break;
            end
        end
        if flag == 0
            ac_symbol = [ac_symbol,d_k_tuple];
            ac_times = [ac_times,1];
        end
    end
end
ac_prob = ac_times/AC_NUM;
[ACdict,ac_avglen] = huffmandict(ac_symbol,ac_prob);
for i = 1:size(ACdict,1)
    ACdict{i,2} = num2str(ACdict{i,2});
    ACdict{i,2} = ACdict{i,2}(~isspace(ACdict{i,2}));
end
ac_h = {};
for i = 1:length(ac_d)
    for j = 1:length(ac_d{i})
        temp = [ac_d{i}(j),ac_k{i}(j)];
        for l = 1:size(ACdict,1)
            if ACdict{l,1} == temp;
                ac_h{i}{j} = ACdict{l,2};
            end
        end
    end
end


for i = 1:length(ac_t)
    for j = 1:length(ac_t{i})
        if j ~=length(ac_t{i})
            temp = dec2bin(ac_t{i}(j));
            l = length(temp);
            k  = ac_k{i}(j);
            if k == 1
                ac_m{i}{j} = '-';
            else
                morebit = k - 1 - l;
                for q = 1:morebit
                    temp = ['0',temp];
                end
                ac_m{i}{j} = temp;
            end
        else
            ac_m{i}{j} = ac_t{i}(j);
        end
        
    end
end

% get DC bitstream
dc_bitstream = [];
ac_bitstream = [];
for i = 1:length(dc_h)
    if dc_m{i} == '-'
        dc_bitstream = [dc_bitstream,dc_h{i},dc_s{i}];
    else
        dc_bitstream = [dc_bitstream,dc_h{i},dc_s{i},dc_m{i}];
    end
    
end


for i = 1:length(ac_h)
    for j = 1:length(ac_h{i});
        if ac_m{i}{j} == 9999
            ac_bitstream = [ac_bitstream,ac_h{i}{j}];
        elseif ac_m{i}{j} == '-'
            ac_bitstream = [ac_bitstream,ac_h{i}{j},ac_s{i}{j}];
        else
            ac_bitstream = [ac_bitstream,ac_h{i}{j},ac_s{i}{j},ac_m{i}{j}];
        end
    end
end
fid1 = fopen('dc_bitstream.txt','w+');
fprintf(fid1,dc_bitstream);
fid2 = fopen('ac_bitstream.txt','w+');
fprintf(fid2,ac_bitstream);


% --------------------
% test image effect
% back1 = blkproc(applyquan,blocksize,@dequantization,scalar*quanMatrix);
% back2 = blkproc(back1,blocksize,@idct2);
% back3 = round(back2 + 128);
% imagesc(back3);
% colormap(gray);
% SNR = snr(image,back3);
% MSE = mse(image,back3);
% --------------------
header.DCHuffuman = DCdict;
header.ACHuffuman = ACdict;
data.DCStream.h = dc_h;
data.DCStream.s = dc_s;
data.DCStream.m = dc_m;
data.ACStream.h = ac_h;
data.ACStream.s = ac_s;
data.ACStream.m = ac_m;
data.DCbitstream = dc_bitstream;
data.ACbitstream = ac_bitstream;



end

function output = quantization(mat,qmat)
output = round(mat./qmat);
end

function output = dequantization(mat,qmat)
output = round(mat.*qmat);
end

function output = getDCterm(mat)
output = mat(1,1);
end

function output = getACterm(mat)
temp = zigzag(mat);
output = temp;
end

function [k,t,s,times] = quantizeDC(dc,dcquantizer)
k = [];
t = [];
s = [];
times = zeros(1,12);
for i = 1:size(dc,2);
    %     if dc(1,i) < 0
    %         s = [s,0];
    %     else
    %         s = [s,1];
    %     end
    if dc(1,i) < 0
        s{1,i} = '0';
    else
        s{1,i} = '1';
    end
    temp = abs(dc(1,i));
    for j = 1:size(dcquantizer,2) - 1
        if (temp >= dcquantizer(1,j) && temp < dcquantizer(1,j+1))
            k = [k,j-1];
            tempt = temp - dcquantizer(1,j);
            t = [t,tempt];
            times(j) = times(j) + 1;
        end
    end
end
end
function num = getACNum(ac_good_pair)
num = 0;
for i = 1:length(ac_good_pair)
    num = num + length(ac_good_pair{i});
end
end
function [d,k,t,s] = quantizeAC(ac_good_pair,acquantizer)
d = {};
k = {};
s = {};
t = {};
for i = 1:length(ac_good_pair)
    for j = 1:length(ac_good_pair{i})
        if j == length(ac_good_pair{i})
            s{i}{j} = 9999;
        elseif ac_good_pair{i}{j}(2) < 0
            s{i}{j} = '0';
        else
            s{i}{j} = '1';
        end
        if ac_good_pair{i}{j}(1) ~= 9999
            temp = abs(ac_good_pair{i}{j}(2));
            d{i}(j) = ac_good_pair{i}{j}(1);
            for m = 1:length(acquantizer)-1
                if temp >= acquantizer(m) && temp < acquantizer(m + 1)
                    k{i}(j) = m;
                    tempt = temp - acquantizer(m);
                    t{i}(j) =  tempt;
                end
                
            end
        else
            k{i}(j) = 9999;
            t{i}(j) = 9999;
            d{i}(j) = 9999;
        end
    end
    
end
end
