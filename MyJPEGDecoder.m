function [output] = MyJPEGDecoder(header,data,imgsize,scalar)
DCquantizer = [0,1,2^1,2^2,2^3,2^4,2^5,2^6,2^7,2^8,2^9,2^10,2^11];
ACquantizer = [1,2^1,2^2,2^3,2^4,2^5,2^6,2^7,2^8,2^9,2^10];
dc_huff = header.DCHuffuman;
dc_h = data.DCStream.h;
dc_s = data.DCStream.s;
dc_m = data.DCStream.m;
ac_huff = header.ACHuffuman;
ac_h = data.ACStream.h;
ac_s = data.ACStream.s;
ac_m = data.ACStream.m;
quanMatrix = [
    16	11	10	16	24	40	51	61
    12	12	14	19	26	58	60	55
    14	13	16	24	40	57	69	56
    14	17	22	29	51	87	80	62
    18	22	37	56	68	109	103	77
    24	35	55	64	81	104	113	92
    49	64	78	87	103	121	120	101
    72	92	95	98	112	100	103	99];
% first rebuild DC terms
for i = 1:length(dc_h)
    temp = dc_h{i};
    dc_s{i} = bin2dec(dc_s{i});
    for j = 1:size(dc_huff,1)
        if  strcmpi(temp,dc_huff{j,2})
            dc_k{i} = dc_huff{j,1};
        end
    end
    if strcmpi(dc_m{i},'-')
        dc_t{i} = 0;
    else
        dc_t{i} = bin2dec(dc_m{i});
    end
    if dc_s{i} == 0
        dcresidual(i) = -1 * (DCquantizer(dc_k{i}+1) + dc_t{i});
    else
        dcresidual(i) = 1 * DCquantizer(dc_k{i}+1) + dc_t{i};
    end
    
    
end

dcterm(1) = dcresidual(1);
for i = 2:length(dcresidual)
    dcterm(i) = dcterm(i - 1) + dcresidual(i);
end
cell_row = imgsize(1)/8;
cell_col = imgsize(2)/8;
dc_matrix = reshape(dcterm,[cell_col,cell_row])';
% then rebuild ac terms
for i = 1:length(ac_h)
    for j = 1:length(ac_h{i})
        temp = ac_h{i}{j};
        if strcmpi(ac_s{i}{j},'1')
            ac_s{i}{j} = 1;
        elseif strcmpi(ac_s{i}{j},'0')
            ac_s{i}{j} = 0;
        end
        if strcmpi(ac_m{i}{j},'-')
            ac_t{i}{j} = 0;
        elseif ac_m{i}{j} == 9999
            ac_t{i}{j} = 9999;
        else
            ac_t{i}{j} = bin2dec(ac_m{i}{j});
        end
        for l = 1:size(ac_huff,1)
            if strcmpi(ac_h{i}{j},ac_huff(l,2))
                ac_k{i}{j} = ac_huff{l,1}(2);
                ac_d{i}{j} = ac_huff{l,1}(1);
            end
        end
        if ac_s{i}{j} == 0
            pair{i}{j} = [ac_d{i}{j},-1*(ACquantizer(ac_k{i}{j})+ac_t{i}{j})];
        elseif ac_s{i}{j} == 1
            pair{i}{j} = [ac_d{i}{j},1*(ACquantizer(ac_k{i}{j})+ac_t{i}{j})];
        else
            pair{i}{j} = 9999;
        end
    end
    
end
mid_pair = {};
for i = 1:length(pair)
    mid_pair{i} = [];
    for j = 1:length(pair{i})
        if pair{i}{j} == 9999
            mid_pair{i} = [mid_pair{i},9999];
        else
            for l = 1:(pair{i}{j}(1))
                mid_pair{i} = [mid_pair{i},0];
            end
            mid_pair{i} = [mid_pair{i},pair{i}{j}(2)];
        end
        
    end
end
for i = 1:length(mid_pair)
    templen = length(mid_pair{i});
    for j = templen:63
        mid_pair{i}(j) = 0;
    end
end
ac_matrix = reshape(mid_pair,[imgsize(2)/8,imgsize(1)/8])';
good_matrix_cell = {};
for i = 1:cell_row
    for j = 1:cell_col
        good_matrix_cell{i,j} = zeros(8);
        good_matrix_cell{i,j}(1,1) = dc_matrix(i,j);
        good_matrix_cell{i,j}(1,2) = ac_matrix{i,j}(1);
        good_matrix_cell{i,j}(2,1) = ac_matrix{i,j}(2);
        good_matrix_cell{i,j}(3,1) = ac_matrix{i,j}(3);
        good_matrix_cell{i,j}(2,2) = ac_matrix{i,j}(4);
        good_matrix_cell{i,j}(1,3) = ac_matrix{i,j}(5);
        good_matrix_cell{i,j}(1,4) = ac_matrix{i,j}(6);
        good_matrix_cell{i,j}(2,3) = ac_matrix{i,j}(7);
        good_matrix_cell{i,j}(3,2) = ac_matrix{i,j}(8);
        good_matrix_cell{i,j}(4,1) = ac_matrix{i,j}(9);
        good_matrix_cell{i,j}(5,1) = ac_matrix{i,j}(10);
        good_matrix_cell{i,j}(4,2) = ac_matrix{i,j}(11);
        good_matrix_cell{i,j}(3,3) = ac_matrix{i,j}(12);
        good_matrix_cell{i,j}(2,4) = ac_matrix{i,j}(13);
        good_matrix_cell{i,j}(1,5) = ac_matrix{i,j}(14);
        good_matrix_cell{i,j}(1,6) = ac_matrix{i,j}(15);
        good_matrix_cell{i,j}(2,5) = ac_matrix{i,j}(16);
        good_matrix_cell{i,j}(3,4) = ac_matrix{i,j}(17);
        good_matrix_cell{i,j}(4,3) = ac_matrix{i,j}(18);
        good_matrix_cell{i,j}(5,2) = ac_matrix{i,j}(19);
        good_matrix_cell{i,j}(6,1) = ac_matrix{i,j}(20);
        good_matrix_cell{i,j}(7,1) = ac_matrix{i,j}(21);
        good_matrix_cell{i,j}(6,2) = ac_matrix{i,j}(22);
        good_matrix_cell{i,j}(5,3) = ac_matrix{i,j}(23);
        good_matrix_cell{i,j}(4,4) = ac_matrix{i,j}(24);
        good_matrix_cell{i,j}(3,5) = ac_matrix{i,j}(25);
        good_matrix_cell{i,j}(2,6) = ac_matrix{i,j}(26);
        good_matrix_cell{i,j}(1,7) = ac_matrix{i,j}(27);
        good_matrix_cell{i,j}(1,8) = ac_matrix{i,j}(28);
        good_matrix_cell{i,j}(2,7) = ac_matrix{i,j}(29);
        good_matrix_cell{i,j}(3,6) = ac_matrix{i,j}(30);
        good_matrix_cell{i,j}(4,5) = ac_matrix{i,j}(31);
        good_matrix_cell{i,j}(5,4) = ac_matrix{i,j}(32);
        good_matrix_cell{i,j}(6,3) = ac_matrix{i,j}(33);
        good_matrix_cell{i,j}(7,2) = ac_matrix{i,j}(34);
        good_matrix_cell{i,j}(8,1) = ac_matrix{i,j}(35);
        good_matrix_cell{i,j}(8,2) = ac_matrix{i,j}(36);
        good_matrix_cell{i,j}(7,3) = ac_matrix{i,j}(37);
        good_matrix_cell{i,j}(6,4) = ac_matrix{i,j}(38);
        good_matrix_cell{i,j}(5,5) = ac_matrix{i,j}(39);
        good_matrix_cell{i,j}(4,6) = ac_matrix{i,j}(40);
        good_matrix_cell{i,j}(3,7) = ac_matrix{i,j}(41);
        good_matrix_cell{i,j}(2,8) = ac_matrix{i,j}(42);
        good_matrix_cell{i,j}(3,8) = ac_matrix{i,j}(43);
        good_matrix_cell{i,j}(4,7) = ac_matrix{i,j}(44);
        good_matrix_cell{i,j}(5,6) = ac_matrix{i,j}(45);
        good_matrix_cell{i,j}(6,5) = ac_matrix{i,j}(46);
        good_matrix_cell{i,j}(7,4) = ac_matrix{i,j}(47);
        good_matrix_cell{i,j}(8,3) = ac_matrix{i,j}(48);
        good_matrix_cell{i,j}(8,4) = ac_matrix{i,j}(49);
        good_matrix_cell{i,j}(7,5) = ac_matrix{i,j}(50);
        good_matrix_cell{i,j}(6,6) = ac_matrix{i,j}(51);
        good_matrix_cell{i,j}(5,7) = ac_matrix{i,j}(52);
        good_matrix_cell{i,j}(4,8) = ac_matrix{i,j}(53);
        good_matrix_cell{i,j}(5,8) = ac_matrix{i,j}(54);
        good_matrix_cell{i,j}(6,7) = ac_matrix{i,j}(55);
        good_matrix_cell{i,j}(7,6) = ac_matrix{i,j}(56);
        good_matrix_cell{i,j}(8,5) = ac_matrix{i,j}(57);
        good_matrix_cell{i,j}(8,6) = ac_matrix{i,j}(58);
        good_matrix_cell{i,j}(7,7) = ac_matrix{i,j}(59);
        good_matrix_cell{i,j}(6,8) = ac_matrix{i,j}(60);
        good_matrix_cell{i,j}(7,8) = ac_matrix{i,j}(61);
        good_matrix_cell{i,j}(8,7) = ac_matrix{i,j}(62);
        good_matrix_cell{i,j}(8,8) = ac_matrix{i,j}(63);
    end
end
for i = 1:size(good_matrix_cell,1)
    for j = 1:size(good_matrix_cell,2)
        good_matrix_cell{i,j} = good_matrix_cell{i,j}.*(scalar*quanMatrix);
        good_matrix_cell{i,j} = (idct2(good_matrix_cell{i,j}));
        good_matrix_cell{i,j} = good_matrix_cell{i,j} + 128;
    end
end
rebuild = cell2mat(good_matrix_cell);
% imagesc(finally);
% colormap(gray);
output = rebuild;
end

