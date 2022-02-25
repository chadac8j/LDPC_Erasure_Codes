% Reed Solomon Erasure Codes Decoder
% Chad Cole
% Aug 10th, 2021
%
% This helper function builds the Multiply, Add and Inverse lookup tables for GF(2^8)

function [GF_add_lookup, GF_mult_lookup, gf_inv, G] = Build_GF256_Lookup_Tables(m, k, n, prim_poly_m8)

GF_SIZE = 2^m;

prim_poly_m8_number = 0;
for ii = 1:length(prim_poly_m8)
    prim_poly_m8_number = prim_poly_m8_number + prim_poly_m8(length(prim_poly_m8) - ii + 1)*2^(ii-1);
end
% Form log table
gftable(m, prim_poly_m8_number);
x_gf = gf(0:2^m-1, m, prim_poly_m8_number);


% Form log table
gf_log_inv = zeros(1, GF_SIZE-1);
gf_log_inv(1) = 0;
gf_log_inv(2) = 1;
alpha = x_gf(3);
gf_temp = alpha;
for ii = 3:GF_SIZE
    gf_log_inv(ii) = double(gf_temp.x);
    gf_temp = alpha*gf_temp;
end
%Appears to be a bug in the log() func for GF
Log_array = [-inf 0:GF_SIZE-2];
log_lookup(gf_log_inv+1) = Log_array;
% gf_log_inv = gf(gf_log_inv, m, prim_poly_m8_number);

gf_inv = zeros(1, GF_SIZE);
gf_inv(1) = 0;
gf_inv(2) = 1;
for ii = 3:GF_SIZE
	gf_inv(ii) = gf_log_inv(GF_SIZE-log_lookup(ii)+1);
end
gf_inv = gf_inv(2:end);

% Unit test for my GF(256) lookup tables vs Matlab built in implementation
GF_mult_lookup = zeros(GF_SIZE, GF_SIZE);
for row = 2:length(x_gf)
    for col = 2:length(x_gf)
        matlab_imp = x_gf(row)*x_gf(col);
        my_imp = gf_log_inv(mod(log_lookup(double(x_gf.x(row)) + 1) + log_lookup(double(x_gf.x(col)) + 1), 255) + 2);
        GF_mult_lookup(row, col) = my_imp;
        if (matlab_imp.x ~= my_imp)
            my_imp
        end
    end
end

% Unit test for bitxor addition/subtraction logic vs Matlab built in GF implementation
GF_add_lookup = zeros(GF_SIZE, GF_SIZE);
for row = 1:length(x_gf)
    for col = 1:length(x_gf)
        matlab_imp = x_gf(row) + x_gf(col);
        my_imp = bitxor( double(x_gf.x(row)), double(x_gf.x(col)) );
        GF_add_lookup(row, col) = my_imp;
        if (matlab_imp.x ~= my_imp)
            my_imp
        end
    end
end

G = gf(zeros(k, n), m, prim_poly_m8_number);
%% This code takes a long time, so if the G is not used, don't worry about it
% for row=1:k
%     for col=1:n
%         G(row, col) = alpha^(col*row);
% %         G(row, col) = gf_log(mod(col*row, 256)+1);
%     end
% end
