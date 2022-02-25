% Chad Cole
% Oct 21 2021
% 
% This script generates data structures and input test files for
%  the OpenCL erasure codec

% % % % Assume desired Clist/Vlist in work space before running script
% [fid_H,msg] = fopen('n2000_k1000_no6cycle_ldpc_H.h','wt');
%  
% assert(fid_H>=3,msg);
% fprintf(fid_H, 'short parity_check_mat[%d][%d] = \n { \n', length(Vlist(:, 1)), max((Vlist(:, 1))));
% 
% for ii=1:length(Vlist(:, 1))
%     fprintf(fid_H,'{');
%     for jj=2:(Vlist(ii, 1))
%          fprintf(fid_H,'%d, ', Vlist(ii, jj));
%     end
%     fprintf(fid_H,'%d ', Vlist(ii, (Vlist(ii, 1))+1));
%     
%     fprintf(fid_H,'},\n');
% end
% fprintf(fid_H,'};');
% 
% fclose(fid_H);
% 

% % % Assume desired Clist/Vlist in work space before running script
% % % Put the H Matrix in Vlist format since the number of variable node
% % % connections is not uniform, want first value to indicate how many
% % % non-zero v-node connections in that parity check
% [fid_H,msg] = fopen('n2000_k1000_no6cycle_ldpc_Vlist.h','wt');
[fid_H,msg] = fopen('n2040_k1530_no6cycle_ldpc_Vlist.h','wt');
 
assert(fid_H>=3,msg);
len_V = length(Vlist(1, :));
fprintf(fid_H, 'short parity_check_mat_Vlist[%d][%d] = \n { \n', length(Vlist(:, 1)), len_V);
for ii=1:length(Vlist(:, 1))
    row_vals = zeros(1, 20);  % We want to zero pad for up to 19 parity check degree
    row_vals(1:len_V) = Vlist(ii, :);
    fprintf(fid_H,'{');
    for jj=1:19 % hard code to length 20 size of row_vals
         fprintf(fid_H,'%d, ', row_vals(jj));
    end
    fprintf(fid_H,'%d ', row_vals(20));
    
    fprintf(fid_H,'},\n');
end
fprintf(fid_H,'};');

fclose(fid_H);


% % % This code will generate enough test input data for one block for the encoder
% k = 1000;
% max_ushort = 2^16 - 1;
% file_str = sprintf("LDPC_ErasureEncoder_IN_k%d_Shorts.txt", k);
% [fid_H,msg] = fopen(file_str,'wt');
% assert(fid_H>=3,msg);
% 
% for ii=1:k
%     rand_sym = round(rand(1)*max_ushort);
%     fprintf(fid_H,'%d \n', rand_sym);
% end
% fclose(fid_H);

% % % This code will generate corresponding test output data for the above block for the encoder
% % % (see source_encode_vec from LDPCErasureCodes_MessagePassingAlgSim) 
% n = length(source_encode_vec);
% file_str = sprintf("LDPC_ErasureEncoder_OUT_n%d_Shorts.txt", n);
% [fid_H,msg] = fopen(file_str,'wt');
% assert(fid_H>=3,msg);
% 
% for ii=1:n
%     fprintf(fid_H,'%d \n', source_encode_vec(ii));
% end
% fclose(fid_H);

% % % Create the noisy input codeword for the decoder, I verified that
% % % recv_vec contained no '0's before subbing in that value for erasures
% file_str = sprintf("LDPC_ErasureDecoder_IN_n%d_Shorts_%dPercentPER.txt", n, 100*PER);
% [fid_H,msg] = fopen(file_str,'wt');
% assert(fid_H>=3,msg);
% for ii=1:n
%     if (recv_vec(ii) == -1) %erasure, change to '0' for OpenCL since it is supposed to be an unsigned short
%         fprintf(fid_H,'%d \n', 0);
%     else
%         fprintf(fid_H,'%d \n', recv_vec(ii));
%     end
% end
% fclose(fid_H);
% 
