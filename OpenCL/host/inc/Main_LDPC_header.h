//This ROM stores LDPC parameters to be accessed by the constel_index value
//Column 0  NLDPC         - NLDPC
//Column 1: KLDPC         - KLDPC
//Column 2: first position in parity_check_mat_Vlist
//Column 3: last position in parity_check_mat_Vlist
//Column 4: RS_n equivalent
//Column 5: RS_k equivalent


int ldpc_params[2][6] =
{
	{2000, 1000, 0, 999, 250, 125},  // 0.  (2000, 1000)
	{2040, 1530, 1000, 1509, 255, 192},  // 1.  (2040, 1530)
};

