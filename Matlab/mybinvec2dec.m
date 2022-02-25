function [out] = mybinvec2dec(dec, n)

temp = 0;
  for ii = 1:n
      temp = temp + dec(ii)*2^(n-ii);
  end
out = temp;