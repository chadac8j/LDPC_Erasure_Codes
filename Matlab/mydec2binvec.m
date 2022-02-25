%#eml
function [out] = mydec2binvec(jj, n)
out = zeros(1,n);
temp = jj;
  for ii = n:-1:1
      if (temp - 2^(ii-1) >= 0)
          out(n-ii+1) = 1;
          temp = temp - 2^(ii-1);
      end
  end
  