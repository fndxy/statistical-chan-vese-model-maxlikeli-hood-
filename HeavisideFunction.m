function H = HeavisideFunction(phi,sigma)
[m,n] = size(phi);
H = zeros(m,n);

for i =1:m
    for j = 1:n
        H(i,j) = 0.5*(1+2/pi*atan(phi(i,j)/sigma));
    end
end
 