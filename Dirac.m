function D = Dirac(phi,sigma)
[m,n] = size(phi);
D = zeros(m,n);

for i =1:m
    for j = 1:n
        D(i,j) = 1/pi*(sigma/(sigma^2+phi(i,j)^2));
    end
end     