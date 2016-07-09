function phi = pixel_wise(Img,phi0,epsilon,mu,nu,timestep)
phi=NeumannBound(phi0);
K=curvature(phi); 
H=HeavisideFunction(phi,epsilon);
Delta = Dirac(phi,epsilon);


I = Img;
u1 = sum(sum(H.*I))/sum(sum(H));
u2 = sum(sum((1-H).*I))/sum(sum(1-H));

sigma1 = sum(sum((I - u1).^2.*H))/sum(sum(H));
sigma2 = sum(sum((I - u1).^2.*(1-H)))/sum(sum(1-H));
sigma1 = sigma1 + eps;
sigma2 = sigma2 + eps;


sub1 = (u1 - Img).^2;
sub2 = (u2 - Img).^2;
e1 = log(sqrt(2*pi*sigma1))+(sub1/sigma1);
e2 = log(sqrt(2*pi*sigma2))+(sub2/sigma2);

laplacian_term = del2(phi);

[hist_f,hist_b] = histForeBack(Img,phi);

[m,n] = size(phi);
P_Mf = zeros(m,n);
P_Mb = zeros(m,n);
for i = 1:m
    for j = 1:n
        hist_y_f = hist_f(Img(i,j)+1);
        hist_y_b = hist_b(Img(i,j)+1);
        P_Mf(i,j) = hist_y_f/(hist_y_f+hist_y_b);
        P_Mb(i,j) = hist_y_b/(hist_y_f+hist_y_b);
    end
end
P_Mf = P_Mf + 0.01;
P_Mb = P_Mb + 0.01;
A = Delta.*(e2-e1+P_Mf-P_Mb);


% A = Delta.*(localForce);
P = 1/(50)*(laplacian_term-K);
phi = phi + timestep*25*(A+P);

end
function g = NeumannBound(f)
[nrow,ncol] = size(f);
g = f;
g([1 nrow],[1 ncol]) = g([3 nrow-2],[3 ncol-2]);
g([1 nrow],2:end-1) = g([3 nrow-2],2:end-1);
g(2:end-1,[1 ncol]) = g(2:end-1,[3 ncol-2]);
end

function [hist_f, hist_b] = histForeBack(Img,phi)  %求 前景、背景 直方图
[m,n] = size(phi);
hist_f = zeros(1,256);
hist_b = zeros(1,256);
for i = 1:m
    for j = 1:n
        temp = Img(i,j);
        if phi(i,j)>0 %前景
            hist_f(temp+1) = hist_f(temp+1) + 1;
        elseif(phi(i,j)<0) %背景
            hist_b(temp+1) = hist_b(temp+1) + 1;
        end
    end
end
end

