function f = caculus(I,phi,sigma,Delta)
phi=NeumannBoundCond(phi);
[phi_x , phi_y] = gradient(phi);
[m n] = size(phi_x);
grad_abs = zeros(m,n);
for i = 1:m
    for j = 1:n
        grad_abs(i,j) = sqrt(phi_x(i,j)^2+phi_y(i,j)^2);
    end
end
[Pf , Pb] = apperaenceModel(I,phi,sigma);
% 
% [P_Mf,P_Mb] = appearence(I,phi);

P_x = HeavisideFunction(phi,sigma).*Pf + (1-HeavisideFunction(phi,sigma)).*Pb; %P（xi|phi,p,yi）出现0子项
P = Pf - Pb;
D = Dirac(phi,sigma);
T1 = D.*P./P_x; %公式12前面那一项  重点！！！！！！！！ 


% 散度计算
Div = div(phi_x./grad_abs,phi_y./grad_abs); 
% laplacian算子
s = [0 1 0;1 -4 1;0 1 0];
L = imfilter(phi,s); %Laplacian 算子
E_phi = T1 + 1/Delta^2*(L - Div); %公式12


function f = div(nx,ny)
    [nxx,junk]=gradient(nx);  
    [junk,nyy]=gradient(ny);
    f=nxx+nyy;

function f = Dirac(x, sigma)
f=(1/2/sigma)*(1+cos(pi*x/sigma));
b = (x<=sigma) & (x>=-sigma);
f = f.*b;

function [Pf, Pb] = apperaenceModel(I,phi,sigma) %求 Pf Pb
[m, n] = size(phi);
Pf = zeros(m,n);
Pb = zeros(m,n);
[nf, nb] = eta(phi,sigma);
[P_Mf, P_Mb]= appearence(I,phi);
Pf = P_Mf./(nf*P_Mf+nb*P_Mb);
Pb = P_Mb./(nf*P_Mf+nb*P_Mb); %计算 Pf 及 Pb


function [P_Mf, P_Mb] = appearence(I,phi) 
% 计算P（y|Mj）概率
[m_I, n_I] = size(I);
Ib_x = zeros(m_I*n_I);
Ib_y = zeros(m_I*n_I);
If_x = zeros(m_I*n_I);
If_y = zeros(m_I*n_I);
k=1;
l=1;
P_Mf = zeros(m_I,n_I);
P_Mb = zeros(m_I,n_I);
for i = 1:1:m_I
    for j = 1:1:n_I
        temp = phi(i,j);
        if temp>0
            Ib_x(k) = i;
            Ib_y(k) = j;
            k = k+1;
        else
            If_x(l) = i;
            If_y(l) = j;
            l = l+1;
        end
    end
end
Ib_x = Ib_x(Ib_x>0);
Ib_y = Ib_y(Ib_y>0);
If_x = If_x(If_x>0);
If_y = If_y(If_y>0);
Ib = [Ib_x Ib_y];
If = [If_x If_y];  % Ib If 为属于背景及前静的像素点 
P_b = size(Ib,1)/(size(Ib,1)+size(If,1));
P_f = size(If,1)/(size(Ib,1)+size(If,1));% P_b P_f 为 P（Mb）P(Mf)先验概率;



hist_Ib = zeros(1,256);
hist_If = zeros(1,256); %创建背景 前景 概率分布（灰度直方图）
% 
for i = 1:size(Ib,1)
    temp = I(Ib(i,1),Ib(i,2));
    hist_Ib(temp) = hist_Ib(temp)+1;
end
for i = 1:size(If,1)
    temp = I(If(i,1),If(i,2));
    hist_If(temp) = hist_If(temp)+1;
end %计算各灰度的密度

sum_f = sum(hist_If);
sum_b = sum(hist_Ib);
% P_Mf = hist_Ib(y)/sum_b;
% P_Mb = hist_If(y)/sum_f;
for i = 1:m_I  %计算 P（yi Mb）及 P（yi Mf）
    for j = 1:n_I
        P_Mb(i,j) = hist_Ib(I(i,j))/sum_b;
        P_Mf(i,j) = hist_If(I(i,j))/sum_f;
    end
end 
P_Mb = P_Mb / P_b;
P_Mf = P_Mf / P_f;


function [nf nb] = eta(phi,sigma) 
    f = HeavisideFunction(phi,sigma);
    nf = sum(sum(f));
    neg_f = 1 - f;
    nb = sum(sum(neg_f));
function g = NeumannBoundCond(f)
% Make a function satisfy Neumann boundary condition
[nrow,ncol] = size(f);
g = f;
g([1 nrow],[1 ncol]) = g([3 nrow-2],[3 ncol-2]);  
g([1 nrow],2:end-1) = g([3 nrow-2],2:end-1);          
g(2:end-1,[1 ncol]) = g(2:end-1,[3 ncol-2]);  


function f = HeavisideFunction(phi,sigma) %Heaviside函数 
[m n] = size(phi);
f = zeros(m,n);
for i = 1:m
    for j = 1:n
        temp = phi(i,j);
        if temp > sigma
            f(i,j) = 1;
        elseif temp < -1*sigma
            f(i,j) = 0;
        else
            f(i,j) = 0.5*(1+phi(i,j)/sigma+1/pi*sin(pi*phi(i,j)/sigma));
        end
    end
end

function f = div(nx,ny)
    [nxx,junk]=gradient(nx);  
    [junk,nyy]=gradient(ny);
    f=nxx+nyy;
