%% Adapted from Lyuu, 2004
%% if L equals S, the value of the Down and In Call Option equals that of
%% an equivalent vanilla call option with the same strike
r=0.0043;                            %% risk-free rate per year
sig=0.11;
S=160;
K=150;                              %% Strike
L=145;                              %% Lower barrier
T=0.5;                              %% lifetime in years
n = floor(T*365);                   %% Steps
u= exp (sig*sqrt(T/n)) ;
d=1/u;
p = (exp(r*T/n)-d)/(u-d);
dr = exp(-r*T);                     %% Discount Rate := r * T
                                    %% log:=ln
m = ceil(log(K/(S*d^n))/log(u/d));  %% at least m up steps are required
                                    %% for the call option to be in the money
h= floor(log(L/(S*d^n))/log (u/d));                                    

tic;
j=m;
b = nchoosek (n,(n-2*h+j)) * (p^j)*(1-p)^(n-j); 
D = S* u^j* d^(n-j);
C= b * (D-K);
for j = m+1: n
    b= b * (p * (n-j+1)/((1-p)*j));
    D= D * (u/d);
    C = C + b* (D-K);
end
C*dr
toc;