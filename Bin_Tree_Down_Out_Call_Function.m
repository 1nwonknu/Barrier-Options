function f = Bin_Tree_Down_Out_Call_Function( r,sig,S,K,T,L,n)
    u= exp (sig*sqrt(T/n)) ;
    d=1/u;
    p = (exp(r*T/n)-d)/(u-d);
                                        %% log:=ln
    m = ceil(log(K/(S*d^n))/log(u/d));  %% at least m up steps are required
                                        %% for the call option to be in the money
    h= floor(log(L/(S*d^n))/log (u/d));
    j=m;
    b = ( nchoosek (n,j)-    nchoosek (n,(n-2*h+j)) ) * (p^j)*(1-p)^(n-j); 
    D = S* u^j* d^(n-j);
    C= b * (D-K);
    for j = m+1: n
        b= b * (p * (n-j+1)/((1-p)*j));
        D= D * (u/d);
        C = C + b* (D-K);
    end 
    f=C* exp(-r*T);
end