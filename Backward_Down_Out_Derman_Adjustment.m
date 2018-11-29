function f = Backward_Down_Out_Derman_Adjustment_v8(r,sig,S,K,T,L,n)
    u= exp (sig*sqrt(T/n)) ;
    d=1/u;
    p = (exp(r*T/n)-d)/(u-d);
    D = S* u^0* d^n;   
    A = D * (u/d).^(n:-1:0);    
    P = max ((A-K),0);                
    P2= P;
    A2=A;                      
    intpol=0;
    for i= n:-1:1  
        if i>1
            E = (A>L);         
            [~,I] = max (E(1:end-1) -E(2:end));
            k = I(1,1);         
            A2 = u * A2(2:end);                
            E = (A2>L);                  
            [~,I] = max (E(1:end-1) -E(2:end)); %% 
            j = I(1,1);        
            if A(1,k) > A2(1,j)               
                                
                P (1,k) = (A(1,k)-L)/(A(1,k)-A(1,k+1))*P2(1,k) + ... 
                            (L-A(1,k+1))/(A(1,k)-A(1,k+1))*P2(1,k+1);  
                intpol=1;            
            end     
        end
        A = u * A(2:end);
        P =(p * P (1:end-1)+(1-p) * P(2:end));     
        P = P.*(A> L);   
        P2 =(p * P2 (1:end-1)+(1-p) * P2(2:end));  
        P2 = P2.*(A> L);                                   
        if intpol==1
            P(1,k) = P2(1,k);
            intpol=0;
        end        
    end
    f = P * exp(-r*T);
end