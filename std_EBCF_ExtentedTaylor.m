function [Value,Error] = std_EBCF_ExtentedTaylor(Zeta,Diff_Order,Error_Abs,Error_Rel)
H_term=1i*sqrt(pi)*exp(-Zeta.^2).*(-1).^Diff_Order.*std_hermite(Diff_Order,Zeta);
if mod(Diff_Order,2)==0
    [KMValue,Error]=SUB_kummerM(1+Diff_Order/2,1.5,-Zeta^2,Error_Abs,Error_Rel);
    M_term=((-2)^(Diff_Order/2))*SUB_Double_Factorial(Diff_Order)*Zeta*KMValue;
elseif mod(Diff_Order,2)==1
    [KMValue,Error]=SUB_kummerM((1+Diff_Order)/2,0.5,-Zeta^2,Error_Abs,Error_Rel);
    M_term=((-2)^((Diff_Order-1)/2))*SUB_Double_Factorial(Diff_Order-1)*KMValue;
end

Value=H_term-2*M_term;

end

function [Value,Error] =SUB_kummerM(a,c,z,Error_Abs,Error_Rel)
An1=1+(a.*c+2*a-2*c).*z./2./c./(c+1)+(a-1).*(a-2).*z.^2./12./c./(c+1);
Bn1=1-(a+2).*z./2./(c+1)+(a+1).*(a+2).*z.^2./12./c./(c+1);
An2=1+(a-1).*z./2./c;
Bn2=1-(a+1).*z./2./c;
An3=1;
Bn3=1;
Value1=An1/Bn1;
Value2=An2/Bn2;
Value3=1;
Delta1=Value1-Value2;
Delta2=Value2-Value3;

Flag_Rel=true;    
Value_0=1e-16;   

for n=3:1000
    F1=1-(n-a-2)./2./(2*n-3)./(n+c-1).*z;
    F2=(2*(2*n-1).*(n-c-1)+(n+a).*z)./4./(2*n-1)./(2*n-3)./(n+c-2)./(n+c-1).*(n+a-1).*z;
    F3=(n+a-2).*(n+a-1).*(n-a-2)./8./(2*n-3).^2./(2*n-5)./(n+c-3)./(n+c-2)./(n+c-1).*z.^3;
    An=F1.*An1+F2.*An2+F3.*An3;
    Bn=F1.*Bn1+F2.*Bn2+F3.*Bn3;
    Delta=-((Bn2*F2+Bn3*F3)*Delta1+Bn3*F3*Delta2)/Bn;
    
    Value=An./Bn;
    Value_Abs=abs(Value);
    Error=abs(Delta);
    if n>2
        if Flag_Rel
            if Error_Rel>Error/Value_Abs
                break;
            end
        else
            if Error_Abs>Error
                break;
            end
        end
    end
    
    if Value_Abs>Value_0
        Flag_Rel=true;
    else
        Flag_Rel=false;
    end
    
    Delta2=Delta1;
    Delta1=Delta;
    An3=An2;
    Bn3=Bn2;
    An2=An1;
    Bn2=Bn1;
    An1=An;
    Bn1=Bn;
    
end
Error=abs(Delta);

end

function Value=SUB_Double_Factorial(n)
if n==0
    Value=1;
else
    Value=prod(2:2:n);
end

end
