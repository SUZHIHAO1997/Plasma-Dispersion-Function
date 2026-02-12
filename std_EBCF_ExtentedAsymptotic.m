function [Value,Error] = std_EBCF_ExtentedAsymptotic(Zeta,Diff_Order,Error_Abs,Error_Rel)
Residue=std_PDF_Residue(Zeta,Diff_Order);

if mod(Diff_Order,2)==0
    [KummerU_Value,Error]=SUB_kummerU(1+Diff_Order/2,3/2,-Zeta.^2,Error_Abs,Error_Rel);
    KummerU_Value=(-1).^(Diff_Order/2).*prod(1:Diff_Order).*Zeta.*KummerU_Value;
else
    [KummerU_Value,Error]=SUB_kummerU((Diff_Order+1)/2,1/2,-Zeta.^2,Error_Abs,Error_Rel);
    KummerU_Value=(-1).^((Diff_Order+1)/2).*prod(1:Diff_Order).*KummerU_Value;
end

Value=KummerU_Value+Residue;

end

function [Value,Error] =SUB_kummerU(a,c,z,Error_Abs,Error_Rel)
b=1+a-c;
An1=1+6*(z-a.*b)./(a+1)./(b+1).*(1+2*z./(a+2)./(b+2))+6*a.*b./(a+2)./(b+2);
Bn1=1+6*z./(a+1)./(b+1)+12*z.^2./(a+1)./(a+2)./(b+1)./(b+2);
An2=1+2*(z-a.*b)./(a+1)./(b+1);
Bn2=1+2*z./(a+1)./(b+1);
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
    F1=(3*n.^2+(a+b-6).*n+2-a.*b-2*(a+b)+2*(2*n-3).*z)./(2*n-3)./(n+a)./(n+b).*(2*n-1);
    F2=-(3*n.^2-(a+b+6).*n+2-a.*b-2*(2*n-1).*z)./(n+a)./(n+b);
    F3=(2*n-1).*(n-a-2).*(n-b-2)./(2*n-3)./(n+a)./(n+b);
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
Value=Value.*z.^(-a);

end


