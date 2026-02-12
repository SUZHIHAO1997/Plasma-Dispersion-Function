function [Value,Error] = std_QEBS_Taylor_Dawson(Zeta,Diff_Order,Error_Abs,Error_Rel)
H_term=(-1).^Diff_Order.*std_hermite(Diff_Order,Zeta).*1i*sqrt(pi).*exp(-Zeta.^2);

if mod(Diff_Order,2)==1
    zz=Zeta.^2;
    
    N_MaxLim=1000; 
    
    Value_Add_Base=(-2).^((Diff_Order+1)./2)./prod(1:2:Diff_Order);
    Value_Base=Value_Add_Base.*std_pochhammer(1,Diff_Order);            
    Value_RefBase=abs(Value_Base);           
    Value_0=1e-16;               
    Flag_Ref=true;                
    

    
    for N=1 : N_MaxLim
        Value_Add_Base=Value_Add_Base.*(-2*zz)./(2*N+Diff_Order);   
        Value_Add=std_pochhammer(2.*N+1,Diff_Order).*Value_Add_Base; 
        Value_Add_ABS=abs(Value_Add);
        if isnan(Value_Add)
            Value=NaN;
            Error=Inf;
            return;
        end
        if Flag_Ref         
            if Error_Rel>Value_Add_ABS/Value_RefBase
                break;
            end
        else                  
            if Error_Abs>Value_Add_ABS
                break;
            end
        end
        
        Value_Base=Value_Base+Value_Add;   
        if abs(Value_Base)>Value_0               
            Value_RefBase=abs(Value_Base);
            Flag_Ref=true;
        else
            Flag_Ref=false;
        end
        
    end
else
    zz=Zeta.^2;
    N_MaxLim=1000; 
    Value_Add_Base=(-2).^(Diff_Order./2+1).*Zeta./prod(1:2:Diff_Order+1);
    Value_Base=Value_Add_Base.*std_pochhammer(2,Diff_Order);   
    Value_RefBase=abs(Value_Base); 
    Value_0=1e-16;   
    Flag_Ref=true;    
    
    for N=1 : N_MaxLim
        Value_Add_Base=Value_Add_Base.*(-2*zz)./(2*N+1+Diff_Order); 
        Value_Add=std_pochhammer(2.*N+2,Diff_Order).*Value_Add_Base;  
        Value_Add_ABS=abs(Value_Add);
        if isnan(Value_Add)
            Value=NaN;
            Error=Inf;
            return;
        end
        if Flag_Ref         
            if Error_Rel>Value_Add_ABS/Value_RefBase
                break;
            end
        else                 
            if Error_Abs>Value_Add_ABS
                break;
            end
        end
        
        Value_Base=Value_Base+Value_Add;   
        if abs(Value_Base)>Value_0                
            Value_RefBase=abs(Value_Base);
            Flag_Ref=true;
        else
            Flag_Ref=false;
        end
        
    end
    
    
end

Value=H_term+Value_Base;
Error=Value_Add_ABS;

end
