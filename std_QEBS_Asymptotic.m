function [Value,Error] = std_QEBS_Asymptotic(Zeta,Diff_Order,Error_Abs,Error_Rel)
Residue=std_PDF_Residue(Zeta,Diff_Order);

zz=1./2./Zeta.^2;
N_MaxLim=1000;   

Value_Add_Base=(-1).^Diff_Order./Zeta.^(1+Diff_Order);
Value_Base=Value_Add_Base.*std_pochhammer(1,Diff_Order);  
Value_RefBase=abs(Value_Base);
Value_0=1e-16;    
Flag_Ref=true;       
Value_Record=Value_RefBase;  

for N=1 : N_MaxLim
    Value_Add_Base=Value_Add_Base.*(2*N-1).*zz;
    Value_Add=Value_Add_Base.*std_pochhammer(2*N+1,Diff_Order);
    Value_Add_ABS=abs(Value_Add);
    %¡¾Ê§Ð§ÖÐ¶Ï¡¿
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
    if Value_Record<Value_Add_ABS
        break;
    end
    
    
    Value_Record=Value_Add_ABS;         
    Value_Base=Value_Base+Value_Add;   
    if abs(Value_Base)>Value_0                
        Value_RefBase=abs(Value_Base);
        Flag_Ref=true;
    else
        Flag_Ref=false;
    end
    
end

Error=Value_Add_ABS;

Value=Residue-Value_Base;

end
