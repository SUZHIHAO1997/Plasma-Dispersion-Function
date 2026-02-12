function [Value,Error] = std_PDF(Zeta,Diff_Order,Error_Abs,Error_Rel)
%   Value£º     function value
%   Error£º     estimation of error
%   Zeta£º      independent variable
%   Diff_Order£ºorder of derivative
%   Error_Abs£º absolute error tolerance
%   Error_Re£º  relative error tolerance


if nargin<4
    Error_Rel=1e-15;
    if nargin<3
        Error_Abs=1e-15;
        if nargin<2
            Diff_Order=0;
        end
    end
end

x=abs(real(Zeta));
y=abs(imag(Zeta));
if y>=1
    [Value,Error]=std_EBCF_ExtentedAsymptotic(Zeta,Diff_Order,Error_Abs,Error_Rel);
elseif x>8
    [Value,Error]=std_QEBS_Asymptotic(Zeta,Diff_Order,Error_Abs,Error_Rel);
elseif x<1 && y<1
    [Value,Error]=std_QEBS_Taylor_Dawson(Zeta,Diff_Order,Error_Abs,Error_Rel);
else
    [Value,Error]=std_EBCF_ExtentedTaylor(Zeta,Diff_Order,Error_Abs,Error_Rel);
end

end


