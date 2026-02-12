function Value = std_PDF_Residue(Zeta,varargin)
y=imag(Zeta);

if y>0
    Value=0;
else
    if y<0
        Value=2i*sqrt(pi)*exp(-Zeta^2);
    else
        Value=1i*sqrt(pi)*exp(-Zeta^2);
    end
end

if nargin==2
    Diff_Order=varargin{1};
    Value=Value*(-1)^Diff_Order*std_hermite(Diff_Order,Zeta);
end

end

