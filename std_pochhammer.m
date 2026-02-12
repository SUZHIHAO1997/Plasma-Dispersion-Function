function Value = std_pochhammer(ValueBase,Levels)
if Levels==0
    Value=1;
else
    Array=ValueBase:ValueBase+Levels-1;
    Value=prod(Array);
end

end

