function val = SatisfiesParachute(hkm,vf)

v = [438.5,487];
h = [6,7.98];

val = (hkm >= 6) && (vf >= 310) && (vf <= v(1));
if val
    return
elseif vf > v(1) && vf < v(2) && hkm >= interp1(v,h,vf,'linear')
    v2 = [476.4,v(2)];
    h2 = [16.73,h(2)];
    if vf > v2(1)
        val = hkm <= interp1(v2,h2,vf,'linear');
    else
        val = true;
    end
else
    val = false;
end

end