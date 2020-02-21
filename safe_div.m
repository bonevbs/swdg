% safe_div function as in the spherical case. 
function v = safe_div(r,q,altv)
    v = r./q;
    v(isnan(v)) = altv;
end