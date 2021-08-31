function loc_s = getSpecificLocations(cat_s, loc)
    % Getting a specific state's catastrophe locations
    global BORROWING
    if BORROWING
        loc_s = struct('j', loc.j(1:cat_s.j), ...
                       'h', loc.h(1:cat_s.h), ...
                       'p', loc.p(1:cat_s.p));
    else
        loc_s = [];
    end
end
