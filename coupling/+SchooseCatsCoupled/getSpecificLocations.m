function loc_s = getSpecificLocations(cat_s, loc)
    % Getting a specific state's catastrophe locations
    global BORROWING
    if BORROWING
        loc_s = struct('i', loc.i(1:cat_s.i), ...
                       'j', loc.j(1:cat_s.j), ...
                       'k', loc.k(1:cat_s.k));
    else
        loc_s = [];
    end
end
