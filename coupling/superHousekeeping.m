function s2 = superHousekeeping(s1, s2)
    % Update sibling information
    global ANST ROOT;

    nodes = find(ismember([s1.type], [ANST, ROOT]));
    if ~all(nodes == find(ismember([s2.type], [ANST, ROOT])))
        error('Node sets do not match');
    end

    for i = nodes
        % If children are the same but in the opposite order in s2 then flip
        if isequal(s1(i).child, fliplr(s2(i).child))
            s2(i).child = fliplr(s2(i).child);
            [s2(s2(i).child).sibling] = deal(s2(fliplr(s2(i).child)).sibling);
        end
    end
end
