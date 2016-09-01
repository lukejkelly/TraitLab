function s=timeconvert(s,rt)

% s=timeconvert(s,rt)
% given a tree s with times at the nodes representing branch length
% timeconvert returns a tree with times representing time below
% the nominated root, rt.  Note that if rt is an internal node
% the times will be adjusted for the subtree below rt
% Note also that an adjustment will have to be made so the most recent
% node will have time zero (or at least some non-negative number)

s(rt).time = s(s(rt).parent).time - s(rt).time;
if ~isempty(s(rt).child)
    for i=1:length(s(rt).child)
        s=timeconvert(s,s(rt).child(i));
        if i>2
            disp(['Error in timeconvert - node ' num2str(rt) ' has >2 children'])
            keyboard
        end
    end
end
