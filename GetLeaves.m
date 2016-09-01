function y=GetLeaves(s,langs)

% function y=GetLeaves(s,langs)

GlobalSwitches;

y=[]; 
for k=1:size(langs,2) 
    fi=find(strcmp(langs(k),{s.Name}));
    %if length(fi)~=1, disp('One of the Languages was not found or was not unique'); end
    y=[y,fi];
end
