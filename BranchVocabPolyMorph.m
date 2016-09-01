function [nextwords,nextL,q]=BranchVocabPolyMorph(prevwords,p,min_wpm)

%function [nextwords,nextL,q]=BranchVocabPolyMorph(prevwords,p)
%simulate cladgenic loss conditional on at least one surviving

prevL=length(prevwords);

if min_wpm==1
    %q is the conditional probability for n words to
    %survive (each survives wp p) conditioned on at least one surviving
    for n=1:prevL
        q(n)=nchoosek(prevL,n)*p^n*(1-p)^(prevL-n);
    end
    
    %prob non survive unconditioned thinning is 1-(1-p)^(n words)
    q=q./(1-(1-p)^prevL);
    
    %how many survive?
    c=disample(q);
    
    %conditional on c words surviving choose the c survivors at
    %random from the set of prevwords
    v=randperm(prevL);
    nextwords=prevwords(v(1:c));
else
    %thin with survival probability p no requirement at least one survives
    nextwords=prevwords(rand(1,prevL)<p);
end

nextL=length(nextwords);