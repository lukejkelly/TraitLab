
if ~exist('noutput','var')
    noutput=getcladetimes(data,output,model.prior);
end

clade=DefineClades;
figure;
pstr={'r--','b'};
for t=1:17, 
    hist(noutput.cladetimes(t,:),20); 
    hold on;
    title(model.prior.clade{t}.name); 
    xy=axis;
    if ~isempty(clade{t}.rootrange)
        for b=1:2
            if isfinite(clade{t}.rootrange(b))
                plot([clade{t}.rootrange(b),clade{t}.rootrange(b)],xy(3:4),pstr{b});
            end
        end
    end
    clade{t}, pause; 
    hold off;
end


