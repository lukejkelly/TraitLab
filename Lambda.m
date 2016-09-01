function lambda_sample=Lambda(state)


if state.L==0
    lambda_sample=0;
else
    [llkd,X]=LogLkd(state);
    lambda_sample=randG(state.L,X);
end

