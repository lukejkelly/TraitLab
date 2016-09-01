function [state,truestate,model]=initSTATE(mcmc,model,data,MCMCINITTREE)

global TRUSTART EXPSTART OLDSTART LEAF

disp(sprintf('\n***Preparing data for analysis'));
treedata = data.content;
if data.initial.is_column_mask
    disp(sprintf('\n**Applying column (trait) mask'));
    treedata = ObserveData(treedata,[],data.initial.column_mask,0);
end
if data.initial.masking
    disp(sprintf('\n**Applying row (taxon) mask'));
    treedata = ObserveData(treedata,data.initial.mask,[],0); 
end
if (~data.initial.masking & ~data.initial.is_column_mask) | data.initial.lost>0
    disp(sprintf('\n**Removing traits displayed by less than %d taxa',data.initial.lost+1));
    treedata = ObserveData(treedata,[],[],data.initial.lost);
end
if model.prior.isclade
    [bad,verybad,verybadn]=badclades(model.prior.clade,treedata.language);
    %keyboard;
    %dcl=[];
    if ~isempty(verybadn)
        %dcl=verybad; %#ok<NASGU>
        disp(sprintf('\n*******************************************************************'));
        disp(sprintf('The following clades have no taxa in the data:'));
        disp([sprintf('Names: '), sprintf('%s ',verybad{:})]);
        disp([sprintf('Index: '), sprintf('%d ',verybadn)]);
        disp(sprintf('*******************************************************************\n'));
    end
    somemissing=setdiff([bad([bad.index]>0).index],verybadn);
    if ~isempty(somemissing)
        disp(sprintf('\n*******************************************************************'));
        disp(sprintf('The following clades have taxa absent in the data.\nFor each clade, you can either completely ignore this issue, or delete the absent taxa from the clade, depending on what is appropriate.\nYou may prefer to stop this run and add the problematic clades to the list of clades to ignore.\nYou should also follow this course of action if your choices lead to unexplained errors.'));
        for k=somemissing
            disp([sprintf('CLADE: %d %s MISSING TAXA:',k,model.prior.clade{k}.name),sprintf(' %s',bad(k).langs{:})]);
            %remlan=setdiff(model.prior.clade{k}.language,bad(k).langs);
            %disp([sprintf('CLADE: %d %s REMAINS LANGS:',k,model.prior.clade{k}.name),sprintf(' %s',remlan{:})]);
            disp(model.prior.clade{k});
            cst=input('delete taxa/ignore? 1/0 ');
            if cst==1
                model.prior.clade{k}.language=setdiff(model.prior.clade{k}.language,bad(k).langs);
            %elseif cst==2
            %    dcl=[dcl,{model.prior.clade{k}.name}];
            end
        end
        %k=1;
        %while k<=size(model.prior.clade,2)
        %    if any(strcmp(model.prior.clade{k}.name,dcl))
        %       model.prior.clade={model.prior.clade{1:(k-1)},model.prior.clade{(k+1):end}};
        %        k=k-1;
        %    end
        %    k=k+1
        %end
        disp(sprintf('*******************************************************************\n'));
    end
end
disp(sprintf('***Data preparation complete: %g distinct traits in %g taxa.',treedata.L,treedata.NS));

disp(sprintf('\n***Preparing MCMC start state'));

switch mcmc.initial.setup
    case TRUSTART
        if data.initial.masking         
            itree=BuildSubTree(data.true.state.tree,data.initial.mask,data.content.language);
        else
            itree= data.true.state.tree ;
        end
        state=makestate(model.prior,mcmc.initial.mu,mcmc.initial.lambda,mcmc.initial.p,mcmc.initial.rho,mcmc.initial.kappa,treedata,itree,data.true.cat, mcmc.initial.beta); % LUKE 05/10/2013
        disp('MCMC initial tree is the synthetic true');
    case EXPSTART
        if ~model.prior.isclade
            itree=ExpTree(treedata.NS,mcmc.initial.theta);
            [itree([itree.type]==LEAF).Name] = deal(treedata.language{:});
            state=makestate(model.prior,mcmc.initial.mu,mcmc.initial.lambda,mcmc.initial.p,mcmc.initial.rho,mcmc.initial.kappa,treedata,itree, [], mcmc.initial.beta); % LUKE 05/10/2013
            disp('MCMC initial tree simulated from ExpTree');
        else
            disp('Constructing MCMC initial tree using RandCladeTree');
            topage=-inf;
            for cld=[model.prior.clade{:}]
                topage=max([topage,min(cld.rootrange),min(cld.adamrange)]);
            end
            if isempty(model.prior.rootmax)
                rootmax = 2*sum(1./(2:(treedata.NS-1)))/mcmc.initial.theta;
                %according to me this is 2*mean depth of a exp-branching tree
                %with NS leaves
                if rootmax<topage
                    rootmax=topage*2;
                    disp(sprintf('rootmax estimated from branching rate is more recent than oldest constraint: setting rootmax to %g',rootmax));
                end
            else
                rootmax = model.prior.rootmax;
                if rootmax<topage
                    rootmax=topage*2;
                    disp(sprintf('Imposed rootmax is more recent than oldest constraint: setting rootmax to %g',rootmax'));
                end
            end
            itree = RandCladeTree(model.prior.clade,treedata.language,rootmax,mcmc.initial.rho);
            itree = mergetreedata(itree,treedata.array,treedata.language);
            ok = checkcladetreematch(itree,model.prior.clade);
            if ~ok
                error('tree and data incompatible')
            end
            state=makestate(model.prior,mcmc.initial.mu,mcmc.initial.lambda,mcmc.initial.p,mcmc.initial.rho,mcmc.initial.kappa,treedata,itree,[], mcmc.initial.beta); % LUKE
        end
    case OLDSTART
        if isempty(MCMCINITTREE)
            disp('empty MCMCINITTREE passed to initSTATE() for MCMC initial setup');
        else
            if data.initial.masking
                itree=BuildSubTree(MCMCINITTREE.tree,data.initial.mask,data.content.language);          
                disp(sprintf('\nData in analysis contains %g distinct traits in %g taxa.\n',treedata.L,treedata.NS));
            else
                itree=MCMCINITTREE.tree;
            end
            ok = checkcladetreematch(itree,model.prior.clade);
            if ~ok
                error('tree and data incompatible')
            end
            state = makestate(model.prior,MCMCINITTREE.mu,MCMCINITTREE.lambda,MCMCINITTREE.p,MCMCINITTREE.rho,MCMCINITTREE.kappa,treedata,itree,MCMCINITTREE.cat, MCMCINITTREE.beta); % LUKE 05/10/2013
            disp('MCMC initial tree loaded from output file ');
        end
    otherwise
        error('no starting tree type recognised in initSTATE')
end
truestate=data.true.state;
if ~isempty(data.true.state)
    disp(sprintf('\n***Preparing true state for comparison purposes'));
    if data.initial.masking
        ttree=BuildSubTree(data.true.state.tree,data.initial.mask,data.content.language);
    else
        ttree=data.true.state.tree;
    end
    truestate=makestate(model.prior,data.true.state.mu,data.true.state.lambda,data.true.state.p,data.true.state.rho, data.true.state.kappa,treedata,ttree,data.true.state.cat, data.true.state.beta); % Luke 04/05/2015
end
