function [true,content]=truthN(prior,initial,true)

%function [true,content]=truthN(prior,initial,true)
%
%prior       prior model                  (model.prior)
%content     synthetic data               (data.content)
%initial     data initialisation          (data.initial)
%true        true state initialisation    (data.true)
         
global OLDTRE NEWTRE OFF ON

disp(sprintf('Synthesizing data using parameter values'));
disp(sprintf('mu = %g [trait deaths/year]\ntheta = %g [taxon splits/year]\nmean num traits = %g [traits]\n',true.mu,true.theta,true.vocabsize));

%(A) SIMULATE OR LOAD THE TREE ON WHICH THE COGNATE DATA WILL BE SIMULATED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch initial.synth
case OLDTRE
    if ~isempty(initial.tree)
        stree = initial.tree;
    else
        stree=nexus2stype(initial.treefile); 
        disp(sprintf(['Simulating trait data on tree loaded from ',initial.treefile,'.\n']));
    end
    % convert the names on the leaves to numbers
    % leaves = find([stree.type]==LEAF);
    %if ~prior.isclade
    %   langnames=num2cell(strjust(num2str((1:length(leaves))'),'left'),2);
    %   [stree(leaves).Name]=deal(langnames{:});
    %   disp('Leaves of OLDTRE renamed by node number (no clading applied)');
    %else
    %   disp('Leaves of OLDTRE retain old names (so clades work)');
    %end
case NEWTRE
    if true.NS==2
        stree=DebugTree(true.NS,true.theta);
        disp('True tree simulated using DebugTree.m has two leaves.');
    else
        stree=ExpTree(true.NS,true.theta); %TODO: synthetic style might depend on model.prior.type
        disp('True tree simulated using ExpTree.m has exp-branch lengths');
    end
end

%(B) SIMULATE SYNTHETIC COGNATE DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%(B.1) simulate cognates on tree
%TODO - change this so only use the PolyMorph() function
%       and make true.br=0 for no borrowing etc. 
%       Reserve Mutations() and GTRmut() for testing.
% TODO missing data handling is messy
if initial.polymorph, polymorph_string='polymorphic '; else polymorph_string=''; end
if initial.borrow, borrow_string=' with borrowing,'; else borrow_string=' with no borrowing,'; end
if true.p<1, cladagenic_string=' with cladagenic loss (p<1),'; else cladagenic_string=' with no cladagenic loss (p=1),'; end
if true.mu>0, anagenic_string=' with anagenic loss (mu>0),'; else anagenic_string=' with no anagenic loss (mu=0),'; end
if initial.cats, cat_string=' with catastrophic events'; else cat_string=' with no catastrophic events'; end
if initial.nmeaningclass > 1, nmc_string = sprintf(' for %1g observation classes',initial.nmeaningclass); else nmc_string = ''; end
if initial.borrow && initial.localborrow
    locborr_string = sprintf('\nBorrowing is local: taxa with a common ancestor in last %1g years may borrow',initial.maxdist);
else
    locborr_string = '';
end
if initial.missing, missing_string=' and with missing data.'; else missing_string=' and with no missing data.'; end
if initial.rhbranchvar > 0, rhbranch_str = sprintf('\nHeterogenity in anagenic loss rate across branches with std dev %g.',initial.rhbranchvar); 
else rhbranch_str = sprintf('\nNo rate heterogenity across branches.');end
if initial.rhclassvar > 0, rhclass_str = sprintf('\nHeterogenity in anagenic loss rate across observation classes with std dev %g.',initial.rhbranchvar); 
else rhclass_str = sprintf('\nNo rate heterogenity across observation classes.');end
    
    content = pop('content');
       

% Modified by RJR 2011-03-02    
if initial.polymorph || ( ( initial.borrow || initial.cats ) && ( initial.rhbranchvar || initial.rhclassvar ) )
    disp('Simulating synthetic data using PolyMorph()');
    if initial.cats
        [true.wordset,content.language,content.cognate,content.NS,content.L,content.L_vals,stree,true.cat]=PolyMorph(true.br,true.mu,true.p,true.lambda,stree,initial.polymorph,initial.nmeaningclass,initial.localborrow,initial.maxdist,initial.rhbranchvar,initial.rhclassvar,1,true.rho,true.kappa,initial.knowcats);
    else
        [true.wordset,content.language,content.cognate,content.NS,content.L,content.L_vals]=PolyMorph(true.br,true.mu,true.p,true.lambda,stree,initial.polymorph,initial.nmeaningclass,initial.localborrow,initial.maxdist,initial.rhbranchvar,initial.rhclassvar,0);
    end

    
elseif (initial.borrow || initial.cats)
    disp('Simulating synthetic data using SimBo()');
    [true.wordset, content.language, content.cognate, content.NS, content.L, true.cat, stree]=SimBo(true.br,true.mu,true.p,true.lambda,true.rho,true.kappa,true.nu,stree);
    
else
    disp('Simulating synthetic data using MutationsRH()');
    [true.wordset,content.language,content.cognate,content.NS,content.L]=MutationsRH(true.mu,true.p,true.lambda,stree,initial.rhbranchvar);
end

% switch (initial.borrow || initial.cats) 
%     case OFF
%         switch initial.polymorph
%             case OFF
%            %     if true.lambda/true.mu>2000, error('required mean number of traits very large'); end
%                 disp('Simulating synthetic data using MutationsRH()');
%                 [true.wordset,content.language,content.cognate,content.NS,content.L]=MutationsRH(true.mu,true.p,true.lambda,stree,initial.rhbranchvar);
%             case ON
%                 disp('Simulating synthetic data using PolyMorph()');
%                 [true.wordset,content.language,content.cognate,content.NS,content.L,content.L_vals]=PolyMorph(true.br,true.mu,true.p,true.lambda,stree,initial.polymorph,initial.nmeaningclass,initial.localborrow,initial.maxdist,initial.rhbranchvar,initial.rhclassvar);
%         end
%     case ON
%         %disp('Simulating synthetic data using PolyMorph()');
%         %[true.wordset,content.language,content.cognate,content.NS,content.L]=PolyMorph(true.br,true.mu,true.p,true.lambda,stree,initial.polymorph,initial.nmeaningclass,initial.localborrow,initial.maxdist,initial.rhbranchvar,initial.rhclassvar);
%         disp('Simulating synthetic data using SimBo()');
%         [true.wordset, content.language, content.cognate, content.NS, content.L, true.cat, stree]=SimBo(true.br,true.mu,true.p,true.lambda,true.rho,true.kappa,true.nu,stree);
% end
disp([polymorph_string,'traits simulated',nmc_string,borrow_string,cladagenic_string,anagenic_string,cat_string,locborr_string,cat_string,missing_string,rhbranch_str,rhclass_str]);





%(B.2) convert the sets of cognate labels into a NSxL data array
%if initial.borrow || ~initial.polymorph %Used SimBo
    content.array = words2array(true.wordset,stree,content.L);
%else % Used PolyMorph
%    content.array = words2array(true.wordset,stree,content.L,true.missingwords);
%end
    
%(B.2a) data go missing %RJR 29/05/08
if initial.missing
    if initial.polymorph % Used PolyMorph
        content.array=SimMis(content.array,stree,content.L,content.L_vals);
    else % Used SimBo or MutationsRH; content.L_vals is not defined
        content.array=SimMis(content.array,stree,content.L);
    end
end


%(B.3) drop cognates that are present in initial.lost languages
content = ObserveData(content,[],[],initial.lost);


%(C) INSTALL DATA IN TREE, RETURN COMPLETE STATE
true.state = makestate(prior,true.mu,true.lambda,true.p,true.rho,true.kappa,content,stree,true.cat, true.beta); % LUKE 05/10/2013

disp(sprintf('\nSimulated data contains %g distinct traits in %g taxa.\n',content.L,content.NS));
disp(sprintf('simulated on a tree of depth %g years\n',true.state.tree(true.state.root).time));

save outSD;

