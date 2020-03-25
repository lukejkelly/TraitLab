function ok = writerunparfile(parfilename,fsu,numtree)

global EXPSTART OLDSTART TRUSTART FLAT YULE VARYMU GAP MIX OUT TOPO LABHIST
ok = 1;
%
fid = fopen(parfilename,'w');
if fid > 0
    % write header
    fprintf(fid,'%% Run parameters accompanying .nex and .txt file for TraitLab MCMC run on %s\n\n',datestr(clock));
    % fsu.VOCABSIZE         = VS    ;
    %
    % fsu.DATAFILE          = DFS   ;
    fprintf(fid,'%% FULL PATH OF DATA FILE INCLUDE .NEX EXTENSION\n');
    fprintf(fid,'Data_file_name = %s\n',fsu.DATAFILE);
    % fsu.DATASOURCE        = DS    ;


    % fsu.MCMCINITTREESTYLE = MI    ;
    fprintf(fid,'\n%% ONE OF THE FOLLOWING THREE OPTIONS MUST BE ONE, THE OTHERS 0\n');
    fprintf(fid,'Start_from_rand_tree = %g\n',fsu.MCMCINITTREESTYLE == EXPSTART);
    fprintf(fid,'Start_from_tree_in_output = %g\n',fsu.MCMCINITTREESTYLE ==OLDSTART);
    fprintf(fid,'Start_from_true_tree = %g\n',fsu.MCMCINITTREESTYLE ==TRUSTART);
    % fsu.MCMCINITTHETA     = IT    ;
    fprintf(fid,'\n%% VALUE OF THETA IGNORED UNLESS Start_from_rand_tree == 1\n');
    fprintf(fid,'Theta = %f\n',fsu.MCMCINITTHETA);

    % fsu.MCMCINITTREEFILE  = MIF   ;
    fprintf(fid,'\n%% NEXT TWO FIELDS IGNORED UNLESS Start_from_tree_in_output == 1\n');
    fprintf(fid,'%% FULL PATH OF OLD OUTPUT FILE INCLUDE .NEX EXTENSION\n');
    fprintf(fid,'Tree_file_name = %s\n',fsu.MCMCINITTREEFILE);
    fprintf(fid,'Use_tree = %g\n\n',numtree);

    % fsu.MASKING           = MK    ;
    fprintf(fid,'Omit_taxa = %g\n',fsu.MASKING);
    % fsu.DATAMASK          = DM    ;
    fprintf(fid,'%% LIST IS IGNORED UNLESS Omit_taxa == 1 CAN USE MATLAB VECTOR NOTATION\n');
    fmtstr = ['Omit_taxa_list =' repmat(' %g',1,length(fsu.DATAMASK)) '\n\n'];
    fprintf(fid,fmtstr,fsu.DATAMASK);
    % fsu.ISCOLMASK         = ICM   ;end
    fprintf(fid,'Omit_traits = %g\n',fsu.ISCOLMASK);
    % fsu.COLUMNMASK        = CM    ;
    fprintf(fid,'%% LIST IS IGNORED UNLESS Omit_traits == 1 CAN USE MATLAB VECTOR NOTATION\n');
    fmtstr = ['Omit_trait_list =' repmat(' %g',1,length(fsu.COLUMNMASK)) '\n'];
    fprintf(fid,fmtstr,fsu.COLUMNMASK);

    % fsu.TREEPRIOR     = TP   ;
    fprintf(fid,'\n%% ONE OF THE FOLLOWING TWO OPTIONS MUST BE 1, THE OTHER 0\n');
    fprintf(fid,'Yule_prior_on_tree = %g\n',fsu.TREEPRIOR == YULE);
    fprintf(fid,'Flat_prior_on_tree = %g\n',fsu.TREEPRIOR == FLAT);

    %GKN 18 Mar 2011 - added topology prior record
    % fsu.TOPOLOGYPRIOR         = TOPOLOGYPRIOR    ;
    fprintf(fid,'%% ONE OF THE FOLLOWING TWO OPTIONS MUST BE 1, THE OTHER 0\n');
    fprintf(fid,'Uniform_prior_on_tree_topologies = %g\n',fsu.TOPOLOGYPRIOR == TOPO);
    fprintf(fid,'Uniform_prior_on_labelled_histories = %g\n',fsu.TOPOLOGYPRIOR == LABHIST);

    % fsu.ROOTMAX           = RM    ;
    fprintf(fid,'%% FOLLOWING IS IGNORED UNLESS Flat_prior_on_tree == 1\n');
    fprintf(fid,'Max_root_age = %g\n\n',fsu.ROOTMAX);
    % fsu.MCMCVARYTOP       = VT    ;
    fprintf(fid,'Vary_topology = %g\n',fsu.MCMCVARYTOP);
    % fsu.LOSTONES          = LO    ;
    fprintf(fid,'Account_rare_traits = %g\n\n',fsu.LOSTONES);
    % fsu.ISCLADE           = IC    ;
    fprintf(fid,'Impose_clades = %g\n',fsu.ISCLADE);
    % fsu.CLADEMASK         = CLM   ;
    fprintf(fid,'%% LIST IS IGNORED UNLESS Impose_clades == 1 CAN USE MATLAB VECTOR NOTATION\n');
    fmtstr = ['Omit_clade_list =' repmat(' %g',1,length(fsu.CLADEMASK)) '\n'];
    fprintf(fid,fmtstr,fsu.CLADEMASK);
    fmtstr = ['Omit_clade_ages_list =' repmat(' %g',1,length(fsu.CLADEAGESMASK)) '\n\n'];
    fprintf(fid,fmtstr,fsu.CLADEAGESMASK);

    % VARYMU
    fprintf(fid,'Vary_loss_rate = %g\n',VARYMU);
    % fsu.LOSSRATE          = LR    ;
    fprintf(fid,'%% FOLLOWING IS IGNORED WHEN Random_initial_loss_rate == 1\n');
    fprintf(fid,'Initial_loss_rate = %g\n',fsu.LOSSRATE);
    % fsu.ISLRRANDOM
    fprintf(fid,'Random_initial_loss_rate = %g\n\n',fsu.ISLRRANDOM);

    % Lateral transfer LJK 23/3/20
    fprintf(fid, 'Account_for_lateral_transfer = %g\n', fsu.BORROWING);
    fprintf(fid, '%% FOLLOWING IS IGNORED WHEN Account_for_lateral_transfer == 0\n');
    fprintf(fid, 'Vary_borrowing_rate = %g\n', fsu.VARYBETA);
    fprintf(fid, 'Random_initial_borrowing_rate = %g\n', fsu.ISBETARANDOM);
    fprintf(fid, '%% NEXT LINE IS IGNORED WHEN Random_initial_borrowing_rate == 1\n');
    fprintf(fid, 'Initial_borrowing_rate = %g\n\n', fsu.MCMCINITBETA);

    fprintf(fid,'Include_catastrophes = %g\n',fsu.MCMCCAT);
    fprintf(fid,'%% NEXT 6 LINES ARE IGNORED WHEN Include_catastrophes = 0\n');
    fprintf(fid,'%% FOLLOWING IS IGNORED WHEN Random_initial_cat_death_prob = 1\n');
    fprintf(fid,'Initial_cat_death_prob = %g\n',fsu.MCMCINITKAPPA);
    fprintf(fid,'Random_initial_cat_death_prob = %g\n',fsu.RANDOMKAPPA);
    fprintf(fid,'%% FOLLOWING IS IGNORED WHEN Random_initial_cat_rate = 1\n');
    fprintf(fid,'Initial_cat_rate = %g\n',fsu.MCMCINITRHO);
    fprintf(fid,'Random_initial_cat_rate = %g\n\n',fsu.RANDOMRHO);

    fprintf(fid,'Model_missing = %g\n\n',fsu.MCMCMISS);

    % fsu.RUNLENGTH         = RL    ;
    fprintf(fid,'Run_length = %g\n',fsu.RUNLENGTH);
    % fsu.SUBSAMPLE         = SS    ;
    fprintf(fid,'Sample_interval = %g\n\n',fsu.SUBSAMPLE);

    % fsu.SEEDRAND          = SR    ;
    fprintf(fid,'Seed_random_numbers = %g\n',fsu.SEEDRAND);

    fprintf(fid,'%% FOLLOWING IS IGNORED UNLESS Seed_random_numbers == 1\n');
    % fsu.SEED              = SE    ;
    fprintf(fid,'With_seed = %g\n\n',fsu.SEED);

    % fsu.OUTFILE           = OF    ;
    fprintf(fid,'%% OUTPUT FILE NAME OMITTING PATH AND ANY EXTENSIONS\n');
    fprintf(fid,'Output_file_name = %s\n',fsu.OUTFILE);
    % fsu.OUTPATH           = OP    ;
    fprintf(fid,'%% FULL PATH FOR DIRECTORY TO OUTPUT FILES\n');
    fprintf(fid,'Output_path_name = %s\n\n',fsu.OUTPATH);

    %GAP
    if ~isempty(GAP)
        fprintf(fid,'%%Gaps are treated as ');
        if GAP==MIX, fprintf(fid,'missing data. '); elseif GAP==OUT, fprintf(fid,'absence of trait. '); end
        fprintf(fid,'To change this, edit GlobalValues.m.');
    end



    fclose(fid);
else
    disp(sprintf('problem opening %s while trying to write .par file',parfilename))
    ok = 0;
end
