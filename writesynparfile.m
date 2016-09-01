function ok = writesynparfile(parfilename,fsu,numtree)

global NEWTRE
ok = 1;
%
fid = fopen(parfilename,'w');
if fid > 0
    % write header
    fprintf(fid,'%% Parameters for accompanying .nex file created by TraitLab Synthesize GUI on %s\n\n',datestr(clock));
    % fsu.VOCABSIZE         = VS    ;
    %
    fprintf(fid,'Mean_number_of_traits = %g\n',fsu.VOCABSIZE);
    % fsu.LOSSRATE          = LR    ;
    %
    fprintf(fid,'Loss_Rate = %g\n',fsu.LOSSRATE);
    % fsu.LOSSRATEBRANCHVAR = RHB   ;
    %
    if fsu.LOSSRATEBRANCHVAR > 0
        fprintf(fid,'Rate_het_across_branches = 1\n');
        fprintf(fid,'Branch_ratehet_stddev = %g\n',fsu.LOSSRATEBRANCHVAR);
    else
        fprintf(fid,'Rate_Het_across_branches = 0\n');
    end
    % fsu.POLYMORPH         = PM    ;
    %
    if fsu.POLYMORPH
        fprintf(fid,'Make_polymorphic_data = 1\n');
        % fsu.NMEANINGCLASS     = NMC   ;
        %
        fprintf(fid,'Observation_classes = %g\n',fsu.NMEANINGCLASS);
        % fsu.LOSSRATECLASSVAR  = RHC   ;
        %
        if fsu.LOSSRATECLASSVAR > 0
            fprintf(fid,'Rate_het_across_classes = 1\n');
            fprintf(fid,'Class_ratehet_stddev = %g\n',fsu.LOSSRATECLASSVAR);
        else
            fprintf(fid,'Rate_Het_across_classes = 0\n');
        end
    else
        fprintf(fid,'Make_polymorphic_data = 0\n');
    end
    % fsu.LOST              = LT    ;
    %
    % LT = 1;   %delete cognates present in LOST or less of the languages
    fprintf(fid,'Remove_rare_traits = %g\n',fsu.LOST);
    % fsu.BORROW            = BO    ;
    %
    fprintf(fid,'Allow_borrowing = %g\n',fsu.BORROW);
    if fsu.BORROW
        % fsu.BORROWFRAC        = BF    ;
        %
        fprintf(fid,'Borrowing_rate = %g\n',fsu.BORROWFRAC);
        % fsu.LOCALBORROW       = LB    ;
        %
        fprintf(fid,'Local_borrowing = %g\n',fsu.LOCALBORROW);
        if fsu.LOCALBORROW
            % fsu.MAXDIST           = MD    ;
            fprintf(fid,'Max_borrow_dist = %g\n',fsu.MAXDIST);
        end
    end
    % %SYNTHSTYLE {OLDTRE NEWTRE} if DATASOURCE=BUILD
    % fsu.SYNTHSTYLE        = ST    ;
    %
    if fsu.SYNTHSTYLE == NEWTRE
        fprintf(fid,'Synthesize_on_random_tree = 1\n');
        % fsu.THETA             = TH    ;
        %
        fprintf(fid,'Branching_rate = %g\n',fsu.THETA);
        % fsu.NUMSEQ            = NS    ;
        %
        fprintf(fid,'Taxa = %g\n',fsu.NUMSEQ);
    else
        fprintf(fid,'Synthesize_on_tree_from_file = 1\n');
        % fsu.SYNTHTREFILE      = STF   ;
        fprintf(fid,'Tree_file = %s\n',fsu.SYNTHTREFILE);
        % fsu.SYNTHTRE          = STR   ;
        fprintf(fid,'Use_tree = %g\n',numtree);
    end
    fclose(fid);
else
    disp(sprintf('problem opening %s while trying to write .par file',parfilename))
    ok = 0;
end