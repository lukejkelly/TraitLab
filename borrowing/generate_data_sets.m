GlobalSwitches; GlobalValues;

list_L = 4:2:12;
root_time = 1e3;
list_lambda = [0.05, 0.1, 0.2];
mu = 5e-4;
beta = 0; % 5e-4;

for L = list_L
    for lambda = list_lambda
        actually_sample_data(L, root_time, lambda, mu, beta);
    end
end

make_submit_file(list_L, list_lambda);

%%

function actually_sample_data(L, root_time, lambda, mu, beta)

    global ROOT ANST

    % Parameters of SDLT process
    kappa = NaN;
    tr = [lambda; mu; beta; kappa];

    output_stem = sprintf('L%d_r%.0e_l%.0e_m%.0e_b%.0e', L, root_time, tr(1:3));
    output_nex = sprintf('%s.nex', output_stem);

    % Generate tree
    THETA = 1;
    s = ExpTree(L, THETA);
    old_root_time = s([s.type] == ROOT).time;
    for j = find([s.type] == ANST | [s.type] == ROOT)
       s(j).time = s(j).time * root_time / old_root_time;
    end

    % Read events and generate data
    [tEvents, rl] = stype2Events(s);
    D = simBorCatDeath(tEvents, tr);

    % Removing empty site-patterns
    D = D(sum(D, 2) > 0, :);

    % Adding data to leaves
    for l = 1:L
      s(rl(l)).dat = D(:, L + 1 - l)';
    end

    % % Make some clades
    % clades = {struct('name', 'root_clade', ...
    %                  'rootrange', [750, 1250], ...
    %                  'adamrange', [], ...
    %                  'language', {{'1', '2', '3', '4', '5'}})};

    % Write to Nexus file
    % sFile = stype2nexus(s, '', 'BOTH', '', clades);
    sFile = stype2nexus(s, 'Estimate coupling distributions', 'DATA', '');
    fid = fopen(fullfile('data', output_nex), 'w');
    fprintf(fid, sFile);
    fclose(fid);

    make_par_file(output_stem, mu, beta, 2.5e5, 1e1)
end

function make_par_file(FILE_NAME, MU, BETA, RUN_LENGTH, SAMPLE_INTERVAL)

    fid = fopen(fullfile('pars', [FILE_NAME, '.par']), 'w');

    fprintf(fid, '%% FULL PATH OF DATA FILE INCLUDE .NEX EXTENSION\n');
    fprintf(fid, 'Data_file_name = data/%s.nex\n', FILE_NAME);
    fprintf(fid, '\n');
    fprintf(fid, '%% ONE OF THE FOLLOWING THREE OPTIONS MUST BE ONE, THE OTHERS 0\n');
    fprintf(fid, 'Start_from_rand_tree = 1\n');
    fprintf(fid, 'Start_from_tree_in_output = 0\n');
    fprintf(fid, 'Start_from_true_tree = 0\n');
    fprintf(fid, '\n');
    fprintf(fid, '%% VALUE OF THETA IGNORED UNLESS Start_from_rand_tree == 1\n');
    fprintf(fid, 'Theta = 0.010000\n');
    fprintf(fid, '\n');
    fprintf(fid, '%% NEXT TWO FIELDS IGNORED UNLESS Start_from_tree_in_output == 1\n');
    fprintf(fid, '%% FULL PATH OF OLD OUTPUT FILE INCLUDE .NEX EXTENSION\n');
    fprintf(fid, 'Tree_file_name =\n');
    fprintf(fid, 'Use_tree = 0\n');
    fprintf(fid, '\n');
    fprintf(fid, 'Omit_taxa = 0\n');
    fprintf(fid, '%% LIST IS IGNORED UNLESS Omit_taxa == 1 CAN USE MATLAB VECTOR NOTATION\n');
    fprintf(fid, 'Omit_taxa_list =\n');
    fprintf(fid, '\n');
    fprintf(fid, 'Omit_traits = 0\n');
    fprintf(fid, '%% LIST IS IGNORED UNLESS Omit_traits == 1 CAN USE MATLAB VECTOR NOTATION\n');
    fprintf(fid, 'Omit_trait_list =\n');
    fprintf(fid, '\n');
    fprintf(fid, '%% ONE OF THE FOLLOWING TWO OPTIONS MUST BE 1, THE OTHER 0\n');
    fprintf(fid, 'Yule_prior_on_tree = 0\n');
    fprintf(fid, 'Flat_prior_on_tree = 1\n');
    fprintf(fid, '%% ONE OF THE FOLLOWING TWO OPTIONS MUST BE 1, THE OTHER 0\n');
    fprintf(fid, 'Uniform_prior_on_tree_topologies = 1\n');
    fprintf(fid, 'Uniform_prior_on_labelled_histories = 0\n');
    fprintf(fid, '%% FOLLOWING IS IGNORED UNLESS Flat_prior_on_tree == 1\n');
    fprintf(fid, 'Max_root_age = 10000\n');
    fprintf(fid, '\n');
    fprintf(fid, 'Vary_topology = 1\n');
    fprintf(fid, 'Account_rare_traits = 0\n');
    fprintf(fid, '\n');
    fprintf(fid, 'Impose_clades = 0\n');
    fprintf(fid, '%% LIST IS IGNORED UNLESS Impose_clades == 1 CAN USE MATLAB VECTOR NOTATION\n');
    fprintf(fid, 'Omit_clade_list =\n');
    fprintf(fid, 'Omit_clade_ages_list =\n');
    fprintf(fid, '\n');
    fprintf(fid, 'Vary_loss_rate = 0\n');
    fprintf(fid, '%% FOLLOWING IS IGNORED WHEN Random_initial_loss_rate == 1\n');
    fprintf(fid, 'Initial_loss_rate = %g\n', LossRate(MU));
    fprintf(fid, 'Random_initial_loss_rate = 0\n');
    fprintf(fid, '\n');
    fprintf(fid, 'Account_for_lateral_transfer = 0\n');
    fprintf(fid, '%% FOLLOWING IS IGNORED WHEN Account_for_lateral_transfer == 0\n');
    fprintf(fid, 'Vary_borrowing_rate = 1\n');
    fprintf(fid, 'Random_initial_borrowing_rate = 1\n');
    fprintf(fid, '%% NEXT LINE IS IGNORED WHEN Random_initial_borrowing_rate == 1\n');
    fprintf(fid, 'Initial_borrowing_rate = %g\n', BETA);
    fprintf(fid, '\n');
    fprintf(fid, 'Include_catastrophes = 0\n');
    fprintf(fid, '%% NEXT 6 LINES ARE IGNORED WHEN Include_catastrophes = 0\n');
    fprintf(fid, '%% FOLLOWING IS IGNORED WHEN Random_initial_cat_death_prob = 1\n');
    fprintf(fid, 'Initial_cat_death_prob = 0.5\n');
    fprintf(fid, 'Random_initial_cat_death_prob = 0\n');
    fprintf(fid, '%% FOLLOWING IS IGNORED WHEN Random_initial_cat_rate = 1\n');
    fprintf(fid, 'Initial_cat_rate = 0.00015\n');
    fprintf(fid, 'Random_initial_cat_rate = 0\n');
    fprintf(fid, '\n');
    fprintf(fid, 'Model_missing = 0\n');
    fprintf(fid, '\n');
    fprintf(fid, 'Run_length = %g\n', RUN_LENGTH);
    fprintf(fid, 'Sample_interval = %g\n', SAMPLE_INTERVAL);
    fprintf(fid, '\n');
    fprintf(fid, 'Seed_random_numbers = 0\n');
    fprintf(fid, '%% FOLLOWING IS IGNORED UNLESS Seed_random_numbers == 1\n');
    fprintf(fid, 'With_seed = 0\n');
    fprintf(fid, '\n');
    fprintf(fid, '%% OUTPUT FILE NAME OMITTING PATH AND ANY EXTENSIONS\n');
    fprintf(fid, 'Output_file_name = %s\n', FILE_NAME);
    fprintf(fid, '%% FULL PATH FOR DIRECTORY TO OUTPUT FILES\n');
    fprintf(fid, 'Output_path_name = output/\n');
    fprintf(fid, '\n');
    fprintf(fid, '%% Gaps are treated as missing data. To change this, edit GlobalValues.m.\n');

    fclose(fid);

end

function make_submit_file(Ls, LAMBDAs)

    fname = sprintf('submit-job-%s.pbs', datestr(datetime, 'yy-mm-dd'));
    fid = fopen(fullfile('~', 'Desktop', fname), 'w');

    fprintf(fid, '#!/bin/bash\n');
    fprintf(fid, '');
    fprintf(fid, 'for L in');
    fprintf(fid, ' %d', Ls);
    fprintf(fid, '; do\n');

    fprintf(fid, '    for LAMBDA in');
    fprintf(fid, ' %g', LAMBDAs);
    fprintf(fid, '; do\n');

    fprintf(fid, '        qsub -N ${L}-${LAMBDA} -v L=${L},LAMBDA=${LAMBDA} job-sd.pbs\n');
    fprintf(fid, '    done\n');
    fprintf(fid, 'done\n');

    fclose(fid);

end
