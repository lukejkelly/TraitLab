function makeJobFile(dest_dir)

    fid = fopen(fullfile(dest_dir, 'job.pbs'), 'w');

    fprintf(fid, '#!/bin/bash\n\n');

    fprintf(fid, '#PBS -M kelly@ceremade.dauphine.fr\n');
    fprintf(fid, '#PBS -m e\n');
    fprintf(fid, '#PBS -l nodes=clust2:ppn=1\n');
    fprintf(fid, '#PBS -o log/kelly.o$PBS_JOBID\n');
    fprintf(fid, '#PBS -j oe\n\n');

    fprintf(fid, 'sleep $[ ($RANDOM %% 60) + 1 ]s\n\n');

    fprintf(fid, 'PAR_FILE=$(LC_NUMERIC="en_GB.UTF-8" \\\n');
    fprintf(fid, '           printf "L%%d_r%%.0e_l%%.0e_m%%.0e_b%%.0e_n%%.0e_s%%.0e" \\\n');
    fprintf(fid, '                  "$L" "$ROOT_TIME" "$LAMBDA" "$MU" "$BETA" "$RUN_LENGTH" \\\n');
    fprintf(fid, '                  "$SAMPLE_INTERVAL")\n\n');

    fprintf(fid, 'cd ~/TraitLabSDLT-coupled\n');
    fprintf(fid, '/usr/local/bin/matlab -nodesktop -nodisplay \\\n');
    fprintf(fid, '    -r "batchTraitLabCoupled(''%s/${PAR_FILE}.par'', ${PBS_ARRAYID}); exit"\n', ...
            fullfile(dest_dir, 'pars'));

end
