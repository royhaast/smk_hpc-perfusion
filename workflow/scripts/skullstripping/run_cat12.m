function run_cat12(in,out)

addpath('./spm12','./spm12/toolbox/cat12');

load('./batch_cat12.mat')
matlabbatch{1}.spm.tools.cat.estwrite.data = {in};
out_mat = [out '/batch_cat12.mat'];
save(out_mat,'matlabbatch');

nrun = 1; % enter the number of runs here
jobfile = {out_mat};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
end

spm('defaults', 'FMRI');
cat12('expert');

spm_jobman('run', jobs, inputs{:});

end