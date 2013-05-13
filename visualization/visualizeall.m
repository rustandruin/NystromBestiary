% Generates all the plots for the paper
% and truncates them if you have perl and pdfcrop installed
 
% change this to point to the directory/ies containing the outputs from
% experiments you want to visualize

dirlist = {'outputs'};
printflag = true; % generate pdfs of the plots
outdir = 'plots'; % location to store the pdfs 

for diridx = 1:length(dirlist)
    curdir = dirlist{diridx};
    datafnames = dir(fullfile(curdir, '*.mat'));
    
    for fileidx = 1:length(datafnames)
        curfile = fullfile(curdir, datafnames(fileidx).name);
        try 
            % the try is here because the tall thin datasets will fail in 
            % visualize because they don't have all the expected data
            visualize(curfile, printflag, outdir);
        catch
        end
    end
end

% point to the particular datasets you want to compare the timing performance of the
% tallthin method on
visualize_tallthin('outputs/SNPSrank5.mat', 'outputs/tallthinSNPSrank5.mat', true, outdir);
visualize_tallthin('outputs/Proteinrank10.mat', 'outputs/tallthinProteinrank10.mat', true, outdir);

wd = cd(outdir);
[~, result] = system('./truncateplots.pl'); % crops all output plots
cd(wd);

disp(result)
