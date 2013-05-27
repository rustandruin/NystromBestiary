% Generates all the plots for the paper
% and truncates them if you have perl and pdfcrop installed

thesisq = false; % produce plots appropriate for a thesis
printflag = true; % generate pdfs of the plots
colorflag = false; % generate color plots
outdir = 'plots'; % location to store the pdfs 

% basic SPSD sketches
dirlist = {'outputs'};
for diridx = 1:length(dirlist)
    curdir = dirlist{diridx};
    datafnames = dir(fullfile(curdir, '*.mat'));
    
    for fileidx = 1:length(datafnames)
        curfile = fullfile(curdir, datafnames(fileidx).name);
        try 
            % the try is here because the tall thin datasets will fail in 
            % visualize because they don't have all the expected data
            if thesisq
                visualizeforthesis(curfile, printflag, colorflag, outdir);
            else
                visualize(curfile, printflag, colorflag, outdir);
            end
        catch
        end
    end
end

% compare SPSD sketches to projection-based approximants
dirlist = {'densified'};
for diridx = 1:length(dirlist)
    curdir = dirlist{diridx};
    datafnames = dir(fullfile(curdir, '*.mat'));
    
    for fileidx = 1:length(datafnames)
        curfile = fullfile(curdir, datafnames(fileidx).name);
        if thesisq
            visualize_randeigforthesis(curfile, curfile, printflag, colorflag, outdir);
        else
            visualize_randeig(curfile, curfile, printflag, colorflag, outdir);
        end
    end
end

if ~thesisq
    % point to the particular datasets you want to compare the timing performance of the
    % tallthin method on
    visualize_tallthin('outputs/SNPSrank5.mat', 'outputs/tallthinSNPSrank5.mat', true, false, outdir);
    visualize_tallthin('outputs/Proteinrank10.mat', 'outputs/tallthinProteinrank10.mat', true, false, outdir);
end

wd = cd(outdir);
[~, result] = system('./truncateplots.pl'); % crops all output plots
cd(wd);

disp(result)
