function visualize_tallthin(datasetfname1, datasetfname2, printflag, ...
                            colorflag, outdir)
% visualize(datasetfname1, datasetfname2, printflag, colorflag, outdir)
%
% Visualizes the information on Nystrom approximation stored in the file
% datasetfname1 and randEig sketches and Nystrom extensions stored in datasetfname2. If printflag 
% is true, saves the generated graphs to pdfs
% with appropriate names. If colorflag is true, uses color, otherwise uses
% grayscale. If outdir is present and printflag is true, the 
% pdfs are saved to this directory.
%
% datasetfname1 is assumed to contain the information on all Nystrom
% extensions other than the one using the tall thin approximation of 
% leverage scores
%
% datasetfname2 is assumed to contain the information on all
% randEig sketching mehtods
%
% Graphs generated:
%  - spectral, frobenius, and trace norm errors of the sketching methods 
%   for both the fixed and nonfixed-rank approximants

%% setup plotting
load(datasetfname1) % so basename can be defined

if printflag
    make_it_tight = true;
    subplot = @(m,n,p) subtightplot (m, n, p, [0.075 0.05], [0.045 0.03], [0.1 0.04]);
    if ~make_it_tight
        clear subplot;  
    end 

    fontsize = 15;
    width = 6.2;
    height = 12.4;
    basename = fullfile(outdir, savedata.in.datasetbasename);
    printfig = @(figname) printcf([basename figname '.pdf'], fontsize, width, height);
    printpanel = @(panelname, figname) printcf([basename figname '-' panelname '.pdf'], fontsize, width, width);
    
    timingwidth = 6.2;
    timingheight = 6.2;
    printtiming = @(figname) printcf([basename figname '.pdf'], fontsize, timingwidth, timingheight);
end


%% calculate the mean errors and timings

% Load the Nystrom info
load(datasetfname1);

for lidx = 1:length(savedata.in.lvals)
  % simple
    means = mean(savedata.simpleData(lidx).specerr,2);
    nonfixed_simple_specerr(lidx) = means(1);
    fixed_simple_specerr(lidx) = means(2);
    
    means = mean(savedata.simpleData(lidx).froerr,2);
    nonfixed_simple_froerr(lidx) = means(1);
    fixed_simple_froerr(lidx) = means(2);
    
    means = mean(savedata.simpleData(lidx).trerr,2);
    nonfixed_simple_trerr(lidx) = means(1);
    fixed_simple_trerr(lidx) = means(2);
    
    means = mean(savedata.simpleData(lidx).timings, 2);
    nonfixed_simple_timing(lidx) = means(1);
    fixed_simple_timing(lidx) = means(2);
    
    % srft
    means = mean(savedata.srftData(lidx).specerr,2);
    nonfixed_srft_specerr(lidx) = means(1);
    fixed_srft_specerr(lidx) = means(2);
    
    means = mean(savedata.srftData(lidx).froerr,2);
    nonfixed_srft_froerr(lidx) = means(1);
    fixed_srft_froerr(lidx) = means(2);
    
    means = mean(savedata.srftData(lidx).trerr,2);
    nonfixed_srft_trerr(lidx) = means(1);
    fixed_srft_trerr(lidx) = means(2);

    means = mean(savedata.srftData(lidx).timings, 2);
    nonfixed_srft_timing(lidx) = means(1);
    fixed_srft_timing(lidx) = means(2);
    
    % gaussian
    means = mean(savedata.gaussianData(lidx).specerr,2);
    nonfixed_gaussian_specerr(lidx) = means(1);
    fixed_gaussian_specerr(lidx) = means(2);
    
    means = mean(savedata.gaussianData(lidx).froerr,2);
    nonfixed_gaussian_froerr(lidx) = means(1);
    fixed_gaussian_froerr(lidx) = means(2);
    
    means = mean(savedata.gaussianData(lidx).trerr,2);
    nonfixed_gaussian_trerr(lidx) = means(1);
    fixed_gaussian_trerr(lidx) = means(2);
    
    means = mean(savedata.gaussianData(lidx).timings, 2);
    nonfixed_gaussian_timing(lidx) = means(1);
    fixed_gaussian_timing(lidx) = means(2);
    
    % levscores
    means = mean(savedata.levscoreData(lidx).specerr,2);
    nonfixed_levscore_specerr(lidx) = means(1);
    fixed_levscore_specerr(lidx) = means(2);
    
    means = mean(savedata.levscoreData(lidx).froerr,2);
    nonfixed_levscore_froerr(lidx) = means(1);
    fixed_levscore_froerr(lidx) = means(2);
    
    means = mean(savedata.levscoreData(lidx).trerr,2);
    nonfixed_levscore_trerr(lidx) = means(1);
    fixed_levscore_trerr(lidx) = means(2);
    
    means = mean(savedata.levscoreData(lidx).timings, 2);
    nonfixed_levscore_timing(lidx) = means(1);
    fixed_levscore_timing(lidx) = means(2);
end

% Load the randEig sketch and Nystrom extension info
load(datasetfname2);

for lidx = 1:length(savedata.in.lvals)
    % simple
    means = mean(savedata.eig_simpleData(lidx).specerr,2);
    eig_nonfixed_simple_specerr(lidx) = means(1);
    eig_fixed_simple_specerr(lidx) = means(2);
    
    means = mean(savedata.eig_simpleData(lidx).froerr,2);
    eig_nonfixed_simple_froerr(lidx) = means(1);
    eig_fixed_simple_froerr(lidx) = means(2);
    
    means = mean(savedata.eig_simpleData(lidx).trerr,2);
    eig_nonfixed_simple_trerr(lidx) = means(1);
    eig_fixed_simple_trerr(lidx) = means(2);
    
    means = mean(savedata.eig_simpleData(lidx).timings, 2);
    eig_nonfixed_simple_timing(lidx) = means(1);
    eig_fixed_simple_timing(lidx) = means(2);
    
    % srft
    means = mean(savedata.eig_srftData(lidx).specerr,2);
    eig_nonfixed_srft_specerr(lidx) = means(1);
    eig_fixed_srft_specerr(lidx) = means(2);
    
    means = mean(savedata.eig_srftData(lidx).froerr,2);
    eig_nonfixed_srft_froerr(lidx) = means(1);
    eig_fixed_srft_froerr(lidx) = means(2);
    
    means = mean(savedata.eig_srftData(lidx).trerr,2);
    eig_nonfixed_srft_trerr(lidx) = means(1);
    eig_fixed_srft_trerr(lidx) = means(2);

    means = mean(savedata.eig_srftData(lidx).timings, 2);
    eig_nonfixed_srft_timing(lidx) = means(1);
    eig_fixed_srft_timing(lidx) = means(2);
    
    % gaussian
    means = mean(savedata.eig_gaussianData(lidx).specerr,2);
    eig_nonfixed_gaussian_specerr(lidx) = means(1);
    eig_fixed_gaussian_specerr(lidx) = means(2);
    
    means = mean(savedata.eig_gaussianData(lidx).froerr,2);
    eig_nonfixed_gaussian_froerr(lidx) = means(1);
    eig_fixed_gaussian_froerr(lidx) = means(2);
    
    means = mean(savedata.eig_gaussianData(lidx).trerr,2);
    eig_nonfixed_gaussian_trerr(lidx) = means(1);
    eig_fixed_gaussian_trerr(lidx) = means(2);
    
    means = mean(savedata.eig_gaussianData(lidx).timings, 2);
    eig_nonfixed_gaussian_timing(lidx) = means(1);
    eig_fixed_gaussian_timing(lidx) = means(2);
    
    % levscores
    means = mean(savedata.eig_levscoreData(lidx).specerr,2);
    eig_nonfixed_levscore_specerr(lidx) = means(1);
    eig_fixed_levscore_specerr(lidx) = means(2);
    
    means = mean(savedata.eig_levscoreData(lidx).froerr,2);
    eig_nonfixed_levscore_froerr(lidx) = means(1);
    eig_fixed_levscore_froerr(lidx) = means(2);
    
    means = mean(savedata.eig_levscoreData(lidx).trerr,2);
    eig_nonfixed_levscore_trerr(lidx) = means(1);
    eig_fixed_levscore_trerr(lidx) = means(2);
    
    means = mean(savedata.eig_levscoreData(lidx).timings, 2);
    eig_nonfixed_levscore_timing(lidx) = means(1);
    eig_fixed_levscore_timing(lidx) = means(2);   
end


for lidx = 1:length(savedata.in.lvals)
    % simple
    means = mean(savedata.eignys_simpleData(lidx).specerr,2);
    eignys_nonfixed_simple_specerr(lidx) = means(1);
    eignys_fixed_simple_specerr(lidx) = means(2);
    
    means = mean(savedata.eignys_simpleData(lidx).froerr,2);
    eignys_nonfixed_simple_froerr(lidx) = means(1);
    eignys_fixed_simple_froerr(lidx) = means(2);
    
    means = mean(savedata.eignys_simpleData(lidx).trerr,2);
    eignys_nonfixed_simple_trerr(lidx) = means(1);
    eignys_fixed_simple_trerr(lidx) = means(2);
    
    means = mean(savedata.eignys_simpleData(lidx).timings, 2);
    eignys_nonfixed_simple_timing(lidx) = means(1);
    eignys_fixed_simple_timing(lidx) = means(2);
    
    % srft
    means = mean(savedata.eignys_srftData(lidx).specerr,2);
    eignys_nonfixed_srft_specerr(lidx) = means(1);
    eignys_fixed_srft_specerr(lidx) = means(2);
    
    means = mean(savedata.eignys_srftData(lidx).froerr,2);
    eignys_nonfixed_srft_froerr(lidx) = means(1);
    eignys_fixed_srft_froerr(lidx) = means(2);
    
    means = mean(savedata.eignys_srftData(lidx).trerr,2);
    eignys_nonfixed_srft_trerr(lidx) = means(1);
    eignys_fixed_srft_trerr(lidx) = means(2);

    means = mean(savedata.eignys_srftData(lidx).timings, 2);
    eignys_nonfixed_srft_timing(lidx) = means(1);
    eignys_fixed_srft_timing(lidx) = means(2);
    
    % gaussian
    means = mean(savedata.eignys_gaussianData(lidx).specerr,2);
    eignys_nonfixed_gaussian_specerr(lidx) = means(1);
    eignys_fixed_gaussian_specerr(lidx) = means(2);
    
    means = mean(savedata.eignys_gaussianData(lidx).froerr,2);
    eignys_nonfixed_gaussian_froerr(lidx) = means(1);
    eignys_fixed_gaussian_froerr(lidx) = means(2);
    
    means = mean(savedata.eignys_gaussianData(lidx).trerr,2);
    eignys_nonfixed_gaussian_trerr(lidx) = means(1);
    eignys_fixed_gaussian_trerr(lidx) = means(2);
    
    means = mean(savedata.eignys_gaussianData(lidx).timings, 2);
    eignys_nonfixed_gaussian_timing(lidx) = means(1);
    eignys_fixed_gaussian_timing(lidx) = means(2);
    
    % levscores
    means = mean(savedata.eignys_levscoreData(lidx).specerr,2);
    eignys_nonfixed_levscore_specerr(lidx) = means(1);
    eignys_fixed_levscore_specerr(lidx) = means(2);
    
    means = mean(savedata.eignys_levscoreData(lidx).froerr,2);
    eignys_nonfixed_levscore_froerr(lidx) = means(1);
    eignys_fixed_levscore_froerr(lidx) = means(2);
    
    means = mean(savedata.eignys_levscoreData(lidx).trerr,2);
    eignys_nonfixed_levscore_trerr(lidx) = means(1);
    eignys_fixed_levscore_trerr(lidx) = means(2);
    
    means = mean(savedata.eignys_levscoreData(lidx).timings, 2);
    eignys_nonfixed_levscore_timing(lidx) = means(1);
    eignys_fixed_levscore_timing(lidx) = means(2);   
end

%% Plot the errors of the simple, gaussian, srft, lev score methods
% display the uniform and QR leverage sampling errors for calibration
lw = 1.5;
ms = 10;

if colorflag
    simple_style = 'ks-';
    simple_color = 'k';
    eig_simple_style = 'ks:';
    eig_simple_color = 'k';
    eignys_simple_style = 'ks--';
    eignys_simple_color = 'k';
    
    srft_style = 'bo-';
    srft_color = 'b';
    eig_srft_style = 'bv:';
    eig_srft_color = 'b';
    eignys_srft_style = 'bv--';
    eignys_srft_color = 'b';
    
    gaussian_style = 'go-';
    gaussian_color = 'g';
    eig_gaussian_style = 'go:';
    eig_gaussian_color = 'g'l
    eignys_gaussian_style = 'go--';
    eignys_gaussian_color = 'g';
    
    levscore_style = 'rd-';
    levscore_color = 'r';
    eig_levscore_style = 'rd:';
    eig_levscore_color = 'r';
    eignys_levscore_style = 'rd--';
    eignys_levscore_color = 'r';
else
    simple_style = 's-';
    simple_color = [.1 .1 .1];
    eig_simple_style = 's:';
    eig_simple_color = [.1 .1 .1];
    eignys_simple_style = 's--';
    eignys_simple_color = [.1 .1 .1];
    
    gaussian_style = 'o-';
    gaussian_color = [.2 .2 .2];
    eig_gaussian_style = 'o:';
    eig_gaussian_color = [.2 .2 .2];
    eignys_gaussian_style = 'o--';
    eignys_gaussian_color = [.2 .2 .2];
    
    srft_style = 'v-';
    srft_color = [.3 .3 .3];
    eig_srft_style = 'v:';
    eig_srft_color = [.3 .3 .3];
    eignys_srft_style = 'v--';
    eignys_srft_color = [.3 .3 .3];
    
    levscore_style = 'd-';
    levscore_color = [.4 .4 .4];
    eig_levscore_style = 'd:';
    eig_levscore_color = [.4 .4 .4];
    eignys_levscore_style = 'd--';
    eignys_levscore_color = [.4 .4 .4];
end
legendloc = 'Northeast';

exactbase_struct = struct(...
    'seriesnames', {{'unif', 'unif_eig', 'unif_eignys', 'gaussian', ...
       'gaussian_eig','gaussian_eignys', 'srft', 'srft_eig', ...
       'srft_eignys'}}, ...
    'xlabel', 'l', ...
    'clipaxis', true, ...
    'legendloc', legendloc, ...
    'x', savedata.in.lvals, ...
    'styles', {{simple_style, eig_simple_style, eignys_simple_style, ...
        gaussian_style, eig_gaussian_style, eignys_gaussian_style, ...
        srft_style, eig_srft_style, eignys_srft_style}}, ...
    'lw', {{lw, lw, lw, lw, lw, lw, lw, lw, lw}}, ...
    'ms', {{ms, ms, ms, ms, ms, ms, ms, ms, ms}}, ...
    'colors', {{simple_color, eig_simple_color, eignys_simple_color, ...
        gaussian_color, eig_gaussian_color, eignys_gaussian_color, ...
        srft_color, eig_srft_color, eignys_srft_color}}, ...
    'mcolors', {{simple_color, eig_simple_color, eignys_simple_color, ...
        gaussian_color, eig_gaussian_color, eignys_gaussian_color, ...
        srft_color, eig_srft_color, eignys_srft_color}});
    
%% Exact method, nonfixed rank plots

exactmethods_nonfixedrank_spec_struct = exactbase_struct;
exactmethods_nonfixedrank_frob_struct = exactbase_struct;
exactmethods_nonfixedrank_trace_struct = exactbase_struct;

exactmethods_nonfixedrank_spec_struct.series = { ...
    nonfixed_simple_specerr/savedata.optspecerr, ...
    eig_nonfixed_simple_specerr/savedata.optspecerr, ...
    eignys_nonfixed_simple_specerr/savedata.optspecerr, ...
    nonfixed_gaussian_specerr/savedata.optspecerr, ...
    eig_nonfixed_gaussian_specerr/savedata.optspecerr, ...
    eignys_nonfixed_gaussian_specerr/savedata.optspecerr, ...
    nonfixed_srft_specerr/savedata.optspecerr, ...
    eig_nonfixed_srft_specerr/savedata.optspecerr, ...
    eignys_nonfixed_srft_specerr/savedata.optspecerr};

exactmethods_nonfixedrank_spec_struct.plotname = 'Relative spectral error';

exactmethods_nonfixedrank_frob_struct.series = { ...
    nonfixed_simple_froerr/savedata.optfroerr, ...
    eig_nonfixed_simple_froerr/savedata.optfroerr, ...
    eignys_nonfixed_simple_froerr/savedata.optfroerr, ...
    nonfixed_gaussian_froerr/savedata.optfroerr, ...
    eig_nonfixed_gaussian_froerr/savedata.optfroerr, ...
    eignys_nonfixed_gaussian_froerr/savedata.optfroerr, ...
    nonfixed_srft_froerr/savedata.optfroerr, ...
    eig_nonfixed_srft_froerr/savedata.optfroerr, ...
    eignys_nonfixed_srft_froerr/savedata.optfroerr};

exactmethods_nonfixedrank_frob_struct.plotname = 'Relative Frobenius error';

exactmethods_nonfixedrank_trace_struct.series = { ...
    nonfixed_simple_trerr/savedata.opttrerr, ...
    eig_nonfixed_simple_trerr/savedata.opttrerr, ...
    eignys_nonfixed_simple_trerr/savedata.opttrerr, ...
    nonfixed_gaussian_trerr/savedata.opttrerr, ...
    eig_nonfixed_gaussian_trerr/savedata.opttrerr, ...
    eignys_nonfixed_gaussian_trerr/savedata.opttrerr, ...
    nonfixed_srft_trerr/savedata.opttrerr, ...
    eig_nonfixed_srft_trerr/savedata.opttrerr, ...
    eignys_nonfixed_srft_trerr/savedata.opttrerr};

exactmethods_nonfixedrank_trace_struct.plotname = 'Relative trace error';

figure();
subplot(3,1,1);
plotpanel(exactmethods_nonfixedrank_spec_struct);
subplot(3,1,2);
plotpanel(exactmethods_nonfixedrank_frob_struct);
subplot(3,1,3);
plotpanel(exactmethods_nonfixedrank_trace_struct);

% if storing plots, store the plot with all three norms at once, then
% separate plots for the individual norms
if printflag
    printfig('exact-methods-nonfixed-rank-errors-with-eig');

	figure();
	plotpanel(exactmethods_nonfixedrank_spec_struct);
	printpanel('specerr', 'exact-methods-nonfixed-rank-errors-with-eig');

	figure();
	plotpanel(exactmethods_nonfixedrank_frob_struct);
	printpanel('froberr', 'exact-methods-nonfixed-rank-errors-with-eig');

	figure();
	plotpanel(exactmethods_nonfixedrank_trace_struct);
	printpanel('trerr', 'exact-methods-nonfixed-rank-errors-with-eig');
end

%% close all figures
if printflag
    close all
end

end