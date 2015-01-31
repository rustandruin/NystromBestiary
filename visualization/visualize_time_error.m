function visualize_time_error(datasetfname, printflag, colorflag, outdir)
% visualize(datasetfname, printflag, colorflag, outdir)
%
% Visualizes the information on Nystrom approximation stored in the file
% datasetfname in the form of time-error curves. If printflag is true, saves
% the generated graphs to pdfs % with appropriate names. If colorflag is true,
% uses colors, otherwise uses grayscale. If outdir is present and printflag is
% true, the pdfs are saved to this directory.
%
% Graphs generated:
%  - spectral, frobenius, and trace norm errors of the exact sampling
%  methods for both the fixed and nonfixed-rank approximants
%  - spectral, frobenius, and trace norm errors of the inexact leverage
%  score sampling methods for both the fixed and nonfixed-rank approximants

%% load the dataset

%%% TODO: when calculating the mean errors and timings, need to account for
% fact that each method may have a different number of sample points in new
% datasets
load(datasetfname);

%% setup plotting

if printflag
    make_it_tight = true;
    subplot = @(m,n,p) subtightplot (m, n, p, [0.082 0.07], [0.04 0.04], [0.1 0.06]);
    if ~make_it_tight
        clear subplot;  
    end 

    fontsize = 15;
    width = 6.2;
    height = 12.4;
    basename = fullfile(outdir, savedata.in.datasetbasename);
    printfig = @(figname) printcf([basename figname '.pdf'], fontsize, width, height);
    printpanel = @(panelname, figname) printcf([basename figname '-' panelname '.pdf'], fontsize, width, width);

    parameterswidth = 6.2;
    parametersheight = 6.2;
    printparameters = @(figname) printcf([basename figname '.pdf'], fontsize, parameterswidth, parametersheight);
    
    timingwidth = 6.2;
    timingheight = 6.2;
    printtiming = @(figname) printcf([basename figname '.pdf'], fontsize, timingwidth, timingheight);
end

%% calculate the mean errors and timings

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
    
    % approxlev
    means = mean(savedata.approxlevData(lidx).specerr,2);
    nonfixed_approxlev_specerr(lidx) = means(1);
    fixed_approxlev_specerr(lidx) = means(2);
    
    means = mean(savedata.approxlevData(lidx).froerr,2);
    nonfixed_approxlev_froerr(lidx) = means(1);
    fixed_approxlev_froerr(lidx) = means(2);
    
    means = mean(savedata.approxlevData(lidx).trerr,2);
    nonfixed_approxlev_trerr(lidx) = means(1);
    fixed_approxlev_trerr(lidx) = means(2);
    
    means = mean(savedata.approxlevData(lidx).timings, 2);
    nonfixed_approxlev_timing(lidx) = means(1);
    fixed_approxlev_timing(lidx) = means(2);
    
    % froblev
    means = mean(savedata.froblevData(lidx).specerr,2);
    nonfixed_froblev_specerr(lidx) = means(1);
    fixed_froblev_specerr(lidx) = means(2);
    
    means = mean(savedata.froblevData(lidx).froerr,2);
    nonfixed_froblev_froerr(lidx) = means(1);
    fixed_froblev_froerr(lidx) = means(2);
    
    means = mean(savedata.froblevData(lidx).trerr,2);
    nonfixed_froblev_trerr(lidx) = means(1);
    fixed_froblev_trerr(lidx) = means(2);
    
    means = mean(savedata.froblevData(lidx).timings, 2);
    nonfixed_froblev_timing(lidx) = means(1);
    fixed_froblev_timing(lidx) = means(2);
    
    % speclev
    means = mean(savedata.speclevData(lidx).specerr,2);
    nonfixed_speclev_specerr(lidx) = means(1);
    fixed_speclev_specerr(lidx) = means(2);
    
    means = mean(savedata.speclevData(lidx).froerr,2);
    nonfixed_speclev_froerr(lidx) = means(1);
    fixed_speclev_froerr(lidx) = means(2);
    
    means = mean(savedata.speclevData(lidx).trerr,2);
    nonfixed_speclev_trerr(lidx) = means(1);
    fixed_speclev_trerr(lidx) = means(2);
    
    means = mean(savedata.speclevData(lidx).timings, 2);
    nonfixed_speclev_timing(lidx) = means(1);
    fixed_speclev_timing(lidx) = means(2);
    
end

%% Plot the time-error tradeoffs of the fixed and nonfixed rank exact methods
lw = 1.5;
ms = 10;
if colorflag
    simple_style = 'ks-';
    simple_color = 'k';
    srft_style = 'bo-';
    srft_color = 'b';
    gaussian_style = 'g*-';
    gaussian_color = 'g';
    levscore_style = 'rd-';
    levscore_color = 'r';
else
    simple_style = 's-';
    simple_color = [.1 .1 .1];
    srft_style = 'o-';
    srft_color = [.2 .2 .2];
    gaussian_style = 'v-';
    gaussian_color = [.3 .3 .3];
    levscore_style = 'd-';
    levscore_color = [.4 .4 .4];
end
legendloc = 'Northeast';

exactbase_struct = struct(...
    'seriesnames', {{'unif', 'srft', 'gaussian', 'levscore'}}, ...
    'xlabel', 'time (s)', ...
    'plottype', 'NOTsemilogx', ...
    'clipaxis', true, ...
    'legendloc', legendloc, ...
    'styles', {{simple_style, srft_style, ...
    gaussian_style, levscore_style}}, ...
    'lw', {{lw, lw, lw, lw}}, ...
    'ms', {{ms, ms, ms, ms}}, ...
    'colors', {{simple_color, srft_color, ...
    gaussian_color, levscore_color}}, ...
    'mcolors', {{simple_color, srft_color, ...
    gaussian_color, levscore_color}});

exactbase_struct.x = { ...
    nonfixed_simple_timing, ...
    nonfixed_srft_timing, ...
    nonfixed_gaussian_timing, ...
    nonfixed_levscore_timing};

%% Exact method, nonfixed rank plots

exactmethods_nonfixedrank_spec_struct = exactbase_struct;
exactmethods_nonfixedrank_frob_struct = exactbase_struct;
exactmethods_nonfixedrank_trace_struct = exactbase_struct;

exactmethods_nonfixedrank_spec_struct.series = { ...
    nonfixed_simple_specerr/savedata.optspecerr, ...
    nonfixed_srft_specerr/savedata.optspecerr, ...
    nonfixed_gaussian_specerr/savedata.optspecerr, ...
    nonfixed_levscore_specerr/savedata.optspecerr};
exactmethods_nonfixedrank_spec_struct.plotname =  ...
    ['$\|\mathbf{A} - \mathbf{C} \mathbf{W}^\dagger \mathbf{C}^T\|_2'...
    '/\|\mathbf{A} - \mathbf{A}_k\|_2$'];

exactmethods_nonfixedrank_frob_struct.series = { ...
    nonfixed_simple_froerr/savedata.optfroerr, ...
    nonfixed_srft_froerr/savedata.optfroerr, ...
    nonfixed_gaussian_froerr/savedata.optfroerr, ...
    nonfixed_levscore_froerr/savedata.optfroerr};
exactmethods_nonfixedrank_frob_struct.plotname =  ...
    ['$\|\mathbf{A} - \mathbf{C} \mathbf{W}^\dagger \mathbf{C}^T\|_F'...
    '/\|\mathbf{A} - \mathbf{A}_k\|_F$'];

exactmethods_nonfixedrank_trace_struct.series = { ...
    nonfixed_simple_trerr/savedata.opttrerr, ...
    nonfixed_srft_trerr/savedata.opttrerr, ...
    nonfixed_gaussian_trerr/savedata.opttrerr, ...
    nonfixed_levscore_trerr/savedata.opttrerr};
exactmethods_nonfixedrank_trace_struct.plotname =  ...
    ['$\|\mathbf{A} - \mathbf{C} \mathbf{W}^\dagger \mathbf{C}^T'...
    '\|_\star/\|\mathbf{A} - \mathbf{A}_k\|_\star$'];

figure();
subplot(3,1,1);
plotpanel_multiplex(exactmethods_nonfixedrank_spec_struct, false);
subplot(3,1,2);
plotpanel_multiplex(exactmethods_nonfixedrank_frob_struct, false);
subplot(3,1,3);
plotpanel_multiplex(exactmethods_nonfixedrank_trace_struct, true);

% if storing plots, store the plot with all three norms at once, then
% separate plots for the individual norms
if printflag
    printfig('exact-methods-nonfixed-rank-time-error-curves');

	figure();
	plotpanel_multiplex(exactmethods_nonfixedrank_spec_struct);
	printpanel('specerr', 'exact-methods-nonfixed-rank-time-error-curves');

	figure();
	plotpanel_multiplex(exactmethods_nonfixedrank_frob_struct);
	printpanel('froberr', 'exact-methods-nonfixed-rank-time-error-curves');

	figure();
	plotpanel_multiplex(exactmethods_nonfixedrank_trace_struct);
	printpanel('trerr', 'exact-methods-nonfixed-rank-time-error-curves');
end

%% Exact methods, fixed rank plots 

exactmethods_fixedrank_spec_struct = exactbase_struct;
exactmethods_fixedrank_frob_struct = exactbase_struct;
exactmethods_fixedrank_trace_struct = exactbase_struct;

exactmethods_fixedrank_spec_struct.series = { ...
    fixed_simple_specerr/savedata.optspecerr, ...
    fixed_srft_specerr/savedata.optspecerr, ...
    fixed_gaussian_specerr/savedata.optspecerr, ...
    fixed_levscore_specerr/savedata.optspecerr};
exactmethods_fixedrank_spec_struct.plotname = ...
    ['$\|\mathbf{A} - \mathbf{C} \mathbf{W}_k^\dagger \mathbf{C}^T\|_2'...
    '/\|\mathbf{A} - \mathbf{A}_k\|_2$'];

exactmethods_fixedrank_frob_struct.series = { ...
    fixed_simple_froerr/savedata.optfroerr, ...
    fixed_srft_froerr/savedata.optfroerr, ...
    fixed_gaussian_froerr/savedata.optfroerr, ...
    fixed_levscore_froerr/savedata.optfroerr};
exactmethods_fixedrank_frob_struct.plotname = ...
    ['$\|\mathbf{A} - \mathbf{C} \mathbf{W}_k^\dagger \mathbf{C}^T\|_F'...
    '/\|\mathbf{A} - \mathbf{A}_k\|_F$'];

exactmethods_fixedrank_trace_struct.series = { ...
    fixed_simple_trerr/savedata.opttrerr, ...
    fixed_srft_trerr/savedata.opttrerr, ...
    fixed_gaussian_trerr/savedata.opttrerr, ...
    fixed_levscore_trerr/savedata.opttrerr};
exactmethods_fixedrank_trace_struct.plotname = ...
    ['$\|\mathbf{A} - \mathbf{C} \mathbf{W}_k^\dagger \mathbf{C}^T'...
    '\|_\star/\|\mathbf{A} - \mathbf{A}_k\|_\star$'];

figure();
subplot(3,1,1);
plotpanel_multiplex(exactmethods_fixedrank_spec_struct, false);
subplot(3,1,2);
plotpanel_multiplex(exactmethods_fixedrank_frob_struct, false);
subplot(3,1,3);
plotpanel_multiplex(exactmethods_fixedrank_trace_struct, true);

% if storing plots, store the plot with all three norms at once, then
% separate plots for the individual norms
if printflag
    printfig('exact-methods-fixed-rank-time-error-curves');

	figure();
	plotpanel_multiplex(exactmethods_fixedrank_spec_struct);
	printpanel('specerr', 'exact-methods-fixed-rank-time-error-curves');

	figure();
	plotpanel_multiplex(exactmethods_fixedrank_frob_struct);
	printpanel('froberr', 'exact-methods-fixed-rank-time-error-curves');

	figure();
	plotpanel_multiplex(exactmethods_fixedrank_trace_struct);
	printpanel('trerr', 'exact-methods-fixed-rank-time-error-curves');
end

%% Plot the errors of the fixed and nonfixed inexact leverage score methods
% display the uniform and leverage sampling errors for calibration
lw = 1.5;
ms = 10;
if colorflag
    simple_style = 'ks-';
    simple_color = 'k';
    froblev_style = 'bo-';
    froblev_color = 'b';
    levscore_style = 'rd-';
    levscore_color = 'r';
    approxlev_style = 'g*-';
    approxlev_color = 'g';
    speclev_style = 'm^-';
    speclev_color = 'm';
else
    simple_style = 's-';
    simple_color = [.1 .1 .1];
    froblev_style = 'o-';
    froblev_color = [.2 .2 .2];
    approxlev_style = 'v-';
    approxlev_color = [.3 .3 .3];
    levscore_style = 'd-';
    levscore_color = [.4 .4 .4];
    speclev_style = '^-';
    speclev_color = [.5 .5 .5];
end
legendloc = 'Northeast';

inexactbase_struct = struct(...
    'seriesnames', {{'levscore', 'unif', 'power', 'frob levscore', ...
                     'spec levscore'}}, ...
    'xlabel', 'timing (seconds)', ...
    'plottype', 'NOTsemilogx', ...
    'clipaxis', true, ...
    'legendloc', legendloc, ...
    'styles', {{levscore_style, simple_style, approxlev_style, ...
                froblev_style, speclev_style}}, ...
    'lw', {{lw, lw, lw, lw, lw}}, ...
    'ms', {{ms, ms, ms, ms, ms}}, ...
    'colors', {{levscore_color, simple_color, approxlev_color, ...
                froblev_color, speclev_color}}, ...
    'mcolors', {{levscore_color, simple_color, approxlev_color, ...
                 froblev_color, speclev_color}});

inexactbase_struct.x = { ...
    nonfixed_levscore_timing, ...
    nonfixed_simple_timing, ...
    nonfixed_approxlev_timing, ...
    nonfixed_froblev_timing, ...
    nonfixed_speclev_timing};

%% Inexact methods, nonfixed rank plots

inexactmethods_nonfixedrank_spec_struct = inexactbase_struct;
inexactmethods_nonfixedrank_frob_struct = inexactbase_struct;
inexactmethods_nonfixedrank_trace_struct = inexactbase_struct;

inexactmethods_nonfixedrank_spec_struct.series = { ...
    nonfixed_levscore_specerr/savedata.optspecerr, ...
    nonfixed_simple_specerr/savedata.optspecerr, ...
    nonfixed_approxlev_specerr/savedata.optspecerr, ...
    nonfixed_froblev_specerr/savedata.optspecerr, ...
    nonfixed_speclev_specerr/savedata.optspecerr};
inexactmethods_nonfixedrank_spec_struct.plotname =  ...
    ['$\|\mathbf{A} - \mathbf{C} \mathbf{W}^\dagger \mathbf{C}^T\|_2'...
    '/\|\mathbf{A} - \mathbf{A}_k\|_2$'];

inexactmethods_nonfixedrank_frob_struct.series = { ...
    nonfixed_levscore_froerr/savedata.optfroerr, ...
    nonfixed_simple_froerr/savedata.optfroerr, ...
    nonfixed_approxlev_froerr/savedata.optfroerr, ...
    nonfixed_froblev_froerr/savedata.optfroerr, ...
    nonfixed_speclev_froerr/savedata.optfroerr};    
inexactmethods_nonfixedrank_frob_struct.plotname = ...
    ['$\|\mathbf{A} - \mathbf{C} \mathbf{W}^\dagger \mathbf{C}^T\|_F/'...
    '\|\mathbf{A} - \mathbf{A}_k\|_F$'];

inexactmethods_nonfixedrank_trace_struct.series = { ...
    nonfixed_levscore_trerr/savedata.opttrerr, ...
    nonfixed_simple_trerr/savedata.opttrerr, ...
    nonfixed_approxlev_trerr/savedata.opttrerr, ...
    nonfixed_froblev_trerr/savedata.opttrerr, ...
    nonfixed_speclev_trerr/savedata.opttrerr};    
inexactmethods_nonfixedrank_trace_struct.plotname =  ...
    ['$\|\mathbf{A} - \mathbf{C} \mathbf{W}^\dagger \mathbf{C}^T'...
    '\|_\star/\|\mathbf{A} - \mathbf{A}_k\|_\star$'];

figure();
subplot(3,1,1);
plotpanel_multiplex(inexactmethods_nonfixedrank_spec_struct, false);
subplot(3,1,2);
plotpanel_multiplex(inexactmethods_nonfixedrank_frob_struct, false);
subplot(3,1,3);
plotpanel_multiplex(inexactmethods_nonfixedrank_trace_struct, true);

% if storing plots, store the plot with all three norms at once, then
% separate plots for the individual norms
if printflag
    printfig('inexact-methods-nonfixed-rank-time-error-curves');
    
    figure();
    plotpanel_multiplex(inexactmethods_nonfixedrank_spec_struct);
	printpanel('specerr', 'inexact-methods-nonfixed-rank-time-error-curves');

	figure();
	plotpanel_multiplex(inexactmethods_nonfixedrank_frob_struct);
	printpanel('froberr', 'inexact-methods-nonfixed-rank-time-error-curves');

	figure();
	plotpanel_multiplex(inexactmethods_nonfixedrank_trace_struct);
	printpanel('trerr', 'inexact-methods-nonfixed-rank-time-error-curves');
end

%% Inexact methods, fixed rank plots

inexactmethods_fixedrank_spec_struct = inexactbase_struct;
inexactmethods_fixedrank_frob_struct = inexactbase_struct;
inexactmethods_fixedrank_trace_struct = inexactbase_struct;

inexactmethods_fixedrank_spec_struct.series = { ...
    fixed_levscore_specerr/savedata.optspecerr, ...
    fixed_simple_specerr/savedata.optspecerr, ...
    fixed_approxlev_specerr/savedata.optspecerr, ...
    fixed_froblev_specerr/savedata.optspecerr, ...
    fixed_speclev_specerr/savedata.optspecerr};
inexactmethods_fixedrank_spec_struct.plotname =  ...
    ['$\|\mathbf{A} - \mathbf{C} \mathbf{W}_k^\dagger \mathbf{C}^T\|_2'...
    '/\|\mathbf{A} - \mathbf{A}_k\|_2$'];

inexactmethods_fixedrank_frob_struct.series = { ...
    fixed_levscore_froerr/savedata.optfroerr, ...
    fixed_simple_froerr/savedata.optfroerr, ...
    fixed_approxlev_froerr/savedata.optfroerr, ...
    fixed_froblev_froerr/savedata.optfroerr, ...
    fixed_speclev_froerr/savedata.optfroerr};    
inexactmethods_fixedrank_frob_struct.plotname = ...
    ['$\|\mathbf{A} - \mathbf{C} \mathbf{W}_k^\dagger \mathbf{C}^T\|_F'...
    '/\|\mathbf{A} - \mathbf{A}_k\|_F$'];

inexactmethods_fixedrank_trace_struct.series = { ...
    fixed_levscore_trerr/savedata.opttrerr, ...
    fixed_simple_trerr/savedata.opttrerr, ...
    fixed_approxlev_trerr/savedata.opttrerr, ...
    fixed_froblev_trerr/savedata.opttrerr, ...
    fixed_speclev_trerr/savedata.opttrerr};    
inexactmethods_fixedrank_trace_struct.plotname = ...
    ['$\|\mathbf{A} - \mathbf{C} \mathbf{W}_k^\dagger \mathbf{C}^T'...
    '\|_\star/\|\mathbf{A} - \mathbf{A}_k\|_\star$'];

figure();
subplot(3,1,1);
plotpanel_multiplex(inexactmethods_fixedrank_spec_struct, false);
subplot(3,1,2);
plotpanel_multiplex(inexactmethods_fixedrank_frob_struct, false);
subplot(3,1,3);
plotpanel_multiplex(inexactmethods_fixedrank_trace_struct, true);

% if storing plots, store the plot with all three norms at once, then
% separate plots for the individual norms
if printflag
    printfig('inexact-methods-fixed-rank-time-error-curves');
    
    figure();
    plotpanel_multiplex(inexactmethods_fixedrank_spec_struct);
	printpanel('specerr', 'inexact-methods-fixed-rank-time-error-curves');

	figure();
	plotpanel_multiplex(inexactmethods_fixedrank_frob_struct);
	printpanel('froberr', 'inexact-methods-fixed-rank-time-error-curves');

	figure();
	plotpanel_multiplex(inexactmethods_fixedrank_trace_struct);
	printpanel('trerr', 'inexact-methods-fixed-rank-time-error-curves');
end

%% close all figures
if printflag
  close all
end

end
