function visualize_tallthin(datasetfname1, datasetfname2, printflag, outdir)
% visualize(datasetfname1, datasetfname2, printflag, outdir)
%
% Visualizes the information on Nystrom approximation stored in the files
% datasetfname1 and datasetfname2. If printflag is true, saves the generated graphs to pdfs
% with appropriate names. If outdir is present and printflag is true, the 
% pdfs are saved to this directory.
%
% datasetfname1 is assumed to contain the information on all inexact
% leverage score methods other than the one using the tall thin
% approximation of leverage scores
%
% datasetfname2 is assumed to contain only the information on the tall
% thin approximation algorithm
%
% Graphs generated:
%  - spectral, frobenius, and trace norm errors of the inexact leverage
%  score sampling methods for both the fixed and nonfixed-rank approximants
%  - timings for the inexact sampling methods (only the nonfixed-rank case)

%% setup plotting
load(datasetfname1) % so basename can be defined

if printflag
    make_it_tight = true;
    subplot = @(m,n,p) subtightplot (m, n, p, [0.082 0.07], [0.04 0.04], [0.1 0.04]);
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
    
    % leverage
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
end

load(datasetfname2);

for lidx = 1:length(savedata.in.lvals)
    % tall thin lev score approx alg
    means = mean(savedata.tallthinData(lidx).specerr,2);
    nonfixed_tallthin_specerr(lidx) = means(1);
    fixed_tallthin_specerr(lidx) = means(2);
    
    means = mean(savedata.tallthinData(lidx).froerr,2);
    nonfixed_tallthin_froerr(lidx) = means(1);
    fixed_tallthin_froerr(lidx) = means(2);
    
    means = mean(savedata.tallthinData(lidx).trerr,2);
    nonfixed_tallthin_trerr(lidx) = means(1);
    fixed_tallthin_trerr(lidx) = means(2);
    
    means = mean(savedata.tallthinData(lidx).timings, 2);
    nonfixed_tallthin_timing(lidx) = means(1);
    fixed_tallthin_timing(lidx) = means(2);
    
    % QR-computed levscores
    means = mean(savedata.qrlevscoreData(lidx).specerr, 2);
    nonfixed_qrlevscore_specerr(lidx) = means(1);
    fixed_qrlevscore_specerr(lidx) = means(2);
    
    means = mean(savedata.qrlevscoreData(lidx).froerr, 2);
    nonfixed_qrlevscore_froerr(lidx) = means(1);
    fixed_qrlevscore_froerr(lidx) = means(2);
    
    means = mean(savedata.qrlevscoreData(lidx).trerr, 2);
    nonfixed_qrlevscore_trerr(lidx) = means(1);
    fixed_qrlevscore_trerr(lidx) = means(2);
    
    means = mean(savedata.qrlevscoreData(lidx).timings, 2);
    nonfixed_qrlevscore_timing(lidx) = means(1);
    fixed_qrlevscore_timing(lidx) = means(2);
    
end


%% Plot the errors of the inexact leverage score methods
% display the uniform and QR leverage sampling errors for calibration
lw = 1.5;
ms = 10;
simple_style = 's-';
simple_color = [.1 .1 .1];
froblev_style = 'o-';
froblev_color = [.2 .2 .2];
approxlev_style = 'v-';
approxlev_color = [.3 .3 .3];
qrlevscore_style = 'd-';
qrlevscore_color = [.4 .4 .4];
speclev_style = '^-';
speclev_color = [.5 .5 .5];
tallthin_style = 'x-';
tallthin_color = [.6 .6 .6];
srft_style = '+-'; 
srft_color = [.7 .7 .7];

legendloc = 'Northeast';

inexactbase_struct = struct(...
    'seriesnames', {{'QR lev', 'unif', 'power', 'frob lev', ...
       'spec lev', 'srft', 'Alg 1'}}, ...
    'xlabel', '$\ell$', ...
    'clipaxis', true, ...
    'legendloc', legendloc, ...
    'x', savedata.in.lvals, ...
    'styles', {{qrlevscore_style, simple_style, approxlev_style, ...
        froblev_style, speclev_style, srft_style, tallthin_style}}, ...
    'lw', {{lw, lw, lw, lw, lw, lw, lw}}, ...
    'ms', {{ms, ms, ms, ms, ms, ms, ms}}, ...
    'colors', {{qrlevscore_color, simple_color, approxlev_color, ...
        froblev_color, speclev_color, srft_color, tallthin_color}}, ...
    'mcolors', {{qrlevscore_color, simple_color, approxlev_color, ...
        froblev_color, speclev_color, srft_color, tallthin_color}});

%% Inexact method, nonfixed rank plots

inexactmethods_nonfixedrank_spec_struct = inexactbase_struct;
inexactmethods_nonfixedrank_frob_struct = inexactbase_struct;
inexactmethods_nonfixedrank_trace_struct = inexactbase_struct;

inexactmethods_nonfixedrank_spec_struct.series = { ...
    nonfixed_qrlevscore_specerr/savedata.optspecerr, ...
    nonfixed_simple_specerr/savedata.optspecerr, ...
    nonfixed_approxlev_specerr/savedata.optspecerr, ...
    nonfixed_froblev_specerr/savedata.optspecerr, ...
    nonfixed_speclev_specerr/savedata.optspecerr, ...
    nonfixed_srft_specerr/savedata.optspecerr, ...
    nonfixed_tallthin_specerr/savedata.optspecerr};
inexactmethods_nonfixedrank_spec_struct.plotname = ...
    ['$\|\mathbf{A} - \mathbf{C} \mathbf{W}^\dagger \mathbf{C}^T\|_2'...
    '/\|\mathbf{A} - \mathbf{A}_k\|_2$'];

inexactmethods_nonfixedrank_frob_struct.series = { ...
    nonfixed_qrlevscore_froerr/savedata.optfroerr, ...
    nonfixed_simple_froerr/savedata.optfroerr, ...
    nonfixed_approxlev_froerr/savedata.optfroerr, ...
    nonfixed_froblev_froerr/savedata.optfroerr, ...
    nonfixed_speclev_froerr/savedata.optfroerr, ...
    nonfixed_srft_froerr/savedata.optfroerr, ...
    nonfixed_tallthin_froerr/savedata.optfroerr};
inexactmethods_nonfixedrank_frob_struct.plotname = ...
    ['$\|\mathbf{A} - \mathbf{C} \mathbf{W}^\dagger \mathbf{C}^T\|_F'...
    '/\|\mathbf{A} - \mathbf{A}_k\|_F$'];

inexactmethods_nonfixedrank_trace_struct.series = { ...
    nonfixed_qrlevscore_trerr/savedata.opttrerr, ...
    nonfixed_simple_trerr/savedata.opttrerr, ...
    nonfixed_approxlev_trerr/savedata.opttrerr, ...
    nonfixed_froblev_trerr/savedata.opttrerr, ...
    nonfixed_speclev_trerr/savedata.opttrerr, ...
    nonfixed_srft_trerr/savedata.opttrerr, ...
    nonfixed_tallthin_trerr/savedata.opttrerr};   
inexactmethods_nonfixedrank_trace_struct.plotname = ...
    ['$\|\mathbf{A} - \mathbf{C} \mathbf{W}^\dagger \mathbf{C}^T'...
    '\|_\star/\|\mathbf{A} - \mathbf{A}_k\|_\star$'];

figure();
subplot(3,1,1);
plotpanel(inexactmethods_nonfixedrank_spec_struct);
subplot(3,1,2);
plotpanel(inexactmethods_nonfixedrank_frob_struct);
subplot(3,1,3);
plotpanel(inexactmethods_nonfixedrank_trace_struct);

% if storing plots, store the plot with all three norms at once, then
% separate plots for the individual norms
if printflag
    printfig('inexact-methods-plus-tallthin-nonfixed-rank-errors');
    
    figure();
    plotpanel(inexactmethods_nonfixedrank_spec_struct);
	printpanel('specerr', 'inexact-methods-plus-tallthin-nonfixed-rank-errors');

	figure();
	plotpanel(inexactmethods_nonfixedrank_frob_struct);
	printpanel('froberr', 'inexact-methods-plus-tallthin-nonfixed-rank-errors');

	figure();
	plotpanel(inexactmethods_nonfixedrank_trace_struct);
	printpanel('trerr', 'inexact-methods-plus-tallthin-nonfixed-rank-errors');
end

%% Inexact methods, fixed rank plots

inexactmethods_fixedrank_spec_struct = inexactbase_struct;
inexactmethods_fixedrank_frob_struct = inexactbase_struct;
inexactmethods_fixedrank_trace_struct = inexactbase_struct;

inexactmethods_fixedrank_spec_struct.series = { ...
    fixed_qrlevscore_specerr/savedata.optspecerr, ...
    fixed_simple_specerr/savedata.optspecerr, ...
    fixed_approxlev_specerr/savedata.optspecerr, ...
    fixed_froblev_specerr/savedata.optspecerr, ...
    fixed_speclev_specerr/savedata.optspecerr, ...
    fixed_srft_specerr/savedata.optspecerr, ...
    fixed_tallthin_specerr/savedata.optspecerr};
inexactmethods_fixedrank_spec_struct.plotname = ...
    ['$\|\mathbf{A} - \mathbf{C} \mathbf{W}_k^\dagger \mathbf{C}^T\|_2'...
    '/\|\mathbf{A} - \mathbf{A}_k\|_2$'];

inexactmethods_fixedrank_frob_struct.series = { ...
    fixed_qrlevscore_froerr/savedata.optfroerr, ...
    fixed_simple_froerr/savedata.optfroerr, ...
    fixed_approxlev_froerr/savedata.optfroerr, ...
    fixed_froblev_froerr/savedata.optfroerr, ...
    fixed_speclev_froerr/savedata.optfroerr, ...
    fixed_srft_froerr/savedata.optfroerr, ...
    fixed_tallthin_froerr/savedata.optfroerr};
inexactmethods_fixedrank_frob_struct.plotname = ...
    ['$\|\mathbf{A} - \mathbf{C} \mathbf{W}_k^\dagger \mathbf{C}^T\|_F'...
    '/\|\mathbf{A} - \mathbf{A}_k\|_F$'];


inexactmethods_fixedrank_trace_struct.series = { ...
    fixed_qrlevscore_trerr/savedata.opttrerr, ...
    fixed_simple_trerr/savedata.opttrerr, ...
    fixed_approxlev_trerr/savedata.opttrerr, ...
    fixed_froblev_trerr/savedata.opttrerr, ...
    fixed_speclev_trerr/savedata.opttrerr, ...
    fixed_srft_trerr/savedata.opttrerr, ...
    fixed_tallthin_trerr/savedata.opttrerr};   
inexactmethods_fixedrank_trace_struct.plotname =...
    ['$\|\mathbf{A} - \mathbf{C} \mathbf{W}_k^\dagger \mathbf{C}^T'...
    '\|_\star/\|\mathbf{A} - \mathbf{A}_k\|_\star$'];
figure();
subplot(3,1,1);
plotpanel(inexactmethods_fixedrank_spec_struct);
subplot(3,1,2);
plotpanel(inexactmethods_fixedrank_frob_struct);
subplot(3,1,3);
plotpanel(inexactmethods_fixedrank_trace_struct);

% if storing plots, store the plot with all three norms at once, then
% separate plots for the individual norms
if printflag
    printfig('inexact-methods-plus-tallthin-fixed-rank-errors');
    
    figure();
    plotpanel(inexactmethods_fixedrank_spec_struct);
	printpanel('specerr', 'inexact-methods-plus-tallthin-fixed-rank-errors');

	figure();
	plotpanel(inexactmethods_fixedrank_frob_struct);
	printpanel('froberr', 'inexact-methods-plus-tallthin-fixed-rank-errors');

	figure();
	plotpanel(inexactmethods_fixedrank_trace_struct);
	printpanel('trerr', 'inexact-methods-plus-tallthin-fixed-rank-errors');
end

%% Show the timings for the inexact methods
inexacttiming_struct = inexactbase_struct;
inexacttiming_struct.legendloc = legendloc;
inexacttiming_struct.ylabel = 'time (s)';
inexacttiming_struct.plottype = 'semilogy';
rmfield(inexacttiming_struct, 'clipaxis');

inexacttiming_struct.series = { ...
    nonfixed_qrlevscore_timing, ...
    nonfixed_simple_timing, ...
    nonfixed_approxlev_timing, ...
    nonfixed_froblev_timing, ...
    nonfixed_speclev_timing, ...
    nonfixed_srft_timing, ...
    nonfixed_tallthin_timing};

figure();
plotpanel(inexacttiming_struct);

if printflag
    printtiming('inexact-methods-plus-tallthin-timings');
end

%% close all figures
if printflag
   close all
end

end