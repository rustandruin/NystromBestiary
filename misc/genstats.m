function [nystrom, srft, gaussian, levscore] = genstats(saveddata)

% [nystrom, srft, gaussian, levscore] = genstats(dataset)
%
% Generates the min/mean/max stats for the errors of the given dataset
%

computestats = @(array) [min(array), mean(array), max(array)];

optspecerr = saveddata.optspecerr;
optfroerr = saveddata.optfroerr;
opttrerr = saveddata.opttrerr;

numvaluesofsamples = size(saveddata.simpleData, 2);

for lidx = 1:numvaluesofsamples
    nystrom.nonrestricted_spec_stats(lidx, :) = ...
        computestats(saveddata.simpleData(lidx).specerr(1, :)/optspecerr);
    nystrom.nonrestricted_fro_stats(lidx, :) = ...
        computestats(saveddata.simpleData(lidx).froerr(1, :)/optfroerr); 
    nystrom.nonrestricted_tr_stats(lidx, :) = ...
        computestats(saveddata.simpleData(lidx).trerr(1, :)/opttrerr); 
        
    nystrom.restricted_spec_stats(lidx, :) = ...
        computestats(saveddata.simpleData(lidx).specerr(2, :)/optspecerr);
    nystrom.restricted_fro_stats(lidx, :) = ...
        computestats(saveddata.simpleData(lidx).froerr(2, :)/optfroerr); 
    nystrom.restricted_tr_stats(lidx, :) = ...
        computestats(saveddata.simpleData(lidx).trerr(2, :)/opttrerr);
    
    srft.nonrestricted_spec_stats(lidx, :) = ...
        computestats(saveddata.srftData(lidx).specerr(1, :)/optspecerr);
    srft.nonrestricted_fro_stats(lidx, :) = ...
        computestats(saveddata.srftData(lidx).froerr(1, :)/optfroerr); 
    srft.nonrestricted_tr_stats(lidx, :) = ...
        computestats(saveddata.srftData(lidx).trerr(1, :)/opttrerr); 
        
    srft.restricted_spec_stats(lidx, :) = ...
        computestats(saveddata.srftData(lidx).specerr(2, :)/optspecerr);
    srft.restricted_fro_stats(lidx, :) = ...
        computestats(saveddata.srftData(lidx).froerr(2, :)/optfroerr); 
    srft.restricted_tr_stats(lidx, :) = ...
        computestats(saveddata.srftData(lidx).trerr(2, :)/opttrerr); 
    
    gaussian.nonrestricted_spec_stats(lidx, :) = ...
        computestats(saveddata.gaussianData(lidx).specerr(1, :)/optspecerr);
    gaussian.nonrestricted_fro_stats(lidx, :) = ...
        computestats(saveddata.gaussianData(lidx).froerr(1, :)/optfroerr); 
    gaussian.nonrestricted_tr_stats(lidx, :) = ...
        computestats(saveddata.gaussianData(lidx).trerr(1, :)/opttrerr); 
        
    gaussian.restricted_spec_stats(lidx, :) = ...
        computestats(saveddata.gaussianData(lidx).specerr(2, :)/optspecerr);
    gaussian.restricted_fro_stats(lidx, :) = ...
        computestats(saveddata.gaussianData(lidx).froerr(2, :)/optfroerr); 
    gaussian.restricted_tr_stats(lidx, :) = ...
        computestats(saveddata.gaussianData(lidx).trerr(2, :)/opttrerr); 
    
    levscore.nonrestricted_spec_stats(lidx, :) = ...
        computestats(saveddata.levscoreData(lidx).specerr(1, :)/optspecerr);
    levscore.nonrestricted_fro_stats(lidx, :) = ...
        computestats(saveddata.levscoreData(lidx).froerr(1, :)/optfroerr); 
    levscore.nonrestricted_tr_stats(lidx, :) = ...
        computestats(saveddata.levscoreData(lidx).trerr(1, :)/opttrerr); 
                
    levscore.restricted_spec_stats(lidx, :) = ...
        computestats(saveddata.levscoreData(lidx).specerr(2, :)/optspecerr);
    levscore.restricted_fro_stats(lidx, :) = ...
        computestats(saveddata.levscoreData(lidx).froerr(2, :)/optfroerr); 
    levscore.restricted_tr_stats(lidx, :) = ...
        computestats(saveddata.levscoreData(lidx).trerr(2, :)/opttrerr); 
end

end