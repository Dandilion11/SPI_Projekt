function [modelname, paths, opts, prev_ident_pars] = temp_opts_m1_o2()
    modelname = 'SF_Modell1_O2.mat';
    scriptDir = '/home/peter/Dokumente/Master/Semester2/StrukturParameteridentifikation/Projekt/SPI_Projekt/Modelle';
    paths.models    = '/home/peter/Dokumente/Master/Semester2/StrukturParameteridentifikation/Projekt/SPI_Projekt/Modelle/../Modelle';
    paths.results   = fullfile(scriptDir, 'results');
    paths.functions = fullfile(scriptDir, 'functions');
    opts.algorithm   = 1;
    opts.maxLietime  = 120;
    opts.decomp      = 0;
    opts.checkObser  = 1;
    opts.forcedecomp = 0;
    opts.tol         = 1e-6;
    opts.multiexp    = 0;
    opts.nnzDerU     = inf;
    opts.nnzDerW     = 0;
    prev_ident_pars  = [];
end
