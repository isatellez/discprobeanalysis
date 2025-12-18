function dump_repo_to_text(rootDir, outFile)
% Dump all .m files under rootDir into a single text file.

    if nargin < 1 || isempty(rootDir)
        rootDir = pwd;
    end

    if nargin < 2 || isempty(outFile)
        outFile = fullfile(rootDir, 'repo_mfiles_dump.txt');
    end

    files = dir(fullfile(rootDir, '**', '*.m'));

    if isempty(files)
        warning('No .m files found under %s', rootDir);
        return;
    end

    fidOut = fopen(outFile, 'w');
    if fidOut == -1
        error('Could not open output file: %s', outFile);
    end

    cleaner = onCleanup(@() fclose(fidOut));

    for k = 1:numel(files)
        fpath = fullfile(files(k).folder, files(k).name);

        relpath = strrep(fpath, [rootDir filesep], '');

        fprintf(fidOut, '%% ===== FILE: %s =====\n', relpath);
        fprintf(fidOut, '%% ===== START =====\n');

        fidIn = fopen(fpath, 'r');
        if fidIn == -1
            fprintf(fidOut, '%% [Could not open file]\n\n');
            continue;
        end

        tline = fgetl(fidIn);
        while ischar(tline)
            fprintf(fidOut, '%s\n', tline);
            tline = fgetl(fidIn);
        end
        fclose(fidIn);

        fprintf(fidOut, '%% ===== END FILE: %s =====\n\n', relpath);
    end

    fprintf('Wrote %d files into %s\n', numel(files), outFile);
end
