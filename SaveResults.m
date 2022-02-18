function [] = SaveResults(C,M,m,n,resultsDir, name)   
    %%%GUARDAR LAS CONCENTRACIONES
        H = reshape(C(1,:,:), m, n);
        E = reshape(C(2,:,:), m, n);
        stains(:,:,1) = H;
        stains(:,:,2) = E;
        stains = reshape(stains, m, n, 2);
        
    
    save(fullfile(resultsDir,strcat(name,'_M','.mat')),'M');
    outputFileName = fullfile('./', [name '.' 'Results.mat'])
        if ~exist(fileparts(outputFileName), 'dir')
            mkdir(fileparts(outputFileName))
        end
        %save(outputFileName, 'stains', 'M', 'whiteBalanceLevel', '-v7.3');

        % Save stains, M and whiteBalanceLevel
        chunk_size = min([1024 1024 1], size(stains)); % define chunk size
        if ~exist(outputFileName)
            h5create(outputFileName, '/stains', size(stains), 'Datatype', 'single', 'ChunkSize', chunk_size, 'Deflate', 9);
        end
        h5write(outputFileName, '/stains', single(stains));
    
end
