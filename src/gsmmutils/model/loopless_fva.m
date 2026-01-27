function loopless_fva(inputPath, franctionOfOptimum)

[modelDir, modelName, ext] = fileparts(inputPath);

model = readCbModel(inputPath);

disp(optimizeCbModel(model));

[minFlux, maxFlux] = fluxVariability(model, franctionOfOptimum, 'max', [model.rxns], 1, 'fastSNP');

numReactions = length(model.rxns);

for i = 1:numReactions
    model.lb(i) = minFlux(i); % Update the lower bound of the i-th reaction
    model.ub(i) = maxFlux(i); % Update the upper bound of the i-th reaction
end

disp(optimizeCbModel(model));

outputFileName = [modelName, '_fva', ext];
outputPath = fullfile(modelDir, outputFileName);

writeCbModel(modelConstrained, 'format', 'sbml', 'fileName', outputPath);


end