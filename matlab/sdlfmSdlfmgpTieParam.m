function tieInd = sdlfmSdlfmgpTieParam(model, options)

% SDLFMSDLFMGPTIEPARAM Tie parameters for a sdlfmgp model with sdlfm kernel
% FORMAT
% DESC Tie the parameters for a sdlfmgp model that uses a sdlfm kernel.
% RETURN tieInd : cell with elements containing the indexes of parameters
% to tie.
% ARG model : model created
% ARG options : options for the model.
%
% COPYRIGHT : Mauricio A. Alvarez 2010

% SDLFMGP

tieInd = cell(model.nIntervals*(model.nlfPerInt+1),1);
% Tie the parameters for the inverse widths
cont = 0;
for i=1:model.nIntervals
    if model.nlfPerInt == 1
        indx = paramNameRegularExpressionLookup(model, ...
            ['inverse width interval ' num2str(i) '\.']); 
        cont = cont + 1;
        tieInd{cont} = indx;
    else
        for j=1:model.nlfPerInt
            indx = paramNameRegularExpressionLookup(model, ['.* inverse width ' num2str(j) '\.' ...
                ' interval ' num2str(i) '\.']); 
            cont = cont + 1;
            tieInd{cont} = indx;
        end
    end
end
% Tie the parameters of the switching times
for i=1:model.nIntervals
   indx = paramNameRegularExpressionLookup(model, ...
            ['.* switching point interval ' num2str(i) '\.']);  
   cont = cont + 1;
   tieInd{cont} = indx;
end

if options.includeVel || options.includeAccel
    massIndx = paramNameRegularExpressionLookup(model, '.* mass');
    springIndx = paramNameRegularExpressionLookup(model, '.* spring');
    damperIndx = paramNameRegularExpressionLookup(model, '.* damper');
    for i=1:model.numPositions
        if options.includeVel && ~options.includeAccel
            indxmass = [massIndx(i) massIndx(i+model.numPositions)];
            indxspring = [springIndx(i) springIndx(i+model.numPositions)];
            indxdamper = [damperIndx(i) damperIndx(i+model.numPositions)];        
        else
            indxmass = [massIndx(i) massIndx(i+model.nout) massIndx(i+2*model.numPositions)];
            indxspring = [springIndx(i) springIndx(i+model.nout) springIndx(i+2*model.numPositions)];
            indxdamper = [damperIndx(i) damperIndx(i+model.nout) damperIndx(i+2*model.numPositions)];
        end
        cont = cont + 1;
        tieInd{cont} = indxmass;
        cont = cont + 1;
        tieInd{cont} = indxspring;
        cont = cont + 1;
        tieInd{cont} = indxdamper;
    end
end


