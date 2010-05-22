function kernType = sdlfmgpKernComposer(type, numOut, approx, options)

% SDLFMGPKERNCOMPOSER Composes kernel types from some generic options.
% FORMAT
% DESC composes a kernel type to pass to kernel create from a generic
% type and a few options.
% RETURN kernType : an array of cells which indicates the types of kernel
% ARG type : kernel type
% ARG numOut : number of outputs.
% ARG approx : indicates the type of approximation for the GP. 
% ARG options : options for kernel. If introduced it's assumed to be
% 'parametric'.
%
% SEE ALSO : multigpKernComposer.m 
%
% COPYRIGHT :  Mauricio A. Alvarez, 2010

% SEEALSO : sdlfmgpCreate

% SDLFMGP

switch type
    case 'sdsim'
        % Not implemented yet
    case 'sdlfm'
        numPositions = numOut/(1+options.includeVel + options.includeAccel);
        switch approx
            case 'ftc'
                kernType = cell(1, numOut+1);
                kernType{1} = 'multi';
                cont = 1;
            case {'dtc','fitc', 'pitc', 'dtcvar'}
                kernType = cell(1, numOut+2);
                kernType{1} = 'multi';
                kernType{2} = {'parametric', options.kern, 'sdrbf'};
                cont = 2;
        end
        for i = 1:numPositions
            cont = cont + 1;
            kernType{cont} = {'parametric', options.kern, 'sdlfm'};
        end
        if options.includeVel
            if options.includeVel~=options.kern.includeVel
                error('The include velocity option has to be the same in options and options.kern')
            end
            for i = 1:numPositions
                cont = cont + 1;
                kernType{cont} = {'parametric', options.kern, 'sdlfmv'};
            end
        end
        if options.includeAccel
            if options.includeAccel~=options.kern.includeAccel
                error('The include acceleration option has to be the same in options and options.kern')
            end
            for i = 1:numPositions
                cont = cont + 1;
                kernType{cont} = {'parametric', options.kern, 'sdlfma'};
            end
        end
end

