function [alp, obj, model] = solve_svm(tra, par, yyK, alp)
    N = size(tra.X, 1);
    switch lower(par.opt)
        case 'libsvm'
            alp = zeros(N,1);
            H = tra.y*tra.y'.*yyK;
            H = H+H';
            f = -tra.y.*diag(yyK);
            lb = zeros(N,1);
            ub = par.C1 .* ones(N,1);
            Aeq = tra.y';
            beq = 1;
            opt = optimset('quadprog');
            opt.Algorithm = 'interior-point-convex';
            opt.Display = 'off';
            [alp,~, ~, output, ~] = quadprog(H, f, [], [], Aeq, beq, lb, ub,[],opt);
            if (isempty(alp))
                warning('No solution for the SVDD model could be found.');
                alp = zeros(N, 1);
                alp(1, 1) = 1;
            end
            alp = tra.y.*alp;
            alp(alp < par.C1*par.eps) = 0;
            alp(alp > par.C1*(1-par.eps)) = par.C1;
            obj = -0.5 * alp'* yyK * alp + alp'* diag(yyK);
            model.alpha = alp;        
            model.numIterations = output.iterations;
            model.supportVectorIndices = find(abs(alp) > par.alptolerance);
            model.boundary = model.supportVectorIndices(find((alp(model.supportVectorIndices) < ...
                ub(model.supportVectorIndices))&(alp(model.supportVectorIndices) > par.alptolerance)));
            if (size(model.boundary, 1) < 1)
                model.boundary = model.supportVectorIndices;
            end
%             alp(find(abs(alp) < par.alptolerance)) = 0;
            model.supportVectors = tra.X(model.supportVectorIndices, :);
            model.supportVectorAlpha = alp(model.supportVectorIndices);
            model.numSupportVectors = size(model.supportVectorIndices, 1);
            tmp_ = -2*sum((ones(N, 1)*alp').*yyK, 2);
            model.offset = sum(sum((alp*alp').*yyK));
            model.radius = sqrt(mean(diag(yyK))+model.offset+mean(tmp_(model.boundary, :)));
    end
end