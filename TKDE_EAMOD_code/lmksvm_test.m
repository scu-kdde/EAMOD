
function [out,K] = lmksvm_test(tes, mod, model2)
[N,~] = size(tes{1}.X);
    P = length(tes) - 1;
    for m = 1:P
        tes{m}.X = normalize_data(tes{m}.X, mod.nor.dat{m});
    end
    loc = locality(tes{P + 1}.X, mod.par.loc.typ);
    loc = normalize_data(loc, mod.nor.loc);
    out.eta = etas(loc, mod.gat, mod.par.eps, mod.par.gat.typ);
dis = zeros(N,1);
    for m = 1:P
        K = kernel(tes{m}, mod.sup{m}, mod.par.ker{m}, mod.par.nor.ker{m});
        K_ = kernel(tes{m}, tes{m}, mod.par.ker{m}, mod.par.ker{m});
        tmp = -2*sum((out.eta(:,m) * (mod.sup{m}.y.* mod.sup{m}.alp)').* K,2);
        dis = dis + out.eta(:, m).* out.eta(:, m) .* diag(K_) + tmp;
%         out.dis = out.eta(:, m).* out.eta(:, m) .* diag(K_) + out.dis - 2*sum((out.eta(:, m) * (mod.sup{m}.alp .* mod.sup{m}.y .* mod.sup{m}.eta)') .* K, 2);
    end
    dis2 = model2.offset * ones(N, 1);
    out.dis = sqrt(dis+dis2);
    out.predictedLabel = ones(N,1);
    out.predictedLabel(out.dis > model2.radius,1) = 2;
    out.numAlarm = sum(out.predictedLabel == -1,1);
end