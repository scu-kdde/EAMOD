function [mod,model2,aa,bb,eta] = lmksvm_train2(tra, par,alp_former)
rand('twister', par.see); %#ok<RAND>
P = length(tra) - 1;
for m = 1:P
    mod.nor.dat{m} = mean_and_std(tra{m}.X, par.nor.dat{m});
    tra{m}.X = normalize_data(tra{m}.X, mod.nor.dat{m});
end
mod.loc = locality(tra{P + 1}.X, par.loc.typ);
mod.nor.loc = mean_and_std(mod.loc, par.nor.loc);
mod.loc = normalize_data(mod.loc, mod.nor.loc);
mod.gat = gating_initial(mod.loc, P, par.gat.typ);
eta = etas(mod.loc, mod.gat, par.eps, par.gat.typ);
N = size(tra{1}.X, 1);
yyKm = zeros(N, N, P);
for m = 1:P
    yyKm(:, :, m) = (tra{m}.y * tra{m}.y') .* kernel(tra{m}, tra{m}, par.ker{m}, par.nor.ker{m});
end
yyKeta = kernel_eta(yyKm, eta);
mod.yyKold = yyKeta;
alp = Inf(size(tra{1}.X,1),1);
[alp, obj, model2,aa,bb] = solve_svm2(tra{1}, par, yyKeta, alp, alp_former);
%display(sprintf('obj1:%10.6f', obj));
mod.obj = obj;
mod.sol = 1;
while 1 && P > 1
    oldObj = obj;
    [alp, eta, mod, obj, yyKeta,model2] = learn_eta(tra, par, yyKm, alp, eta, mod, obj, yyKeta,model2);
    %display(sprintf('obj: %10.6f', obj));
    mod.obj = [mod.obj, obj];
    if abs(obj - oldObj) <= par.eps * abs(oldObj)
        break;
    end
end
mod = rmfield(mod, 'loc');
for m = 1:P
    sup = find(alp .* eta(:, m) ~= 0);
    mod.sup{m}.ind = tra{m}.ind(sup);
    mod.sup{m}.X = tra{m}.X(sup, :);
    mod.sup{m}.y = tra{m}.y(sup);
    mod.sup{m}.alp = alp(sup);
    mod.sup{m}.eta = eta(sup, m);
end
sup = find(alp ~= 0);
act = find(alp ~= 0 & alp < par.C2);
if isempty(act) == 0
    mod.b =  mean(tra{1}.y(act) .* (1 - yyKeta(act, sup) * alp(sup)));
else
    mod.b = 0;
end
mod.par = par;
end
