function [alp_latter, obj, model,aa,bb] = solve_svm2(tra, par, yyK, alp, alp_former)
N = size(tra.X, 1);
beta = zeros(size(tra.y,1),1);
rr = alp_former'*yyK*alp_former + sum(sum(alp_former.*yyK.*beta')) + 1/4*beta'*yyK*beta;
dist_down = (diag(yyK) - 2*sum(alp_former.*yyK,2) + rr).^(1/2) - (diag(yyK)'*beta/2 + rr - alp_former'*yyK*alp_former).^(1/2);
dist_up = (diag(yyK) - 2*sum(alp_former.*yyK,2) + rr).^(1/2) + (diag(yyK)'*beta/2 + rr - alp_former'*yyK*alp_former).^(1/2);
du1 = sort(dist_down);
R_up = du1(ceil(1/par.C2));
if (floor((1/par.C2))<=0)
    R_down = -1000;
else
    R_down = du1(floor(1/par.C2));
end
alp(dist_down > R_up) = par.C2;
alp(dist_up < R_down) = 0;
R = alp == Inf;
R1 = alp == 0;
R2 = alp == par.C2;
RR = R1|R2;
aa =sum(RR);
bb = aa/N;
if sum(R) == 0
    alp_latter = alp;
else
    yyK1 = yyK(R,R);
    H = tra.y(R)*tra.y(R)'.*yyK1;
    H1 = tra.y*tra.y'.*yyK;
    H = H+H';
    H1 = H1+H1';
    F = -tra.y(R).*diag(yyK1) + 2*H1(R,RR)*alp_former(RR);
    lb = zeros(size(F));
    ub = par.C2.* ones(size(F));
    Aeq = tra.y(R)';
    beq = 1 - tra.y(RR)'*alp_former(RR);
    opt = optimset('quadprog');
    opt.Algorithm = 'interior-point-convex';
    opt.Display = 'off';
    [alp_yu,~, ~, output, ~] = quadprog(2*H, F, [], [], Aeq, beq, lb, ub,[],opt);
    if (isempty(alp_yu))
        warning('No solution for the SVDD model could be found.');
        alp_yu = zeros(size(alp_yu, 1),1);
        alp_yu(1, 1) = 1;
    end
    alp(R) = alp_yu;
    alp(alp < par.C2*1e-3) = 0;
    alp(alp > par.C2*(1-1e-3)) = par.C2;
    alp_latter = alp;
end
obj = -0.5 * alp_latter'* yyK * alp_latter + alp_latter'* diag(yyK);
model.alpha = alp_latter;
model.supportVectorIndices = find(abs(alp) > par.alptolerance);
model.boundary = model.supportVectorIndices;
%             alp(find(abs(alp) < par.alptolerance)) = 0;
model.supportVectors = tra.X(model.supportVectorIndices, :);
model.supportVectorAlpha = alp_latter(model.supportVectorIndices);
model.numSupportVectors = size(model.supportVectorIndices, 1);
tmp_ = -2*sum((ones(N, 1)*alp_latter').*yyK, 2);
model.offset = sum(sum((alp_latter*alp_latter').*yyK));
model.radius = sqrt(mean(diag(yyK))+model.offset+mean(tmp_(model.boundary, :)));
end
