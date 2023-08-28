 function [accuracy sensitivity specificity precision recall f_measure gmean] = Evaluate(ACTUAL,PREDICTED,pos_class)
idx = (ACTUAL()==pos_class);

p = length(ACTUAL(idx));
n = length(ACTUAL(~idx));
N = p+n;

tp = sum(ACTUAL(idx)==PREDICTED(idx));
tn = sum(ACTUAL(~idx)==PREDICTED(~idx));
fp = n-tn;
fn = p-tp;

tp_rate = tp/p;
tn_rate = tn/n;

accuracy = (tp+tn)/N;
sensitivity = tp_rate;
specificity = tn_rate;
precision = tp/(tp+fp);
recall = sensitivity;
f_measure = 2*((precision*recall)/(precision + recall));
gmean = sqrt(precision*recall);

EVAL = [accuracy sensitivity specificity precision recall f_measure gmean];
