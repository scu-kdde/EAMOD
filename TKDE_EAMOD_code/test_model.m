function [a,output] = test_model(test_data,model,model2,kernels,c)
test.X=test_data;
 test_data = cell(1, 4);
 test_data{1} = test;
 test_data{2} = test;
 test_data{3} = test;
 test_data{4} = test;
parameters = lmksvm_parameter();
parameters.C = c;
parameters.ker = kernels;
parameters.nor.dat = {'true', 'true','true'};
parameters.nor.ker = {'true', 'true','true'};

output = lmksvm_test(test_data, model,model2);
a = output.predictedLabel;
end
