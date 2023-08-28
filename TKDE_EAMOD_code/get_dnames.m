 clearvars -except a b res_all data_num path res_acc res_gmean; 
 if data_num==1
    
    dname = 'iris';
    dnames = [path dname];
end

if data_num==6
     
    dname = 'pima';
    dnames = [path dname];
end

if data_num==3
     
    dname = 'bupa';
    dnames = [path dname];
end


if data_num==7
     
    dname = 'germ';
    dnames = [path dname];    
end
if data_num==5
     
    dname = 'aust';
    dnames = [path dname];    
end

if data_num==4
     
    dname = 'japa';
    dnames = [path dname];    
end
if data_num==2
     
    dname = 'heart';
    dnames = [path dname];    
end
if data_num==8
     
    dname = 'park';
    dnames = [path dname];    
end
if data_num==9
     
    dname = 'abal';
    dnames = [path dname];    
end
