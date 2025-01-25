function output = flag(x)
% x = 0, means small region, ewfd2
   % x = 333 means no need to distinguish them
   output = struct();
   if x==0
       output.dset = 'dset1';
       output.normE = 'ewfd2.normE';
       output.std = 'std1';
   else
       output.dset = 'dset2';
       output.normE = 'ewfd1.normE';      
       output.std = 'std2';
   end 

   if x == 333
        output.dset = 'dset1';
        output.normE = 'ewfd1.normE';      
        output.std = 'std1';
   end
end
