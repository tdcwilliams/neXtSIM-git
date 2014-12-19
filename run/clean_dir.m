!mkdir -p ../results

dirs  = dir('test_*_profile_results');
for j=1:length(dirs)
   pdir  = dirs(1).name;
   sss   = strsplit(pdir,'_');
   cnum  = sss{2};
   num   = str2num(cnum);
   odir  = ['test_',cnum,'_outputs'];
   %%
   newdir   = ['../results/test_',sss{2}];
   eval(['!mkdir -p ',newdir]);
   eval(['!mv ',pdir,' ',newdir]);
   eval(['!mv ',odir,'/pdf ',newdir]);
   eval(['!mv ',odir,' ',newdir]);
end
