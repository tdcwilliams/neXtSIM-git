%%define location of /Data/sim:
if exist('/Volumes/sim')
   %% johansen
   data_sim = '/Volumes/sim'
elseif exist('/Volumes/Tim_Ext_HD2/neXtSIM/data')
   %% external hard drive
   data_sim = '/Volumes/Tim_Ext_HD2/neXtSIM/data'
else
   disp('Add paths to neXtSIM data');
   error('- eg connect to johansen (/Data/sim) with cmd+k');
end

%% define paths;
bamg_path      = '../../ISSM-trunk-jpl-svn/lib';
cholmod2_paths = [data_sim,'/MATLAB/SuiteSparse/CHOLMOD/MATLAB'];
johansen_paths = [data_sim,'/data'];%%+all subdirs
nextsim_path   = '../../neXtSim-trunk-Sourcetree';
mmap_path      = '~/Dropbox/MATHS/programs/matlab/custom-packages/m_map';
%% add paths
rmpaths;
dirs  = {nextsim_path,...
         [nextsim_path,'/code'],...
         [nextsim_path,'/tools'],...
         bamg_path,...
         cholmod2_paths,...
         mmap_path};
for loop_i=1:length(dirs)
   addpath(dirs{loop_i});
end

%% get all the tools,data:
%% - need all the sub-directories in these folders
lookin_dirs = {[nextsim_path,'/tools'],johansen_paths};

for j=1:length(lookin_dirs)

   lookin   = lookin_dirs{j};

   subdirs     = dir(lookin);
   isub        = [subdirs(:).isdir]; %# returns logical vector
   nameFolds   = {subdirs(isub).name}';

   nameFolds(ismember(nameFolds,{'.','..'})) = [];

   for loop_i=1:length(nameFolds)
      addpath([lookin '/' nameFolds{loop_i}]);
   end

end

showpaths;
