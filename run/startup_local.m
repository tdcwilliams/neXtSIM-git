%%define location of /Data/sim:
if exist('/Volumes/sim')
   %% johansen
   data_sim = '/Volumes/sim'
elseif exist('/Volumes/Tim_Ext_HD2/WORK/neXtSIM')
   %% external hard drive
   data_sim = '/Volumes/Tim_Ext_HD2/WORK/neXtSIM'
else
   disp('Add paths to neXtSIM data');
   error('- eg connect to johansen (/Data/sim) with cmd+k');
end

%% define paths;
bamg_path      = '../../ISSM-trunk-jpl-svn/lib';
johansen_paths = [data_sim,'/data'];%%+all subdirs
nextsim_path   = '../../neXtSim-trunk-Sourcetree';
topaz_path     = [johansen_paths,'/TOPAZ4/200709_201102'];
amsre_path     = [johansen_paths,'/AMSRE_ice_conc/2008/mar']

%% add paths
rmpaths;
dirs  = {nextsim_path,...
         [nextsim_path,'/code'],...
         [nextsim_path,'/tools'],...
         bamg_path,...
         topaz_path,...
         amsre_path};
for loop_i=1:length(dirs)
   addpath(dirs{loop_i});
end

%% other "MATLAB" paths
mat_path = [data_sim,'/MATLAB'];
addpath([mat_path]);
addpath([mat_path,'/age/']);
addpath([mat_path,'/defo_rgps/']);
addpath([mat_path,'/from An/']);
addpath([mat_path,'/SuiteSparse/CHOLMOD/MATLAB/']);
addpath([mat_path,'/m_map/']);

%% get all the tools,data:
%% - need all the sub-directories in these folders
lookin_dirs = {[nextsim_path,'/tools'],johansen_paths,[mat_path,'/m_map/']};

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
