function [simul_out,gridprams]   = wim2sim(simul_out,gridprams)
%% function to call from perform_simul.m (called at end of EB_model_v2.m)
%% do it round "Time interpolation of the forcings"
%%    - when Vair and Voce are calculated
%% - wave stresses should go into system_assemble_mex.c at "Step1"
%%    - similar to Voce, Vair?
%% - Dmax should maybe go into Step1 too (influence damage?)
%%    - also thermodynamic effect (lat melt - thermo_ow_mex.c)

if ~exist(gridprams,'var')
   %% Get WIM grid from file
   %% NB needs to be a regular grid in x,y (stere proj coords)
   gitdir      = getenv('GIT_REPOS');
   gdir        = [gitdir,'/WIM2d/run/inputs/'];%%directory with grid files
   gridprams   = fn_get_grid(gdir);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Interpolate needed stuff from mesh to WIM grid.
xcent = simul_out.bamg.mesh.Centers(:,1);    % = x centres of  FEM mesh [km] TODO check units!
ycent = simul_out.bamg.mesh.Centers(:,2);    % = y centres of  FEM mesh [km] TODO check units!
xvert = simul_out.bamg.mesh.Vertices(:,1);   % = x vertices of FEM mesh [km] TODO check units!
yvert = simul_out.bamg.mesh.Vertices(:,2);   % = y vertices of FEM mesh [km] TODO check units!

%%data on FEM mesh (at centres)
data        = zeros(length(x),3);%%1 col for each field
data(:,1)   = simul_out.c;       %% conc
data(:,2)   = simul_out.h;       %% thickness
if isfield(simul_out,'Dmax')
   INIT_DMAX   = 0;
   data(:,3)   = simul_out.Dmax;
else
   INIT_DMAX   = 1;
   data(:,3)   = [];
end

%WIM grid info
xmin     = min(gridprams.X(:));  % WIM grid xmin (stere proj, km)
ymax     = max(gridprams.Y(:));  % WIM grid ymax (stere proj, km)
xposting = gridprams.dx;         % res in x dirn (km)
yposting = gridprams.dy;         % res in y dirn (km)
nlines   = gridprams.nx;         % no of rows in WIM grid
ncol     = gridprams.ny;         % no of cols in WIM grid

%% global variable element
index = element.num_node(:,[1 3 2]);

%%call interpolation routine
missing_m2g          = -1000; % missing value (if WIM grid pt is out of FEM mesh)
[xWIM,yWIM,griddata] = ...
   InterpFromMeshToGrid(index,xcent,ycent,data,xmin,ymax,...
                        xposting,yposting,nlines,ncols,missing_m2g);
% xWIM = x points of WIM grid (already have)
% yWIM = y points of WIM grid (already have)
% griddata: columns are fields on WIM grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set conc, thickness, Dmax on the WIM grid.
ice_fields  = struct('cice'      ,0*gridprams.X,...
                     'hice'      ,0*gridprams.X,...
                     'Dmax'      ,[0*gridprams.X...
                     'ICE_MASK'  ,1+0*gridprams.X);

j_msg = find(griddata(:,1)==missing_m2g);
for j=1:3
   griddata(jmsg,j)  = 0;
end
ice_fields.ICE_MASK(jmsg)  = 0;
ice_fields.cice(:)         = griddata(:,1);
ice_fields.hice(:)         = griddata(:,2);
if INIT_DMAX==0
   ice_fields.dmax(:)   = griddata(:,3);
else
   ice_fields.dmax(:)   = 300*ice_fields.ICE_MASK;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define wave forcing on WIM grid
%% - could also interpolate from wave model to WIM grid
wave_fields = struct('Hs'        ,0*gridprams.X,...
                     'Tp'        ,0*gridprams.X,...
                     'mwd'       ,0*gridprams.X,...
                     'WAVE_MASK' ,0*gridprams.X);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Call WIM2d
out_fields  = run_WIM2d_io_mex(ice_fields,wave_fields);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Interpolate
%% taux,tauy
%% onto VERTICES of FEM grid
data        = zeros(length(xWIM),2);
data(:,1)   = out_fields.taux(:);
data(:,2)   = out_fields.tauy(:);
%%
missing_g2m = 0;%missing value - just set to 0 (one of FEM mesh points is out of WIM grid)
data_mesh   = InterpFromGridToMesh(xWIM,yWIM,data,xvert,yvert,missing_g2m);
%%
simul_out.taux_waves = data_mesh(:,1);
simul_out.tauy_waves = data_mesh(:,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Interpolate
%% Dmax
%% onto CENTRES of FEM grid
data(:,1)   = out_fields.Dmax(:);
data(:,2)   = [];
%%
missing_g2m = 0;%missing value - just set to 0 (one of FEM mesh points is out of WIM grid)
data_mesh   = InterpFromGridToMesh(xWIM,yWIM,data,xcent,ycent,missing_g2m);
%%
simul_out.Dmax = data_mesh(:,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
