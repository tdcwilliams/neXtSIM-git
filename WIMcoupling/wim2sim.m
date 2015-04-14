function [simul_out,gridprams]   = wim2sim(simul_out,gridprams)
%% function to call from perform_simul.m (called at end of EB_model_v2.m)
%% do it round "Time interpolation of the forcings"
%%    - when Vair and Voce are calculated
%% - wave stresses should go into system_assemble_mex.c at "Step1"
%%    - similar to Voce, Vair?
%%               coef_Vair=Vair_factor?
%% - Dmax should maybe go into Step1 too (influence damage?)
%%    - also thermodynamic effect (lat melt - thermo_ow_mex.c)

if ~exist('gridprams','var')
   %% Get WIM grid from file
   %% NB needs to be a regular grid in x,y (stere proj coords)
   gitdir      = getenv('GIT_REPOS');
   gdir        = [gitdir,'/WIM2d/fortran/run/inputs/'];%%directory with grid files
   gridprams   = fn_get_grid(gdir);
end

if ~exist('simul_out','var')
   testdir     = 'test_inputs';
   testfile    = [testdir,'/simul_out_squaresmall1km_test2_step0.mat'];
   tf          = load(testfile);
   simul_out   = tf.simul_out;
   clear tf;
end

if ~exist('mesh','var')
   %[mesh, element, simul_in.ind_node_fix_bnd, simul_in.ind_node_free_bnd, simul_in.ind_element_free_bnd] =... 
   [mesh,element] =...
      importbamg(simul_out.bamg.mesh, simul_out.bamg.geom);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Interpolate needed stuff from mesh to WIM grid.
[xvert,yvert,xcent,ycent] = get_centres(simul_out,mesh,element);
% [x,y]vert: vertices of FEM mesh [km]
% [x,y]cent: centres of FEM mesh  [km]
Nn = length(xvert)%%number of nodes
Ne = length(xcent)%%number of elements

%%data on FEM mesh (at centres)
data        = zeros(length(xcent),3);%%1 col for each field
data(:,1)   = simul_out.c;           %% conc
data(:,2)   = simul_out.h;           %% thickness
if isfield(simul_out,'Dmax')
   INIT_DMAX   = 0;
   data(:,3)   = simul_out.Dmax;
else
   INIT_DMAX   = 1;
   data(:,3)   = [];
end

if 0
   % from ISSM, use mex from InterpFromMeshToGrid.cpp
   % to interp from mesh to WIM grid
   % TODO get working

   %WIM grid info
   xmin     = 1e3*min(gridprams.X(:)); % WIM grid xmin (stere proj, km)
   ymax     = 1e3*max(gridprams.Y(:)); % WIM grid ymax (stere proj, km)
   xposting = 1e3*gridprams.dx;        % res in x dirn (km)
   yposting = 1e3*gridprams.dy;        % res in y dirn (km)
   nlines   = gridprams.nx;            % no of rows in WIM grid
   ncols    = gridprams.ny;            % no of cols in WIM grid

   %% global variable element
   index = element.num_node(:,[1 3 2]);

   %%call interpolation routine
   missing_m2g          = -1000.; % missing value (if WIM grid pt is out of FEM mesh)
   size(index)
   size(xcent)
   size(ycent)
   size(data)
                     %xmin,ymax,xposting,yposting,nlines,ncols,missing_m2g
   [xWIM,yWIM,griddata] = ...
      InterpFromMeshToGrid(index,xcent,ycent,data,xmin,ymax,...
                           xposting,yposting,nlines,ncols,missing_m2g);
   %[x_m y_m data_grid]=InterpFromMeshToGrid(elements,x,y,data,xlim(1),ylim(2),
   %                        post,post,round(diff(ylim)/post),round(diff(xlim)/post),NaN);
   %[x_m,y_m,griddata]=InterpFromMeshToGrid(index,x,y,data,xmin,ymax,xposting,yposting,nlines,ncols,default_value)
   % xWIM = x points of WIM grid (already have)
   % yWIM = y points of WIM grid (already have)
   % griddata: columns are fields on WIM grid
else
   % use interp2
   griddata = 0*data;
   gx       = gridprams.X(:,1);
   gy       = gridprams.Y(1,:).';
   for j=1:size(data,2)
      griddata(:,j)  = interp2(xcent,ycent,data(:,j),gx,gy);
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set conc, thickness, Dmax on the WIM grid.
ice_fields  = struct('cice'      ,0*gridprams.X,...
                     'hice'      ,0*gridprams.X,...
                     'Dmax'      ,0*gridprams.X,...
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [node_x,node_y,xc,yc] = get_centres(simul_out,mesh,element);

% Adding displacement
node_x = mesh.node.x' + simul_out.UM(1:2:end)*1e-3;% km: size(Nn,1) Nn=no of nodes
node_y = mesh.node.y' + simul_out.UM(2:2:end)*1e-3;% km: (Nn,1) Nn=no of nodes

% Compute the position of the 3 nodes
xy_tricorner(:,:,1) = node_x(element.num_node)*1000;%m: (Ne,3,1) Ne=no of elements
xy_tricorner(:,:,2) = node_y(element.num_node)*1000;%m: (Ne,3,1) Ne=no of elements

% Compute position of the center
xc = mean(xy_tricorner(:,:,1),2)/1000;%km: (Ne,1)
yc = mean(xy_tricorner(:,:,2),2)/1000;%km: (Ne,1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
