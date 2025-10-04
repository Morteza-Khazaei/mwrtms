function [depth, rho, Tsnow, dia, tau] = load_snowpack(filename)
% Load snowpack description file. 
% In the file, each row describes one layer, from top to bottom.
% Each row contains five number:
% layer_thickness (cm)	density (gm/cc)	Temperature (K)	grain_diameter (cm) stickiness
% 

fid = fopen(filename);
if fid < 0
    error('Could not open file');
end

data = fscanf(fid,'%f',[5 inf]);
fclose(fid);

depth   = data(1,:)';
rho     = data(2,:)';
Tsnow   = data(3,:)';
dia     = data(4,:)';
tau     = data(5,:)';

end
