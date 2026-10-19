% Independent M5 cross-check using VoxScatter's DDA implementation.
% Usage: octave --quiet tools/voxscatter_m5_reference.m VOXSCATTER_ROOT
args = argv();
if numel(args) ~= 1
  error('usage: voxscatter_m5_reference.m VOXSCATTER_ROOT');
end
addpath(genpath(fullfile(args{1}, 'src')));

epsilon_r = 3.0;
refractive_index = sqrt(epsilon_r);
size_parameter = 0.8;
voxels_per_internal_wavelength = 8;
[r, idx, resolution, projected_area, lambda_ext, lambda_int] = ...
  geometry('sphere', refractive_index, size_parameter, voxels_per_internal_wavelength);
[L, M, N, unused] = size(r);
n_domain = L*M*N;
k0 = 2*pi/lambda_ext;
incident_polarization = [0, 0, 1];
incident_direction = [1, 0, 0];
kvector = k0*incident_direction;
incident = PlaneWavePointWise(incident_polarization, kvector, r);
[operator, toeplitz, alpha_ldr] = ...
  getOPERATOR_DDA(r, k0, refractive_index, kvector, incident_polarization, 'off');
identity = ones(L, M, N);
mask = zeros(L, M, N);
mask(idx) = 1;
active = [idx; n_domain+idx; 2*n_domain+idx];
voxel_volume = resolution^3;
apply = @(polarization) MVP(polarization, operator, identity, ...
  voxel_volume*mask, 1/alpha_ldr, 'notransp', active, 0);
[solution, flag, relative_residual] = bicgstab(apply, incident(active), 1e-11, 2000);
physical_dipole = 4*pi*sum(reshape(solution, [], 3), 1)*voxel_volume;
fprintf('VOXSCATTER_CELLS %d\n', numel(idx));
fprintf('VOXSCATTER_FLAG %d\n', flag);
fprintf('VOXSCATTER_RESIDUAL %.17g\n', relative_residual);
fprintf('VOXSCATTER_NORMALIZED_DIPOLE %.17g %.17g %.17g %.17g %.17g %.17g\n', ...
  real(physical_dipole(1)), imag(physical_dipole(1)), ...
  real(physical_dipole(2)), imag(physical_dipole(2)), ...
  real(physical_dipole(3)), imag(physical_dipole(3)));
