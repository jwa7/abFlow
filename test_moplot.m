%% Calculating the MOs

clear, clc, clf

% H2O molecule

atoms = [8 1 1];

charge = 0;

xyz_a0 = [0     0          0.12716; ... 
          0     0.758081  -0.50864; ... 
          0    -0.758081  -0.50864];
      
settings.method = 'RHF';
settings.basisset = 'cc-pVDZ';
settings.tolEnergy = 1e-8;
settings.tolDensity = 1e-8;

out = mocalc(atoms, xyz_a0, charge, settings);

%% Plotting the MOs

iMO = 0.1;
isolevel = 0.5;     % a0^(-3/2)
moplot(atoms, xyz_a0, out, iMO, isolevel);
