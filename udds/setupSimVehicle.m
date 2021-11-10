% function results = setupSimVehicle
% 
%   Sets up and runs a simulation of an electric vehicle. The default
%   parameters roughly describe a first-generation GM/Chevy Volt during
%   charge-depletion mode, but are not validated against the actual
%   product. 
%
%   Inputs: [none]
%
%   Outputs:
%     results: A 1x4 cell array of simulation internal variables and
%              results, one for each of the simulated drive cycles. See
%              text for description of the variables in "results".

% Copyright (c) 2016 by Gregory L. Plett of 
% University of Colorado Colorado Springs (UCCS). 
%
% This work is licensed under a Creative Commons 
% Attribution-NonCommercial-ShareAlike 4.0 Intl. License, v. 1.0
%
% It is provided "as is", without express or implied warranty, for 
% educational and informational purposes only.
%
% This file is provided as a supplement to: Plett, Gregory L., "Battery
% Management Systems, Volume II, Equivalent-Circuit Methods," Artech House, 
% 2015.
function results = setupSimVehicle
  files = {'nycc.txt','udds.txt','us06.txt','hwfet.txt'}; % drive cycles

  % Setup the Chevy Volt-like vehicle
  % set up cell: capacity [Ah], weight [g], (vmax, vnom, vmin) [V]
  cell = setupCell(15,450,4.2,3.8,3.0);
  % set up module: numParallel, numSeries, overhead for this cell by weight
  module = setupModule(3,8,0.08,cell);
  % set up pack: numSeries, overhead by weight, (full SOC, empty SOC) [%],
  % efficiency for this module
  pack = setupPack(12,0.1,75,25,0.96,module);
  % set up motor: Lmax [Nm], (RPMrated, RPMmax) [RPM], efficiency,
  % inertia [kg/m2]
  motor = setupMotor(275,4000,12000,.95,0.2);
  % set up wheel: radius [m], inertia [kg/m2], rollCoef
  wheel = setupWheel(0.35,8,0.0111);
  % set up drivetrain: inverter efficiency, fractional regen torque limit,
  % gear ratio, gear inertia [kg/m2], gear efficiency, for this pack,
  % motor, and wheel 
  drivetrain = setupDrivetrain(0.94,0.9,12,0.05,0.97,pack,motor,wheel);
  % set up vehicle: # wheels, roadForce [N], Cd, frontal area [m2], weight
  % [kg], payload [kg], overhead power [W] for this drivetrain
  vehicle = setupVehicle(4,0,0.22,1.84,1425,75,200,drivetrain);
  
  fprintf('\n\nStarting sims...\n');
  for theCycle = 1:length(files)
    % Load the drive cycle: format = time (s), speed (mph)
    cycle = dlmread(files{theCycle},'\t',2,0);
    % Now, simulate the vehicle, store in "results" structure
    results{theCycle} = simVehicle(vehicle,cycle,0.3); %#ok<AGROW>
    % Estimate vehicle range
    range = (vehicle.drivetrain.pack.socFull - ...
             vehicle.drivetrain.pack.socEmpty) /...
            (vehicle.drivetrain.pack.socFull - ...
             results{theCycle}.batterySOC(end)) * ...
             results{theCycle}.distance(end);
    fprintf('Cycle = %s, range = %6.1f [km]\n',files{theCycle},range);
  end

  % Create a data structure to store cell-specific variables
  function cell = setupCell(capacity,weight,vmax,vnom,vmin)
    cell.capacity = capacity; % ampere hours
    cell.weight = weight; % grams
    cell.vmax = vmax; % volts
    cell.vnom = vnom; % volts
    cell.vmin = vmin; % volts
    cell.energy = vnom * capacity; % Watt-hours
    cell.specificEnergy = 1000*cell.capacity*cell.vnom/cell.weight; % Wh/kg
  end

  % Create a data structure to store module-specific variables
  function module = setupModule(numParallel,numSeries,overhead,cell)
    module.numParallel = numParallel;
    module.numSeries = numSeries;
    module.overhead = overhead;
    module.cell = cell;
    module.numCells = numParallel * numSeries;
    module.capacity = numParallel * cell.capacity;
    module.weight = module.numCells*cell.weight * 1/(1-overhead)/1000; % kg
    module.energy = module.numCells * cell.energy/1000; % kWh
    module.specificEnergy = 1000 * module.energy / module.weight; % Wh/kg
  end

  % Create a data structure to store pack-specific variables
  function pack = setupPack(numSeries,overhead,socFull,socEmpty,...
                            efficiency,module)
    pack.numSeries = numSeries;
    pack.overhead = overhead;
    pack.module = module;
    pack.socFull = socFull;
    pack.socEmpty = socEmpty; % unitless
    pack.efficiency = efficiency; % unitless, captures I*I*R losses
    pack.numCells = module.numCells * numSeries;
    pack.weight = module.weight * numSeries * 1/(1 - overhead); % kg
    pack.energy = module.energy * numSeries; % kWh
    pack.specificEnergy = 1000 * pack.energy / pack.weight; % Wh/kg
    pack.vmax = numSeries*module.numSeries*module.cell.vmax;
    pack.vnom = numSeries*module.numSeries*module.cell.vnom;
    pack.vmin = numSeries*module.numSeries*module.cell.vmin;
  end

  % Create a data structure to store motor-specific variables
  function motor = setupMotor(Lmax,RPMrated,RPMmax,efficiency,inertia)
    motor.Lmax = Lmax; % N-m
    motor.RPMrated = RPMrated;
    motor.RPMmax = RPMmax;
    motor.efficiency = efficiency;
    motor.inertia = inertia; %kg-m2
    motor.maxPower = 2*pi*Lmax*RPMrated/60000; % kW
  end

  % Create a data structure to store wheel-specific variables
  function wheel = setupWheel(radius,inertia,rollCoef)
    wheel.radius = radius; % m
    wheel.inertia = inertia; % km-m2
    wheel.rollCoef = rollCoef;
  end

  % Create a data structure to store drivetrain-specific variables
  function drivetrain = setupDrivetrain(inverterEfficiency,regenTorque,...
        gearRatio,gearInertia,gearEfficiency,pack,motor,wheel)
    drivetrain.inverterEfficiency = inverterEfficiency;
    drivetrain.regenTorque = regenTorque; % fraction of total torque avail.
    drivetrain.pack = pack;
    drivetrain.motor = motor;
    drivetrain.wheel = wheel;
    drivetrain.gearRatio = gearRatio;
    drivetrain.gearInertia = gearInertia; % km-m2, measured on motor side
    drivetrain.gearEfficiency = gearEfficiency;
    drivetrain.efficiency = pack.efficiency * inverterEfficiency * ...
                            motor.efficiency * gearEfficiency;
  end

  % Create a data structure to store all vehicle specifications
  function vehicle = setupVehicle(wheels,roadForce,Cd,A,weight,payload,...
                                  overheadPwr,drivetrain)
    vehicle.drivetrain = drivetrain;
    vehicle.wheels = wheels; % number of them
    vehicle.roadForce = roadForce; % N
    vehicle.Cd = Cd; % drag coeff
    vehicle.A = A; % frontal area, m2
    vehicle.weight = weight; % kg
    vehicle.payload = payload; % kg
    vehicle.overheadPwr = overheadPwr; % W
    vehicle.curbWeight = weight + drivetrain.pack.weight;
    vehicle.maxWeight = vehicle.curbWeight + payload;
    vehicle.rotWeight = ((drivetrain.motor.inertia + ...
                          drivetrain.gearInertia) * ...
                         drivetrain.gearRatio^2 + ...
                         drivetrain.wheel.inertia*wheels)/...
                         drivetrain.wheel.radius^2;
    vehicle.equivMass = vehicle.maxWeight + vehicle.rotWeight;
    vehicle.maxSpeed = 2 * pi * drivetrain.wheel.radius * ...
        drivetrain.motor.RPMmax * 60 / (1000*drivetrain.gearRatio); % km/h
  end
end