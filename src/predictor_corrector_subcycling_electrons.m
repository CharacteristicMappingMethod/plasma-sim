function [fs,params] = predictor_corrector_subcycling_electrons(params,fs)

    assert(params.Ns == 2, "only ions and electrons please");
    iT = params.it+1;
    
    s_ions = find(params.species_name=="ions");
    s_electrons = find(params.species_name=="electrons");

    % timesteps
    dt_ions = params.dt;            % big time step
    dt_e = params.dt/params.Nsubs;    % small time steps

    [Efield] = vPoisson(fs, params.grids, params.charge);
    Efield_s = Efield;
    
    f12 = fs;
    % make half step of ions:
    f12(:, :, s_ions) = Advect(f12(:, :, s_ions), params.charge(s_ions) / params.Mass(s_ions) * Efield, params.grids(s_ions), dt_ions / 2);

    for i =1:params.Nsubs
        f12(:, :, s_electrons) = Advect(f12(:, :, s_electrons), params.charge(s_electrons) / params.Mass(s_electrons) * Efield, params.grids(s_electrons), dt_e / 2);
        [Efield] = vPoisson(f12, params.grids, params.charge);
    end
    

    fs(:, :, s_ions) = Advect(fs(:, :, s_ions), params.charge(s_ions) / params.Mass(s_ions) * Efield, params.grids(s_ions), params.dt);

    % Recompute electric field
    for i =1:params.Nsubs
        fs(:, :, s_electrons) = Advect(fs(:, :, s_electrons), params.charge(s_electrons) / params.Mass(s_electrons) * Efield, params.grids(s_electrons), dt_e);
        [Efield] = vPoisson(f12, params.grids, params.charge);
    end
   
    % save electric field;
    params.Efield = Efield;
    params.Efield_list(:,iT) = Efield;

    
end