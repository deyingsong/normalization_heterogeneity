function [I,Ie,Ii,Ix] = atttRec_compute_input_current(Wrr, Wrf, Jrr, Jrf, r, r1, Ne, Ni, Nx, Ke, Ki, Kx)
    % Vectorized computation of input currents from excitatory and inhibitory neurons

    N = Ne + Ni;
    I = zeros(N, 1);  % Output current

    % Excitatory connections
    idx_e = 1:(Ne*Ke);                             % Indices for excitatory connections
    k_e = repelem((1:Ne)', Ke);                    % Presynaptic excitatory neuron indices
    postsyn_e = Wrr(idx_e);                        % Postsynaptic targets
    contrib_e = Jrr(idx_e) .* r(k_e);              % Contributions to I

    % Accumulate excitatory contributions
    Ie = accumarray(postsyn_e, contrib_e, [N, 1])*1e-3;

    % Inhibitory connections
    idx_i = (Ne*Ke+1):(Ne*Ke+Ni*Ki);               % Indices for inhibitory connections
    k_i = repelem((1:Ni)', Ki);                    % Presynaptic inhibitory neuron indices
    k_i = k_i + Ne;                                % Shift index to global neuron index
    postsyn_i = Wrr(idx_i);                        % Postsynaptic targets
    contrib_i = Jrr(idx_i) .* r(k_i);              % Contributions to I

    % Accumulate inhibitory contributions
    Ii = accumarray(postsyn_i, contrib_i, [N, 1])*1e-3;

    % Feedforward connections
    idx_x = 1:(Nx*Kx);                             % Indices for ffwd connections
    k_x = repelem((1:Nx)', Kx);                    % Presynaptic ffwd neuron indices
    postsyn_x = Wrf(idx_x);                        % Postsynaptic targets
    contrib_x = Jrf(idx_x) .* r1(k_x);              % Contributions to I

    % Accumulate ffwd contributions
    Ix = accumarray(postsyn_x, contrib_x, [N, 1])*1e-3;

    I = Ie + Ii + Ix;

end
