function [n, coeff, f, H] = wls_fir_filter(W_p, W_t, W_s, N, k, f_pass, f_stop,fs)
%   Inputs:
%       W_p        - Weight for the passband (positive scalar)
%       W_t        - Weight for the transition band (positive scalar)
%       W_s        - Weight for the stopband (positive scalar)
%       N          - Desired filter length (must be odd and positive)
%       k          - Number of frequency points (must be > (N-1)/2 + 1)
%       f_pass     - Passband cutoff frequency in Hz (must be < fs/2)
%       f_stop     - Stopband start frequency in Hz (must be > f_pass and < fs/2)
%
%   Outputs:
%       n          - Filter index vector (0 to N-1)
%       coeff      - Filter coefficients (1 x N vector)
%       f          - Frequency vector in Hz (1 x k vector)
%       H          - Frequency response (complex vector, 1 x k)
%
%   Example:
%       [n, coeff, f, H] = ls_fir_filter(100, 1, 1, 301, 500, 5000, 5500);

    %% 1. Initialize and Validate Inputs

    % Default Parameters (in case inputs are missing or invalid)
    default_W_p = 1;           % Default weight for passband
    default_W_t = 1;           % Default weight for transition band
    default_W_s = 1;           % Default weight for stopband
    default_N = 301;           % Default filter length (odd and positive)
    default_f_pass = 5000;     % Default passband cutoff frequency in Hz

    % Validate and Assign Weight for Passband W_p
    if  W_p <= 0
        warning('Invalid or missing passband weight W_p. Using default W_p = %d.', default_W_p);
        W_p = default_W_p;
    end

    % Validate and Assign Weight for Transition Band W_t
    if W_t <= 0
        warning('Invalid or missing transition band weight W_t. Using default W_t = %d.', default_W_t);
        W_t = default_W_t;
    end

    % Validate and Assign Weight for Stopband W_s
    if W_s <= 0
        warning('Invalid or missing stopband weight W_s. Using default W_s = %d.', default_W_s);
        W_s = default_W_s;
    end

    % Validate and Assign Filter Length N
    if N <= 0
        warning('Invalid or missing filter length N. Using default N = %d.', default_N);
        N = default_N;
    end
    if mod(N, 2) == 0
        warning('Filter length N is even. Incrementing N by 1 to make it odd.');
        N = N + 1;
    end

    L = (N-1)/2;               % Half-length for symmetric filter design
    n = 0:1:N - 1;               % Filter index vector

    % Validate and Assign Number of Frequency Points k
    if k <= 0 || k <= L + 1
        warning('Invalid or insufficient number of frequency points k. Using default k = %d.', (N-1)/2 + 1);
        k = (N - 1)/2 + 2;
    end

    % Validate and Assign Passband Cutoff Frequency f_pass
    if f_pass <= 0 || f_pass >= fs/2
        warning('Invalid or missing passband cutoff frequency f_pass. Using default f_pass = %d Hz.', default_f_pass);
        f_pass = default_f_pass;
    end

    % Validate and Assign Stopband Start Frequency f_stop
    if f_stop <= f_pass || f_stop >= fs/2
        warning('Invalid or missing stopband start frequency f_stop. Using default f_stop = %d Hz.', f_pass + 2000);
        f_stop = f_pass + 2000;
    end

    %% 2. Define Frequency Design Grid

    w = pi * linspace(0, 1, k); % Frequency vector from 0 to pi (rad/sample)

    %% 3. Desired Frequency Response (Low-Pass)

    Hd = zeros(1, k);                  % Initialize desired frequency response
    cut_off_pass = 2 * pi * f_pass / fs;   % Normalize passband cutoff frequency (rad/sample)
    cut_off_stop = 2 * pi * f_stop / fs;   % Normalize stopband start frequency (rad/sample)
    Hd(w <= cut_off_pass) = 1;               % Passband: 1 up to passband cutoff, 0 otherwise

    %% 4. Define Weights for Each Band

    % Initialize weight vector
    W = zeros(1, k);

    % Assign weights based on frequency bands
    W(w <= cut_off_pass) = W_p;                     % Passband
    W(w > cut_off_pass & w < cut_off_stop) = W_t;   % Transition band
    W(w >= cut_off_stop) = W_s;                     % Stopband

    % Construct diagonal weight matrix
    W_matrix = diag(W);                              % Diagonal matrix with weights

    %% 5. Construct the F Matrix for Weighted Least Squares

    % Construct cosine terms for the F matrix
    cos_terms = cos(w' * (1:L)); 

    % Assemble F matrix
    F = [ones(k,1), 2*cos_terms];

    %% 6. Solve the Weighted Least Squares Equation Precisely

    % Compute F' * W * F
    FWF = F' * W_matrix * F;

    % Compute F' * W * H
    FWH = F' * W_matrix * Hd';

    % Solve for h_half: (F'WF) h_half = F'WH
    h_half = FWF \ FWH;        
    h_half = h_half';           % Convert to row vector

    % h_half(1) = h(0), h_half(2) = h(1), ..., h_half(L+1) = h(L)

    %% 7. Construct Full Symmetric Filter Coefficients

    % Mirror the coefficients to ensure symmetry: h(k) = h(-k)
    h_full = [fliplr(h_half(2:end)), h_half]; % Length: N
    coeff = h_full;

    %% 8. Compute Frequency Response

    [H, f] = freqz(coeff, 1, k, fs);     % Frequency response and frequency vector in Hz
    
end
