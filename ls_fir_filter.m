function [n, coeff, f, H] = ls_fir_filter(N, k, f_cutoff,fs)
%   Inputs:
%       N          - Desired filter length (must be odd and positive)
%       k          - Number of frequency points (must be > (N-1)/2 + 1)
%       f_cutoff   - Desired cutoff frequency in Hz (must be < fs/2)
%
%   Outputs:
%       n          - Filter index vector (0 to N-1)
%       coeff      - Filter coefficients (1 x N vector)
%       f          - Frequency vector in Hz (1 x k vector)
%       H          - Frequency response (complex vector, 1 x k)
%
%   Example:
%       [n, coeff, f, H] = ls_fir_filter(11, 50, 2000);
    %% 1. Initialize and Validate Inputs

    % Default Parameters
    default_N = 15;            % Default filter length (odd and positive)
    default_f_cutoff = 5000;   % Default cutoff frequency in Hz

    % Validate and Assign Filter Length N
    if N <= 0
        warning('Invalid filter length N. Using default N = %d.', default_N);
        N = default_N;
    end
    if mod(N, 2) == 0
        warning('Filter length N is even. Incrementing N by 1 to make it odd.');
        N = N + 1;
    end

    L = (N-1)/2;               % Half-length for symmetric filter design
    n = 0:1:N-1;               % Filter index vector

    % Validate and Assign Number of Frequency Points k
    if k <= 0 || k <= (L + 1)
        warning('Invalid or insufficient number of frequency points k. Using default k = %d.', (N - 1)/2 + 2);
        k = (N - 1)/2 + 2;
    end

    % Validate and Assign Cutoff Frequency f_cutoff
    if f_cutoff <= 0 || f_cutoff >= fs/2
        warning('Invalid or missing cutoff frequency f_cutoff. Using default f_cutoff = %d Hz.', default_f_cutoff);
        f_cutoff = default_f_cutoff;
    end

    %% 2. Define Frequency Design Grid

    w = pi * linspace(0, 1, k); % Frequency vector from 0 to pi (rad/sample)

    %% 3. Desired Frequency Response (Low-Pass)

    Hd = zeros(1, k);                  % Initialize desired frequency response
    cut_off = 2 * pi * f_cutoff / fs;  % Normalize cutoff frequency (rad/sample)
    Hd(w <= cut_off) = 1;               % Passband: 1 up to cutoff frequency, 0 otherwise

    %% 4. Construct the F Matrix for Least Squares

    % Construct cosine terms for the F matrix
    cos_terms = cos(w' * (1:L)); 

    % Assemble F matrix
    F = [ones(k,1), 2*cos_terms]; 

    %% 5. Solve the Least Squares Problem

    % Solve F * h_half' = Hd' for h_half
    h_half = F \ Hd';            % Least squares solution (column vector)
    h_half = h_half';            % Convert to row vector

    % h_half(1) = h(0), h_half(2) = h(1), ..., h_half(L+1) = h(L)

    %% 6. Construct Full Symmetric Filter Coefficients

    % Mirror the coefficients to ensure symmetry: h(k) = h(-k)
    h_full = [fliplr(h_half(2:end)), h_half]; % Length: N
    coeff = h_full;

    %% 7. Compute Frequency Response
    [H, f] = freqz(coeff, 1, k, fs);     % Frequency response and frequency vector in Hz

end
