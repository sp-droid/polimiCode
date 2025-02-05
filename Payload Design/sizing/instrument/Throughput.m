function [KIPS_a, RAM_a, ROM_a] = Throughput(N, KIPS_typ, Freq, Freq_typ, Code_tot, bpw, Data_tot)


    % Constants
    conv = 8 * 1024 * 1024; % Conversion bytes to MB

    % Element-wise calculations
    KIPS = (N .* KIPS_typ .* Freq) ./ Freq_typ; % KIPS calculation
    RAM = (Code_tot .* bpw) / conv;             % RAM calculation
    ROM = ((Code_tot + Data_tot) .* bpw) / conv; % ROM calculation

    % Summing results
    KIPS_a = sum(KIPS);
    RAM_a = sum(RAM);
    ROM_a = sum(ROM);

    % Display results
    fprintf('Throughput (KIPS): %.4f\n', KIPS_a);
    fprintf('Total RAM (MB): %.4f\n', RAM_a);
    fprintf('Total ROM (MB): %.4f\n', ROM_a);
end


