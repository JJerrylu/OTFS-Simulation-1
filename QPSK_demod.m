function S_demod = QPSK_demod(r_symbol)

    QPSK_constellation = [1+1i, -1+1i, -1-1i, 1-1i];
    bin_map = [0 0; 0 1; 1 1; 1 0];
    num_symbols = length(r_symbol);
    S_demod = zeros(1, num_symbols * 2);

    for i = 1:num_symbols
        [~, symbol_index] = min(abs(r_symbol(i) - QPSK_constellation));
        S_demod((i - 1) * 2 + 1:i * 2) = bin_map(symbol_index, :);
    end
end