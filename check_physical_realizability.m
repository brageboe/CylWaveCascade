function bool = check_physical_realizability(S, error, options)
    % Checks if the sum of each column of S is physical (has value 1),
    % within a certain margin of error.
    %
    % If unphysical: a possible reason (if not due to fault in code) is
    % that operating frequency is below cutoff.
    arguments
        S
        error = 1e-3;
        options.print_warning = false
    end
    N = length(S); % [2 times #modes]
    for j=1:N
        test = sum(abs(S(:,j)).^2);
        if abs(test-1) > error
            if j > 1 
                test0 = 0;
                for k=(j:-1:2) % sum over all modes lower than j
                    test0 = test0 + sum(abs(S(:,j-1)).^2);
                end
                if abs(test0-1) < error
                    % if mode j fails but modes sub j succeeds, then we are
                    % likely below the cutoff of mode j but above cutoff of 1:(j-1). 
                    % We don't want to spam warning msg for these cases.
                    bool = true;
                    return
                end
            end
            bool = false;
            if options.print_warning
                %warning("Unphysical S-matrix (column "+j+")");
                warning("sum(abs(S(:,j)).^2) = "+test+", j = "+j);
            end
            return
        end
    end
    bool = true;
end