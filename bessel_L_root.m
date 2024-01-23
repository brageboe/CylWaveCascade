function [roots] = bessel_L_root(m, amount_roots, r_o, r_i, options)
% Returns the first [amount_roots] of zeroes of the
% superposition-bessel-function L of order [m]. The function L depends on
% the ratio of outer and inner waveguide radii, [r_o] and [r_i].
    arguments 
        m = 0;                  % mode 
        amount_roots = 2;       % amount of roots to find
        r_o = 6;                % outer waveguide radius
        r_i = 3;                % inner waveguide radius
        options.plot_figure = false;
        options.max_root_value = 1e4;
        options.offset = 1e-6;  % search from (0+error) and up, to avoid the inf at zero
        options.dp = 1e5;       % #datapoints
    end

    A = r_o/r_i;
    if A == 1
        error("Inner and outer conductors cannot have equal radii.")
    end

    % the transcendental equation that determines the roots: func1 = func2.
    func1 = @(x) besselj(m, x.*A).*bessely(m, x);
    func2 = @(x) besselj(m, x).*bessely(m, x.*A);

    %% Estimate root values
    x_dp = linspace(0+options.offset, options.max_root_value, options.dp);
    func_difference = func1(x_dp) - func2(x_dp);
    
    root_estimate_value = [];
    for i = 1:(options.dp-1)
        if sign(func_difference(i)) ~= sign(func_difference(i+1))
            root_estimate_value(end+1) = x_dp(i);
        end
        if length(root_estimate_value) == amount_roots % break early if we have #roots we want
            break
        end
    end
    if length(root_estimate_value) < amount_roots
        warning("Loop reached end without finding the requested amount of roots.");            
    end

    %% Calculate root values accurately
    roots = zeros(1, length(root_estimate_value));
    for i = 1:length(root_estimate_value)
        roots(i) = fzero(@(x) func1(x) - func2(x), root_estimate_value(i));
    end

    %% Figure
    if options.plot_figure
        %x = linspace(0, options.max_root_value, options.dp);
        figure()
        hold on
        plot(x_dp, func1(x_dp), 'Linewidth', 1.5);
        plot(x_dp, func2(x_dp), 'Linewidth', 1.5);
        plot(x_dp, func_difference, 'k:', 'Linewidth', 1.5);
        yline(0)
        title("$r_0="+r_o+"; r_i="+r_i+"$",'Interpreter', 'latex')
        legend("$J_"+m+"(xr_o/r_i)Y_"+m+"(x)$", ...
            "$J_"+m+"(x)Y_"+m+"(xr_o/r_i)$", ...
            "$J_"+m+"(xr_o/r_i)Y_"+m+"(x)-J_"+m+"(x)Y_"+m+"(xr_o/r_i)$", ...
            'Interpreter', 'latex', 'location','SE')
        xlabel("$x$", 'Interpreter', 'latex')
        grid on
        hold off
    end
end