function bool = above_cutoff(n, d_z, i)
    arguments
        n   % mode number
        d_z % longitudinal wavenumbers; left-side stored in d_z(:,1), right-side in d_z(:,2)
        i   % 1 or 2, i.e. left-side or right-side of junction
    end
    if real(d_z(n,i)) == 0
        bool = false;
        return
    elseif real(d_z(n,i)) > 0
        bool = true;
        return
    else
        bool = -1;
        error("Negative longitudinal wavenumber."); % Should never happen.
    end
end