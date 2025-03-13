function res = profile_check(mean_U,zpos,Reynolds_number)

    log_law = log(zpos*Reynolds_number)/0.41+5;
    profile_error_vec = mean_U-log_law;
    profile_error_vec(profile_error_vec == Inf) = [];
    len = find(zpos*Reynolds_number>30,3,'first');

    profile_error = min(abs(profile_error_vec(len)));
    if(profile_error>1)
        res = true;
    else
        res = false;
    end
end