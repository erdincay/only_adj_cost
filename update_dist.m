function MuNew = update_dist(Mu, pi, pol, Astate, Kt_ind_hi, Kt_ind_lo, weight_hi)

pol_fl = floor(pol);
pol_ce = ceil(pol);

T_z = pi.z(:,:,Astate);

MuNew = zeros(size(Mu));

[Mu_row, Mu_col] = find(Mu>0);

    for i = 1:length(Mu_row)
        kprime_ind_lo_fl = pol_fl(Mu_row(i), Mu_col(i), Astate, Kt_ind_lo);
        kprime_ind_lo_ce = pol_ce(Mu_row(i), Mu_col(i), Astate, Kt_ind_lo);
        kprime_ind_hi_fl = pol_fl(Mu_row(i), Mu_col(i), Astate, Kt_ind_hi);
        kprime_ind_hi_ce = pol_ce(Mu_row(i), Mu_col(i), Astate, Kt_ind_hi);

        if (kprime_ind_lo_fl == kprime_ind_lo_ce) || kprime_ind_lo_ce == 1
            MuNew(:,kprime_ind_lo_ce) = MuNew(:,kprime_ind_lo_ce) + ...
                (1-weight_hi) * (Mu(Mu_row(i), Mu_col(i)) * T_z(Mu_row(i),:)');
        else
            weight_kprime_hi = pol(Mu_row(i), Mu_col(i), Astate, Kt_ind_lo) - ...
                pol_fl(Mu_row(i), Mu_col(i), Astate, Kt_ind_lo);
            MuNew(:,kprime_ind_lo_fl:kprime_ind_lo_ce) = MuNew(:,kprime_ind_lo_fl:kprime_ind_lo_ce) + ...
                (1-weight_hi) * (Mu(Mu_row(i), Mu_col(i)) * T_z(Mu_row(i),:)') * ...
                [1-weight_kprime_hi weight_kprime_hi];
        end

        if (kprime_ind_hi_fl == kprime_ind_hi_ce) || kprime_ind_hi_ce == 1
            MuNew(:,kprime_ind_hi_ce) = MuNew(:,kprime_ind_hi_ce) + ...
                weight_hi * (Mu(Mu_row(i), Mu_col(i)) * T_z(Mu_row(i),:)');
        else
            weight_kprime_hi = pol(Mu_row(i), Mu_col(i), Astate, Kt_ind_hi) - ...
                pol_fl(Mu_row(i), Mu_col(i), Astate, Kt_ind_hi);
            MuNew(:,kprime_ind_hi_fl:kprime_ind_hi_ce) = MuNew(:,kprime_ind_hi_fl:kprime_ind_hi_ce) + ...
                weight_hi * (Mu(Mu_row(i), Mu_col(i)) * T_z(Mu_row(i),:)') * ...
                [1-weight_kprime_hi weight_kprime_hi];
        end
    end