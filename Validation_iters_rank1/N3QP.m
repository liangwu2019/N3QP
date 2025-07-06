function [z,Iters,Ranks,Gap] = N3QP(H,h,epsilon,alpha,delta)
n = length(h);
z = zeros(n,1);
Iters = [];
Ranks = [];
Gap = [];
h_max_norm = max(abs(h));
if h_max_norm==0
    return
else
    ub = 1+delta;
    lb = 1.0/ub;
    sigma = sqrt(2)*delta*(1+delta)^2*alpha*sqrt((1+alpha)/(1-alpha))+0.5*(1+delta)^2*alpha^2/(1-alpha);
    beta = (alpha-sigma)/(1+alpha/sqrt(2*n));
    decrease_rate = 1.0-beta/sqrt(2*n);
    lambda = alpha/sqrt(2*n);
    H_hat = 2*lambda/h_max_norm * H;
    gamma = 1 - lambda/h_max_norm * h;
    theta = 1 + lambda/h_max_norm * h;
    phi = ones(n,1);
    psi = ones(n,1);
    gamma_hat = gamma;
    theta_hat = theta;
    phi_hat = phi;
    psi_hat = psi;
    tau = 1.0;
    M = inv(H_hat + 2*eye(n));
    d = gamma_hat./phi_hat + theta_hat./psi_hat;
    Idx = zeros(n,1);  
    iters = 0;
    ranks = 0;
    while(1)
        count = 0;
        for i=1:n
            ratio_gamma = gamma_hat(i)/gamma(i);
            ratio_theta = theta_hat(i)/theta(i);
            ratio_phi = phi_hat(i)/phi(i);
            ratio_psi = psi_hat(i)/psi(i);
            flag_gamma = (ratio_gamma<lb || ratio_gamma>ub);
            flag_theta = (ratio_theta<lb || ratio_theta>ub);
            flag_phi = (ratio_phi<lb || ratio_phi>ub);
            flag_psi = (ratio_psi<lb || ratio_psi>ub);
            if flag_gamma
                gamma_hat(i) = gamma(i);
            end
            if flag_theta
                theta_hat(i) = theta(i);
            end
            if flag_phi
                phi_hat(i) = phi(i);
            end
            if flag_psi
                psi_hat(i) = psi(i);
            end
            if(flag_gamma || flag_theta || flag_phi || flag_psi)
                count = count + 1;
                Idx(count) = i;
            end
        end
        temp1 = gamma_hat./phi_hat;
        temp2 = theta_hat./psi_hat;
        temp_r1 = (tau-gamma.*phi)./phi_hat;
        temp_r2 = (tau-theta.*psi)./psi_hat;
        delta_z_temp = temp_r2-temp_r1;
        for j=1:count
            i = Idx(j);
            coeff = -(temp1(i)+temp2(i)-d(i))/(1+(temp1(i)+temp2(i)-d(i))*M(i,i));
            M = M + coeff*M(:,i)*M(:,i)';
            d(i) = temp1(i)+temp2(i);
        end
        delta_z = M*delta_z_temp;        
        z = z + delta_z;
        gamma = gamma + temp_r1 + temp1.*delta_z;
        theta = theta + temp_r2 - temp2.*delta_z;
        phi = phi - delta_z;
        psi = psi + delta_z;
        tau = decrease_rate*tau;
        iters = iters + 1;
        Iters = [Iters; iters];
        ranks = ranks + count;
        Ranks = [Ranks; ranks];
        duality_gap = gamma_hat'*phi + theta'*psi;
        Gap = [Gap;duality_gap];
        if(duality_gap<=epsilon)    
            break
        end
    end
    return
end
