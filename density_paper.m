% Raw data processing
close all; clear; clc; load('density_data.mat');
% Ice Mass Balance (IMB) buoy processing
bot_T66_cor = interp1(convertTo(t_T66,"datenum"),bot_int_T66,convertTo(t_fy,"datenum"),'pchip'); % Bottom interpolation for coring events
frz_T66 = T_T66; % temperature of IMB sesnors within ice; data from Lei et al. (2021), doi:10.1594/PANGAEA.938134
for i = 1:size(T_T66,1) % time
    for j = 1:size(T_T66,2) % depth
        if  j > (0.92 - bot_int_T66(i))*100/2
            frz_T66(i,j) = NaN;
        elseif  j < (0.92 - top_int_T66(i))*100/2
            frz_T66(i,j) = NaN;
        end
    end
end
frz_T62 = T_T62; % temperature of IMB sesnors within ice; data from Lei et al. (2022), doi:10.1594/PANGAEA.940231
for i = 1:size(T_T62,1) % time
    for j = 1:size(T_T62,2) % depth
        if  j > (0.92 - bot_int_T62(i))*100/2
            frz_T62(i,j) = NaN;
        elseif  j < (0.92 - top_int_T62(i))*100/2
            frz_T62(i,j) = NaN;
        end
    end
end
% Coring data processing; data from Oggier et al. (2023), doi:10.1594/PANGAEA.956732 and doi:10.1594/PANGAEA.959830
fb_sy = (h_sy_avg-d_sy_avg)/100; % freeboard
dS_fy = zS_fy; drho_fy = rho_fy; dT_fy = T_fy; dS_sy = S_sy; dT_sy = S_sy; drho_sy = S_sy; hs_sy = hd_sy(:,7); % depth
fbT_fy = hd_fy(:,1)-hd_fy(:,2); fbS_fy = hd_fy(:,3)-hd_fy(:,4); fbrho_fy = hd_fy(:,5)-hd_fy(:,6); % freeboard FYI
fbS_sy = hd_sy(:,3)-hd_sy(:,4); fbT_sy = hd_sy(:,1)-hd_sy(:,2); fbrho_sy = hd_sy(:,5)-hd_sy(:,6); % freeboard SYI
hS = zeros(length(S_fy),1); % preallocation
for i = 1:length(S_fy)
    hS(i) = zS_fy{i}(end) + 0.5*(zS_fy{i}(end)-zS_fy{i}(end-1)); % salinity core length (from S coordinates)
    dS_fy{i,1} = -zS_fy{i,1}*hd_fy(i,3)'/hS(i)/100+fbS_fy(i)/100; % depth of salinity measurements for FYLI
    drho_fy{i,1} = -zrho_fy{i,1}/100+fbrho_fy(i)/100; % depth of salinity measurements for FYLI
    dT_fy{i,1} = -zT_fy{i,1}/100+fbT_fy(i)/100; % depth of temperature measurements for FYLI
end
hS(22) = zS_fy{22}(end-3) + 0.5*(zS_fy{22}(end-3)-zS_fy{22}(end-1-3));
dS_fy{22,1} = -zS_fy{22,1}*hd_fy(22,3)'/hS(22)/100+fbS_fy(22)/100;
hS(21) = zS_fy{21}(end-2) + 0.5*(zS_fy{21}(end-2)-zS_fy{21}(end-1-2));
dS_fy{21,1} = -zS_fy{21,1}*hd_fy(21,3)'/hS(21)/100+fbS_fy(21)/100;
hS(20) = zS_fy{21}(end-1) + 0.5*(zS_fy{20}(end-1)-zS_fy{20}(end-1-1));
dS_fy{20,1} = -zS_fy{20,1}*hd_fy(20,3)'/hS(20)/100+fbS_fy(20)/100;
zT_sy{14,1}(3:29) = zT_sy{14,1}(3:29)*100;
for i = 1:length(S_sy); dT_sy{i,1} = -zT_sy{i,1}/100+fbT_sy(i)/100; end % depth of temperature measurements for SYI
for i = 1:length(S_sy); dS_sy{i,1} = -zS_sy{i,1}/100+fbS_sy(i)/100; end % depth of salinity measurements for SYI
for i = 1:length(S_sy); drho_sy{i,1} = -zrho_sy{i,1}/100+fbrho_sy(i)/100; end % depth of density measurements for SYI
S_fy_bulk = zeros(1,23); for i = 1:23; S_fy_bulk(i) = mean(S_fy{i}); end % bulk FYI salinity
S_sy_bulk = zeros(1,length(S_sy)); for i = 1:length(S_sy); S_sy_bulk(i) = mean(S_sy{i}); end % bulk SYI salinity
T_fy_bulk = zeros(1,23); for i = 1:23; T_fy_bulk(i) = mean(T_fy{i}(3:end),'omitnan'); end % bulk FYI temperature
T_sy_bulk = zeros(1,18); for i = 1:18; T_sy_bulk(i) = mean(T_sy{i}(3:end),'omitnan'); end % bulk SYI temperature
T_fy_b = zeros(1,23); for i = 1:23; T_fy_b(i) = mean(T_fy{i}(end),'omitnan'); end % Water temperature below FYI
% Density calculations
for i = 1:length(S_fy) % FYI density calculations
    T_fy{i} = min(-0.1,T_fy{i});
    T_lab_fy{i} = ones(size(Srho_fy{i})) * hd_fy(i,8); % Lab temperature
    F1_pr_rho_fy{i} = -4.732-22.45*T_lab_fy{i} - 0.6397*T_lab_fy{i}.^2 - 0.01074*T_lab_fy{i}.^3;
    F2_pr_rho_fy{i} = 8.903*10^-2 - 1.763*10^-2*T_lab_fy{i} - 5.33*10^-4*T_lab_fy{i}.^2 - 8.801*10^-6*T_lab_fy{i}.^3;
    vb_pr_rho_fy{i} = rho_fy{i} .* Srho_fy{i} ./ F1_pr_rho_fy{i}; % brine volume for T_lab
    rhoi_pr_fy{i} = (917-1.403*10^-1*T_lab_fy{i}); % pure ice density, Pounder (1965)
    vg_pr_fy{i} = max((1-rho_fy{i}.*(F1_pr_rho_fy{i}-rhoi_pr_fy{i}.*Srho_fy{i}/1000.*F2_pr_rho_fy{i})./(rhoi_pr_fy{i}.*F1_pr_rho_fy{i}))*1000,0,'includenan'); % gas volume for T_lab
    vg_pr_bulk_fy(i) = mean(vg_pr_fy{i},'omitnan'); vb_pr_bulk_fy(i) = mean(vb_pr_rho_fy{i},'omitnan');
    T_sal_fy{i} = interp1((dS_fy{i}(1)-dS_fy{i}(end))/(dT_fy{i}(3)-dT_fy{i}(end))*(dT_fy{i}(3:end)-dT_fy{i}(3))+dS_fy{i}(1),T_fy{i}(3:end),dS_fy{i}(1:end),'pchip');
    T_rho_fy{i} = interp1(dT_fy{i}(3:end),T_fy{i}(3:end),drho_fy{i},'pchip'); % in-situ temperature interpolation for density depth
    T_rho_fy{i} = min(-0.1,T_rho_fy{i});
    F3_pr_fy{i} = rhoi_pr_fy{i}.*Srho_fy{i}/1000./(F1_pr_rho_fy{i}-rhoi_pr_fy{i}.*Srho_fy{i}/1000.*F2_pr_rho_fy{i}); 
    rhoi_fy{i} = (917-1.403*10^-1*T_sal_fy{i}); % pure ice density, Pounder (1965) for T_insitu
    rhoi_rho_fy{i} = (917-1.403*10^-1*T_rho_fy{i}); % pure ice density, Pounder (1965) for T_insitu
    F1_fy{i} = -4.732-22.45*T_sal_fy{i} - 0.6397*T_sal_fy{i}.^2 - 0.01074*T_sal_fy{i}.^3; % Cox and Weeks (1983)  
    F1_fy{i}(T_sal_fy{i}>-2) = -4.1221*10^-2 + -1.8407*10^1*T_sal_fy{i}(T_sal_fy{i}>-2).^1 + 5.8402*10^-1*T_sal_fy{i}(T_sal_fy{i}>-2).^2 + 2.1454*10^-1*T_sal_fy{i}(T_sal_fy{i}>-2).^3; % F1 from Lepparanta and Manninen (1988)
    F1_rho_fy{i} = -4.732-22.45*T_rho_fy{i} - 0.6397*T_rho_fy{i}.^2 - 0.01074*T_rho_fy{i}.^3; % Cox and Weeks (1983)  
    F1_rho_fy{i}(T_rho_fy{i}>-2) = -4.1221*10^-2 + -1.8407*10^1*T_rho_fy{i}(T_rho_fy{i}>-2).^1 + 5.8402*10^-1*T_rho_fy{i}(T_rho_fy{i}>-2).^2 + 2.1454*10^-1*T_rho_fy{i}(T_rho_fy{i}>-2).^3; % F1 from Lepparanta and Manninen (1988)
    F2_fy{i} = 8.903*10^-2 - 1.763*10^-2*T_sal_fy{i} - 5.33*10^-4*T_sal_fy{i}.^2 - 8.801*10^-6*T_sal_fy{i}.^3; % Cox and Weeks (1983)
    F2_fy{i}(T_sal_fy{i}>-2) = 9.0312*10^-2 + -1.6111*10^-2*T_sal_fy{i}(T_sal_fy{i}>-2).^1 + 1.2291*10^-4*T_sal_fy{i}(T_sal_fy{i}>-2).^2 + 1.3603*10^-4*T_sal_fy{i}(T_sal_fy{i}>-2).^3; % F2 from Lepparanta and Manninen (1988)
    F2_rho_fy{i} = 8.903*10^-2 - 1.763*10^-2*T_rho_fy{i} - 5.33*10^-4*T_rho_fy{i}.^2 - 8.801*10^-6*T_rho_fy{i}.^3; % Cox and Weeks (1983)
    F2_rho_fy{i}(T_rho_fy{i}>-2) = 9.0312*10^-2 + -1.6111*10^-2*T_rho_fy{i}(T_rho_fy{i}>-2).^1 + 1.2291*10^-4*T_rho_fy{i}(T_rho_fy{i}>-2).^2 + 1.3603*10^-4*T_rho_fy{i}(T_rho_fy{i}>-2).^3; % F2 from Lepparanta and Manninen (1988)
    F3_rho_fy{i} = rhoi_rho_fy{i}.*Srho_fy{i}/1000./(F1_rho_fy{i}-rhoi_rho_fy{i}.*Srho_fy{i}/1000.*F2_rho_fy{i});
    vb_fy{i} = (1000).*(rhoi_fy{i}.*S_fy{i}/1000)./(F1_fy{i}-rhoi_fy{i}.*S_fy{i}/1000.*F2_fy{i}); % Brine volume for T_insitu, LM
    vb_rho_fy{i} = vb_pr_rho_fy{i} .* F1_pr_rho_fy{i} ./ F1_rho_fy{i}; % Brine volume for T_insitu, CW + LM
    vb_limit = 350;
    vb_bulk_fy(i) = mean(vb_fy{i}(vb_fy{i}>0 & vb_fy{i}<vb_limit),'includenan');
    rho_si_pr_fy{i} = (-vg_pr_fy{i}/1000+1).*rhoi_rho_fy{i}.*F1_rho_fy{i}./(F1_rho_fy{i}-rhoi_rho_fy{i}.*Srho_fy{i}/1000.*F2_rho_fy{i});
    rho_si_pr_bulk_fy(i) = mean(rho_si_pr_fy{i} (vb_rho_fy{i}>0 & vb_rho_fy{i}<vb_limit),'includenan');
end
for i = 20
    T_sal_fy{i} = interp1((dS_fy{i}(1)-dS_fy{i}(end-1))/(dT_fy{i}(3)-dT_fy{i}(end))*(dT_fy{i}(3:end)-dT_fy{i}(3))+dS_fy{i}(1),T_fy{i}(3:end),dS_fy{i}(1:end-1),'pchip');
    rhoi_fy{i} = (917-1.403*10^-1*T_sal_fy{i}); % pure ice density, Pounder (1965) for T_insitu
    F1_fy{i} = -4.732-22.45*T_sal_fy{i} - 0.6397*T_sal_fy{i}.^2 - 0.01074*T_sal_fy{i}.^3; % Cox and Weeks (1983)  
    F1_fy{i}(T_sal_fy{i}>-2) = -4.1221*10^-2 + -1.8407*10^1*T_sal_fy{i}(T_sal_fy{i}>-2).^1 + 5.8402*10^-1*T_sal_fy{i}(T_sal_fy{i}>-2).^2 + 2.1454*10^-1*T_sal_fy{i}(T_sal_fy{i}>-2).^3; % F1 from Lepparanta and Manninen (1988)
    F2_fy{i} = 8.903*10^-2 - 1.763*10^-2*T_sal_fy{i} - 5.33*10^-4*T_sal_fy{i}.^2 - 8.801*10^-6*T_sal_fy{i}.^3; % Cox and Weeks (1983)
    F2_fy{i}(T_sal_fy{i}>-2) = 9.0312*10^-2 + -1.6111*10^-2*T_sal_fy{i}(T_sal_fy{i}>-2).^1 + 1.2291*10^-4*T_sal_fy{i}(T_sal_fy{i}>-2).^2 + 1.3603*10^-4*T_sal_fy{i}(T_sal_fy{i}>-2).^3; % F2 from Lepparanta and Manninen (1988)
    vb_fy{i} = (1000).*(rhoi_fy{i}.*S_fy{i}(1:end-1)/1000)./(F1_fy{i}-rhoi_fy{i}.*S_fy{i}(1:end-1)/1000.*F2_fy{i}); % Brine volume for T_insitu, LM
    vb_bulk_fy(i) = mean(vb_fy{i}(vb_fy{i}>0 & vb_fy{i}<1000),'includenan');
end
for i = 21
    T_sal_fy{i} = interp1((dS_fy{i}(1)-dS_fy{i}(end-2))/(dT_fy{i}(3)-dT_fy{i}(end-1))*(dT_fy{i}(3:end-1)-dT_fy{i}(3))+dS_fy{i}(1),T_fy{i}(3:end-1),dS_fy{i}(1:end-2),'pchip');
    rhoi_fy{i} = (917-1.403*10^-1*T_sal_fy{i}); % pure ice density, Pounder (1965) for T_insitu
    F1_fy{i} = -4.732-22.45*T_sal_fy{i} - 0.6397*T_sal_fy{i}.^2 - 0.01074*T_sal_fy{i}.^3; % Cox and Weeks (1983)  
    F1_fy{i}(T_sal_fy{i}>-2) = -4.1221*10^-2 + -1.8407*10^1*T_sal_fy{i}(T_sal_fy{i}>-2).^1 + 5.8402*10^-1*T_sal_fy{i}(T_sal_fy{i}>-2).^2 + 2.1454*10^-1*T_sal_fy{i}(T_sal_fy{i}>-2).^3; % F1 from Lepparanta and Manninen (1988)
    F2_fy{i} = 8.903*10^-2 - 1.763*10^-2*T_sal_fy{i} - 5.33*10^-4*T_sal_fy{i}.^2 - 8.801*10^-6*T_sal_fy{i}.^3; % Cox and Weeks (1983)
    F2_fy{i}(T_sal_fy{i}>-2) = 9.0312*10^-2 + -1.6111*10^-2*T_sal_fy{i}(T_sal_fy{i}>-2).^1 + 1.2291*10^-4*T_sal_fy{i}(T_sal_fy{i}>-2).^2 + 1.3603*10^-4*T_sal_fy{i}(T_sal_fy{i}>-2).^3; % F2 from Lepparanta and Manninen (1988)
    vb_fy{i} = (1000).*(rhoi_fy{i}.*S_fy{i}(1:end-2)/1000)./(F1_fy{i}-rhoi_fy{i}.*S_fy{i}(1:end-2)/1000.*F2_fy{i}); % Brine volume for T_insitu, LM
    vb_bulk_fy(i) = mean(vb_fy{i}(vb_fy{i}>0 & vb_fy{i}<1000),'includenan');
end
for i = 22
    T_sal_fy{i} = interp1((dS_fy{i}(1)-dS_fy{i}(end-3))/(dT_fy{i}(3)-dT_fy{i}(end-1))*(dT_fy{i}(3:end-1)-dT_fy{i}(3))+dS_fy{i}(1),T_fy{i}(3:end-1),dS_fy{i}(1:end-3),'pchip');
    rhoi_fy{i} = (917-1.403*10^-1*T_sal_fy{i}); % pure ice density, Pounder (1965) for T_insitu
    F1_fy{i} = -4.732-22.45*T_sal_fy{i} - 0.6397*T_sal_fy{i}.^2 - 0.01074*T_sal_fy{i}.^3; % Cox and Weeks (1983)  
    F1_fy{i}(T_sal_fy{i}>-2) = -4.1221*10^-2 + -1.8407*10^1*T_sal_fy{i}(T_sal_fy{i}>-2).^1 + 5.8402*10^-1*T_sal_fy{i}(T_sal_fy{i}>-2).^2 + 2.1454*10^-1*T_sal_fy{i}(T_sal_fy{i}>-2).^3; % F1 from Lepparanta and Manninen (1988)
    F2_fy{i} = 8.903*10^-2 - 1.763*10^-2*T_sal_fy{i} - 5.33*10^-4*T_sal_fy{i}.^2 - 8.801*10^-6*T_sal_fy{i}.^3; % Cox and Weeks (1983)
    F2_fy{i}(T_sal_fy{i}>-2) = 9.0312*10^-2 + -1.6111*10^-2*T_sal_fy{i}(T_sal_fy{i}>-2).^1 + 1.2291*10^-4*T_sal_fy{i}(T_sal_fy{i}>-2).^2 + 1.3603*10^-4*T_sal_fy{i}(T_sal_fy{i}>-2).^3; % F2 from Lepparanta and Manninen (1988)
    vb_fy{i} = (1000).*(rhoi_fy{i}.*S_fy{i}(1:end-3)/1000)./(F1_fy{i}-rhoi_fy{i}.*S_fy{i}(1:end-3)/1000.*F2_fy{i}); % Brine volume for T_insitu, LM
    vb_bulk_fy(i) = mean(vb_fy{i}(vb_fy{i}>0 & vb_fy{i}<1000),'includenan');
end
clearvars F1_fy F2_fy F1_rho_fy F2_rho_fy F3_fy F1_pr_fy F2_pr_fy F3_pr_fy T_lab_fy F1_pr_rho_fy F2_pr_rho_fy F3_rho_fy rhoi_fy
for i = 1:length(S_sy) % SYI density calculations
    T_sy{i} = min(-0.1,T_sy{i});
    T_lab_sy{i} = ones(size(Srho_sy{i})) * hd_sy(i,8); % Lab temperature
    F1_pr_sy{i} = -4.732-22.45*T_lab_sy{i} - 0.6397*T_lab_sy{i}.^2 - 0.01074*T_lab_sy{i}.^3;
    F2_pr_sy{i} = 8.903*10^-2 - 1.763*10^-2*T_lab_sy{i} - 5.33*10^-4*T_lab_sy{i}.^2 - 8.801*10^-6*T_lab_sy{i}.^3;
    vb_pr_sy{i} = rho_sy{i} .* Srho_sy{i} ./ F1_pr_sy{i}; % brine volume for T_lab
    rhoi_pr_sy{i} = (917-1.403*10^-1*T_lab_sy{i}); % pure ice density, Pounder (1965)
    vg_pr_sy{i} = max((1-rho_sy{i}.*(F1_pr_sy{i}-rhoi_pr_sy{i}.*Srho_sy{i}/1000.*F2_pr_sy{i})./(rhoi_pr_sy{i}.*F1_pr_sy{i}))*1000,0,'includenan'); % gas volume for T_lab
    vg_pr_bulk_sy(i) = mean(vg_pr_sy{i},'omitnan');
    T_sal_sy{i} = interp1(dT_sy{i}(3:end),T_sy{i}(3:end),dS_sy{i},'pchip'); % in-situ temperature interpolation for salinity depth
    T_rho_sy{i} = interp1(dT_sy{i}(3:end),T_sy{i}(3:end),drho_sy{i},'pchip'); % in-situ temperature interpolation for density depth
    T_rho_sy{i} = min(-0.1,T_rho_sy{i});
    F3_pr_sy{i} = rhoi_pr_sy{i}.*Srho_sy{i}/1000./(F1_pr_sy{i}-rhoi_pr_sy{i}.*Srho_sy{i}/1000.*F2_pr_sy{i});
    rhoi_sy{i} = (917-1.403*10^-1*T_sal_sy{i}); % pure ice density, Pounder (1965) for T_insitu
    rhoi_rho_sy{i} = (917-1.403*10^-1*T_rho_sy{i}); % pure ice density, Pounder (1965) for T_insitu
    F1_sy{i} = -4.732-22.45*T_sal_sy{i} - 0.6397*T_sal_sy{i}.^2 - 0.01074*T_sal_sy{i}.^3; % Cox and Weeks (1983)  
    F1_sy{i}(T_sal_sy{i}>-2) = -4.1221*10^-2 + -1.8407*10^1*T_sal_sy{i}(T_sal_sy{i}>-2).^1 + 5.8402*10^-1*T_sal_sy{i}(T_sal_sy{i}>-2).^2 + 2.1454*10^-1*T_sal_sy{i}(T_sal_sy{i}>-2).^3; % F1 from Lepparanta and Manninen (1988)
    F2_sy{i} = 8.903*10^-2 - 1.763*10^-2*T_sal_sy{i} - 5.33*10^-4*T_sal_sy{i}.^2 - 8.801*10^-6*T_sal_sy{i}.^3; % Cox and Weeks (1983)
    F2_sy{i}(T_sal_sy{i}>-2) = 9.0312*10^-2 + -1.6111*10^-2*T_sal_sy{i}(T_sal_sy{i}>-2).^1 + 1.2291*10^-4*T_sal_sy{i}(T_sal_sy{i}>-2).^2 + 1.3603*10^-4*T_sal_sy{i}(T_sal_sy{i}>-2).^3; % F2 from Lepparanta and Manninen (1988)
    F1_rho_sy{i} = -4.732-22.45*T_rho_sy{i} - 0.6397*T_rho_sy{i}.^2 - 0.01074*T_rho_sy{i}.^3; % Cox and Weeks (1983)  
    F1_rho_sy{i}(T_rho_sy{i}>-2) = -4.1221*10^-2 + -1.8407*10^1*T_rho_sy{i}(T_rho_sy{i}>-2).^1 + 5.8402*10^-1*T_rho_sy{i}(T_rho_sy{i}>-2).^2 + 2.1454*10^-1*T_rho_sy{i}(T_rho_sy{i}>-2).^3; % F1 from Lepparanta and Manninen (1988)
    F2_rho_sy{i} = 8.903*10^-2 - 1.763*10^-2*T_rho_sy{i} - 5.33*10^-4*T_rho_sy{i}.^2 - 8.801*10^-6*T_rho_sy{i}.^3; % Cox and Weeks (1983)
    F2_rho_sy{i}(T_rho_sy{i}>-2) = 9.0312*10^-2 + -1.6111*10^-2*T_rho_sy{i}(T_rho_sy{i}>-2).^1 + 1.2291*10^-4*T_rho_sy{i}(T_rho_sy{i}>-2).^2 + 1.3603*10^-4*T_rho_sy{i}(T_rho_sy{i}>-2).^3; % F2 from Lepparanta and Manninen (1988)
    F3_sy{i} = rhoi_rho_sy{i}.*Srho_sy{i}/1000./(F1_rho_sy{i}-rhoi_rho_sy{i}.*Srho_sy{i}/1000.*F2_rho_sy{i});
    vb_sy{i} = (1000).*(rhoi_sy{i}.*S_sy{i}/1000)./(F1_sy{i}-rhoi_sy{i}.*S_sy{i}/1000.*F2_sy{i}); % Brine volume for T_insitu, LM
    vb_rho_sy{i} = vb_pr_sy{i} .* F1_pr_sy{i} ./ F1_rho_sy{i}; % Brine volume for T_insitu, CW + LM
    vb_bulk_sy(i) = mean(vb_sy{i} (vb_sy{i}>0 & vb_sy{i}<1000),'includenan');
    rho_si_pr_sy{i} = (-vg_pr_sy{i}/1000+1).*rhoi_rho_sy{i}.*F1_rho_sy{i}./(F1_rho_sy{i}-rhoi_rho_sy{i}.*Srho_sy{i}/1000.*F2_rho_sy{i});
    rho_si_pr_bulk_sy(i) = mean(rho_si_pr_sy{i}(vb_rho_sy{i}>0 & vb_rho_sy{i}<1000),'omitnan');
end
Sb_fy_ct = [39.5 58.5 9.73 6.87 7.34 3.69 11.61]; % in-situ brine salinity from CT bottom sample, FYI
clearvars F1_sy F2_sy F1_rho_sy F2_rho_sy F3_sy F1_pr_sy F2_pr_sy F3_pr_sy T_lab_sy rhoi_sy

% sea-ice density from ALS and ROV
t_ALS = datetime(['30-Jun-2020';'04-Jul-2020';'17-Jul-2020';'22-Jul-2020']); fb_sn_FYI_ALS = [0.206 0.111 0.177 0.1768]; % Values are calculated in Figure 5
sn_FYI_ALS = [0.08 0.048 0.048 0.048]; fb_FYI_ALS = fb_sn_FYI_ALS - sn_FYI_ALS;
t_ROV = datetime(['24-Jun-2020';'28-Jun-2020';'01-Jul-2020';'07-Jul-2020';'14-Jul-2020';'21-Jul-2020';'28-Jul-2020']);
d_FYI_ROV = [-1.4298 -1.4023 -1.3897 -1.3851 -1.0839 -1.019 -0.869]; % Average draft from ROV sonar for FYI ROV site; values are calculated in Figure 5
t_int = convertTo(t_ALS(1),"datenum"):1/24:convertTo(t_ALS(4),"datenum"); td_int = datetime(t_int,'ConvertFrom','datenum');
fb_sn_ALS_int = interp1(convertTo(t_ALS([1 3:4]),"datenum"),fb_sn_FYI_ALS([1 3:4]),t_int,'linear');
fb_ALS_int = interp1(convertTo(t_ALS([1 3:4]),"datenum"),fb_FYI_ALS([1 3:4]),t_int,'linear');
d_ROV_int = interp1(convertTo(t_ROV,"datenum"),d_FYI_ROV,t_int,'linear');
rhow = 1017; rhosn = 400;
rho_si_rov_als = (rhow*-d_ROV_int - rhosn*(fb_sn_ALS_int-fb_ALS_int)) ./ (-d_ROV_int+fb_ALS_int); % density from ALS fb, ROV draft, snow depth from Webster

rhosw_imb = [1017.5 1015.0 1015.5];
h_mp_als = [0.0754 0.0901 0.0996]; % melt pond depth for ALS-ROV area from Fuchs
a_mp_rov_fuchs = [0.2095 0.2256 0.2437]; % melt pond fraction for ALS-ROV area from Fuchs
fb_filled = [0.2019 0.1744 0.1710]; % ALS-ROV area freeboard from ALS with filled gaps; values are calculated in Figure 5
rho_si_rov_als_mp(1) = (rhosw_imb(1)*-d_ROV_int(1) - rhosn*sn_FYI_ALS(1)*(1-a_mp_rov_fuchs(1)) - 1000*h_mp_als(1)*a_mp_rov_fuchs(1)) / ( (-d_ROV_int(1)+fb_FYI_ALS(1))*(1-a_mp_rov_fuchs(1)) + (fb_filled(1)-h_mp_als(1)-d_ROV_int(1))*a_mp_rov_fuchs(1) ); % 30 June
rho_si_rov_als_mp(2) = (rhosw_imb(2)*-d_ROV_int(409) - rhosn*sn_FYI_ALS(3)*(1-a_mp_rov_fuchs(2)) - 1000*h_mp_als(2)*a_mp_rov_fuchs(2)) / ( (-d_ROV_int(409)+fb_FYI_ALS(3))*(1-a_mp_rov_fuchs(2)) + (fb_filled(2)-h_mp_als(2)-d_ROV_int(409))*a_mp_rov_fuchs(2) ); % 17 July
rho_si_rov_als_mp(3) = (rhosw_imb(3)*-d_ROV_int(end) - rhosn*sn_FYI_ALS(4)*(1-a_mp_rov_fuchs(3)) - 1000*h_mp_als(3)*a_mp_rov_fuchs(3)) / ( (-d_ROV_int(end)+fb_FYI_ALS(4))*(1-a_mp_rov_fuchs(3)) + (fb_filled(3)-h_mp_als(3)-d_ROV_int(end))*a_mp_rov_fuchs(3) ); % 22 July

% transect snow, SSL and melt pond thickness
t_tr_old = datetime(['30-Oct-2019';'06-Nov-2019';'13-Nov-2019';'04-Dec-2019';'02-Jan-2020';'08-Jan-2020';'29-Jan-2020';'19-Feb-2020';'25-Feb-2020';'04-Mar-2020';'21-Mar-2020';'08-Apr-2020';'23-Apr-2020';'10-May-2020']);
h_sn_tr_old =       [       0.094          0.10         0.10          0.107         0.111         0.104          0.12         0.144         0.134         0.121         0.124         0.124         0.147         0.157 ]; % snow or SSL depth from Itkin (2022), level ice
t_tr = datetime(['19-Jun-2020';'21-Jun-2020';'27-Jun-2020';'29-Jun-2020';'30-Jun-2020';'02-Jul-2020';'03-Jul-2020';'04-Jul-2020';'05-Jul-2020';'06-Jul-2020';'07-Jul-2020';'08-Jul-2020';'10-Jul-2020';'13-Jul-2020';'14-Jul-2020';'17-Jul-2020';'19-Jul-2020';'20-Jul-2020';'21-Jul-2020';'25-Jul-2020';'26-Jul-2020';'27-Jul-2020']);
h_sn_tr =       [       0.129          0.135       0.0942         0.089          0.10        0.1048        0.0863         0.068         0.077        0.0584         0.052        0.0584        0.0531         0.053         0.054         0.064        0.0531        0.0518         0.061         0.045          0.041         0.044]; % snow or SSL depth from Webster (2022)
clearvars i j vb_limit T_T66 T_T62 vb_sy vb_pr_sy dT_fy dT_sy fbrho_fy fbrho_sy hd_fy hd_sy

%% Figure 1: density historical overview
load("density_data.mat","rho_fy_mosaic","t_mosaic","T_hist","rho_hist","T_mosaic","src_hist","t_hist","c"); % data import
% analytical gas-free density for FYI with a fixed salinity
T_anal = -30:0.05:-0.1; Si = 2; % temperature range and post-melt salinity
F1 = -4.732-22.45*T_anal - 0.6397*T_anal.^2 - 0.01074*T_anal.^3; % Cox and Weeks (1983)
F2 = 8.903*10^-2 - 1.763*10^-2*T_anal - 5.33*10^-4*T_anal.^2 - 8.801*10^-6*T_anal.^3; % Cox and Weeks (1983)
f1t0 = -4.1221*10^-2; f1t1 = -1.8407*10^1; f1t2 = 5.8402*10^-1; f1t3 = 2.1454*10^-1; % F1 coefficients from Lepparanta and Manninen (1988)
f2t0 = 9.0312*10^-2; f2t1 = -1.6111*10^-2; f2t2 = 1.2291*10^-4; f2t3 = 1.3603*10^-4; % F2 coefficients from Lepparanta and Manninen (1988)
F1_lep = f1t0 + f1t1*T_anal.^1 + f1t2*T_anal.^2 + f1t3*T_anal.^3; % F1 from Lepparanta and Manninen (1988)
F2_lep = f2t0 + f2t1*T_anal.^1 + f2t2*T_anal.^2 + f2t3*T_anal.^3; % F2 from Lepparanta and Manninen (1988)
F1(T_anal>-2) = F1_lep(T_anal>-2); F2(T_anal>-2) = F2_lep(T_anal>-2);
rho_i = 916.8 - 0.1403*T_anal; % pure ice density, Pounder (1965)
rho_si_cox = 1./((F1 - rho_i.*Si/1000.*F2)./(rho_i.*F1)); % Cox and Weeks (1983)
clearvars Si F1 f1t0 f1t1 f1t2 f1t3 f2t0 f2t1 f2t2 f2t3 F1_lep F2_lep F2 rho_i
t_fons = datetime(['01-Jan-2020';'01-Mar-2020';'01-Mar-2020';'01-Jun-2020';'01-Jun-2020';'01-Sep-2020';'01-Sep-2020';'01-Dec-2020';'01-Dec-2020';'31-Dec-2020']);
rho_fons = [920 920 915 915 875 875 900 900 920 920]; % from Fons et al. (2023), doi:10.5194/tc-17-2487-2023

figure
tile = tiledlayout(1,2); tile.TileSpacing = 'compact'; tile.Padding = 'compact';
nexttile
p = plot(t_hist([4 6:32 34:end]),rho_hist([4 6:32 34:end]),'kx'); p.MarkerSize = 4; hold on % Timco & Frederking 1996
p = plot(t_hist(38),rho_hist(38),'x','Color',c{3}); p.MarkerSize = 4; % Alexandrov 2010
p = plot(t_hist(34:37),rho_hist(34:37),'x','Color',c{5}); p.MarkerSize = 4; % Forsstrom 2011
p = plot(t_hist(33),rho_hist(33),'x','Color',c{2}); p.MarkerSize = 4; % Wang 2020
p = plot(t_hist(42:43),rho_hist(42:43),'x','Color',c{4}); p.MarkerSize = 4; % Jutila 2022
p = plot(t_fons,rho_fons,':','Color','k','LineWidth',2.5); p.MarkerSize = 4; % Fons et al. (2023), doi:10.5194/tc-17-2487-2023
p = plot(t_mosaic(8:21),rho_fy_mosaic(8:21),'o:','Color',c{1},'LineWidth',3.0); p.MarkerSize = 2; set(p,'markerfacecolor',get(p,'color')); % MOSAiC, winter and spring
p = plot(t_mosaic(1:7),rho_fy_mosaic(1:7),'o:','Color',c{1},'LineWidth',3.0); p.MarkerSize = 2; set(p,'markerfacecolor',get(p,'color')); % MOSAiC, autumn
p = plot(t_mosaic(22:27),rho_fy_mosaic(22:27),'o:','Color',c{1},'LineWidth',3.0); p.MarkerSize = 2; set(p,'markerfacecolor',get(p,'color')); % MOSAiC, leg 5
leg = legend('Timco et al., 1996','Alexandrov, 2010','Forsstrom et al., 2011','Wang et al., 2020','Jutila et al., 2022','Fons et al., 2023','MOSAiC','box','off','NumColumns',1); set(leg,'FontSize',6,'Location','southwest'); leg.ItemTokenSize = [30*0.5,18*0.5];
hYLabel = ylabel('FYI density (kg  m^-^3)'); set([hYLabel gca],'FontSize',8,'FontWeight','normal');
ax = gca; ax.XTick = datetime(['01-Jan-2020';'01-Mar-2020';'01-May-2020';'01-Jul-2020';'01-Sep-2020';'01-Nov-2020';'01-Jan-2021']); datetick('x','mmm','keepticks'); xtickangle(0); % time

nexttile
plot(T_anal,rho_si_cox,'--','Color','k','LineWidth',1); hold on % Analytical gas-free FYI density following Cox and Weeks (1983)
p = plot(T_hist([1:32 34:end]),rho_hist([1:32 34:end]),'kx'); p.MarkerSize = 4; % Timco and Frederking, 1996
p = plot(T_hist(33),rho_hist(33),'x','Color',c{2}); p.MarkerSize = 4; % Wang 2022
p = plot(T_mosaic,rho_fy_mosaic,'o','Color',c{1},'LineWidth',2); p.MarkerSize = 2; set(p,'markerfacecolor',get(p,'color')); % MOSAiC
p = plot(T_mosaic(22:27),rho_fy_mosaic(22:27),'o','Color',c{1},'LineWidth',2); p.MarkerSize = 2; set(p,'markerfacecolor',get(p,'color')); % MOSAiC
leg = legend('Gas-free, Cox and Weeks, 1983','Timco et al., 1996','Wang et al., 2020','MOSAiC','','box','off','NumColumns',1); set(leg,'FontSize',6,'Location','southwest'); leg.ItemTokenSize = [30*0.7,18*0.7];
hXLabel = xlabel('FYI temperature (°C)'); hYLabel = ylabel('FYI density (kg  m^-^3)'); set([hXLabel hYLabel gca],'FontSize',8,'FontWeight','normal');
ylim([780 940]);
annotation('textbox',[0 .51 0.02 .51],'String','(a)','EdgeColor','none','HorizontalAlignment','center','FontSize',9);
annotation('textbox',[0.33 .51 0.33 .51],'String','(b)','EdgeColor','none','HorizontalAlignment','center','FontSize',9);
clearvars ax hXLabel hYLabel i leg p T_anal rho_si_cox tile

%% Figure 2: preparation, density contourplot for FYI
vg_lvl = rho_si_pr_fy'; vg_lvl(7)=[]; vg_lvl(4)=[]; td=datenum(t_fy)';
drho_fy_lvl = drho_fy; drho_fy_lvl(7)=[]; drho_fy_lvl(4)=[]; 
td_lvl = td_rho_fy; td_lvl(7)=[]; td_lvl(4)=[];
td_lvl{16}=td_lvl{16}./td_lvl{16}(1)*datenum(datetime('10-Jun-2020')); td(18)=datenum(datetime('10-Jun-2020')); % fake date for May 4 coring
drho_fy_lvl{19}(31:32)=[]; drho_fy_lvl{20}(30:31)=[]; vg_lvl{19}(31:32)=[]; vg_lvl{20}(30:31)=[]; td_lvl{19}(31:32)=[]; td_lvl{20}(30:31)=[]; % remove false bottom
drho_fy_lvl{21}(18:20)=[]; vg_lvl{21}(18:20)=[]; td_lvl{21}(18:20)=[]; % remove 
drho_fy_lvl{1}(7:8)=[]; vg_lvl{1}(7:8)=[]; td_lvl{1}(7:8)=[]; % remove NaN week 1
drho_fy_lvl{2}(8:9)=[]; vg_lvl{2}(8:9)=[]; td_lvl{2}(8:9)=[]; % remove NaN week 2
drho_fy_lvl{3}(4:5)=[]; vg_lvl{3}(4:5)=[]; td_lvl{3}(4:5)=[]; % remove NaN week 3
drho_fy_lvl{4}(2)=[]; vg_lvl{4}(2)=[]; td_lvl{4}(2)=[]; % remove NaN week 4
vg_vect = cell2mat(vg_lvl);
drho_fy_vect = cell2mat(drho_fy_lvl);
tdvg_fy_vect = cell2mat(td_lvl);
x = tdvg_fy_vect; y = drho_fy_vect; z = vg_vect;
xv = linspace(min(x), max(x), 275*2*2);
yv = linspace(min(y), max(y), 200*2*2);
[X,Y] = meshgrid(xv, yv);
Z = griddata(x,y,z,X,Y);
for i = 1:21; drho_fy_bot_lvl(i) = drho_fy_lvl{i}(end); end % depth of the bottom density sample
drho_bot_int = interp1(td([1:3 5:6 8:23]),drho_fy_bot_lvl,xv,'linear'); % bottom interpolation
for i = 1:length(xv); [~,idx(i)] = min(abs(yv - drho_bot_int(i))); end
for i = 1:length(xv) % removing interpolated values below ice core bottom
    for j = 1:length(yv)
        if j < idx(i)
            Z(j,i) = NaN;
        end
    end
end
for i = 1:21; drho_fy_top_lvl(i) = drho_fy_lvl{i}(1); end % depth of the top surface density sample
drho_top_int = interp1(td([1:3 5:6 8:23]),drho_fy_top_lvl,xv,'linear'); % top surface interpolation
for i = 1:length(xv); [~,idx(i)] = min(abs(yv - drho_top_int(i))); end
for i = 1:length(xv) % removing interpolated values above top surface of ice core
    for j = 1:length(yv)
        if j > idx(i)
            Z(j,i) = NaN;
        end
    end
end
[dem,~] = geotiffread('20200722_01_PS122-4_48-69_orthomosaic_hr_l2.tif'); % from Neckel et al., doi:10.1594/PANGAEA.949433
x = 0:0.5:(size(dem,1)/2-0.5); x = x-4895; x = x/1000;
y = 0:0.5:(size(dem,2)/2-0.5); y = -(y-4972); y = y/1000;
J = imrotate(dem(:,:,1:3),-124.2,'bilinear','crop');

%% Figure 2: map, orthomosaic and density contourplot
figure
tile = tiledlayout(1,3); tile.TileSpacing = 'compact'; tile.Padding = 'compact'; ax = gobjects(1,3);
ax(1) = nexttile; % bathymetry with ice floe drift trajectory
m_proj('lambert','lons',[-16 24],'lat',[76 84]);
[cs,~]=m_etopo2('contourf',-6000:100:0,'edgecolor','none');
m_gshhs_h('patch',[.7 .7 .7],'edgecolor','none'); % coastline
m_grid('linewi',1,'layer','top','FontSize',8); % axis settings
edge0728_85(1,2) = edge0728_85(1,2)-0.1; m_line(smooth(edge0728_85(:,1),5,'sgolay'),smooth(edge0728_85(:,2),5,'sgolay'),'color',[0.5 0.5 0.5],'linewi',0.05); % sea ice edge AMSR2, 85%
m_line(smooth(edge0728_15(:,1),5,'sgolay'),smooth(edge0728_15(:,2),5,'sgolay'),'color','k','linewi',0.05); % sea ice edge AMSR2, 15%
clim([-6000 0]);
colormap(ax(1),m_colmap('blue'));
C = linspecer(length(7500:13165));
m_line(long_gps_T66(9037:13165),lat_gps_T66(9037:13165),'color',c{2},'linewi',2.0); hold on % drift during sonar surveys
for i = 9038:10:13165
    m_line(long_gps_T66(i),lat_gps_T66(i),'marker','o','color',C(i-7499,:),'markerfacecolor',C(i-7499,:),'linewi',0.1,'markersize',2);
end
m_line(long_gps_T66(9037),lat_gps_T66(9037),'marker','o','color','k','markerfacecolor',C(9038-7499,:),'linewi',0.5,'markersize',5); % 22 June
m_line(long_gps_T66(11389),lat_gps_T66(11389),'marker','o','color','k','markerfacecolor',C(11389-7499,:),'linewi',0.5,'markersize',5); % 22 June
m_line(long_gps_T66(12061),lat_gps_T66(12061),'marker','o','color','k','markerfacecolor',C(12061-7499,:),'linewi',0.5,'markersize',5); % 06 July
m_line(long_gps_T66(12397),lat_gps_T66(12397),'marker','o','color','k','markerfacecolor',C(12397-7499,:),'linewi',0.5,'markersize',5); % 13 July
m_line(long_gps_T66(12733),lat_gps_T66(12733),'marker','o','color','k','markerfacecolor',C(12733-7499,:),'linewi',0.5,'markersize',5); % 20 July
m_line(long_gps_T66(13165),lat_gps_T66(13165),'marker','o','color','k','markerfacecolor',C(13165-7499,:),'linewi',0.5,'markersize',5); % 29 July
p = m_text(16,79,'Svalbard'); set(p,'HorizontalAlignment','center','FontSize',7);
p = m_text(18.27,83.90-0.4,'4 May'); set(p,'HorizontalAlignment','center','FontSize',7);
p = m_text(9.96+6,82.04+0.2,'22 June'); set(p,'HorizontalAlignment','center','FontSize',7);
p = m_text(long_gps_T66(12061),lat_gps_T66(12061)+0.4,'6 July'); set(p,'HorizontalAlignment','right','FontSize',7);
p = m_text(long_gps_T66(12397)-2,lat_gps_T66(12397),'13 July'); set(p,'HorizontalAlignment','right','FontSize',7);
p = m_text(long_gps_T66(12733)-2,lat_gps_T66(12733),'20 July'); set(p,'HorizontalAlignment','right','FontSize',7);
p = m_text(-2.71-2,79.36,'29 July'); set(p,'HorizontalAlignment','right','FontSize',7);
p = m_text(-7,78.3,'SIC 80%','color',[0.5 0.5 0.5]); set(p,'HorizontalAlignment','center','FontSize',6); set(p,'Rotation',45);
p = m_text(-8,76.9,'SIC 15%'); set(p,'HorizontalAlignment','center','FontSize',6); set(p,'Rotation',45);
title('Coring locations','FontSize',7,'FontWeight','normal');

nexttile % optical image of CO2 ice floe with coring and ROV sites
imagesc(x,y,J); hold on
values = spcrv([ort_fyi(:,1)'/1000;ort_fyi(:,2)'/1000],3); plot(values(1,:),values(2,:),'-','color',c{1},'LineWidth',3); hold on % FYI area
values = spcrv([ort_fyi(1:80,1)'/1000+0.018;ort_fyi(1:80,2)'/1000],3); plot(values(1,:),values(2,:),'-','color',c{2},'LineWidth',3); % FYI-SYI border
values = spcrv([ort_syi(:,1)'/1000;ort_syi(:,2)'/1000],3); plot(values(1,:),values(2,:),'-','color',c{2},'LineWidth',3); % SYI area
values = spcrv([ort_rov_fyi(:,1)'/1000;ort_rov_fyi(:,2)'/1000],3); plot(values(1,:),values(2,:),'-','color',c{4},'LineWidth',3); % ROV FYI area
p = plot(-0.35,0.24,'o','color','k'); p.MarkerSize = 4; set(p,'markerfacecolor',c{1});
p = plot(+0.38,0.54,'o','color','k'); p.MarkerSize = 4; set(p,'markerfacecolor',c{2});
leg = legend('FYI CO2','','SYI CO2','FYI ROV','FYI coring','SYI coring','box','on'); set(leg,'FontSize',6,'Location','northwest'); leg.ItemTokenSize = [30*0.2,18*0.2];
leg.BoxFace.ColorType='truecoloralpha'; leg.BoxFace.ColorData=uint8(255*[1 1 1 0.7]'); % transparent legend
hXLabel = xlabel('x (km)'); hYLabel = ylabel('y (km)'); set([hXLabel hYLabel gca],'FontSize',8,'FontWeight','normal');
xticks(-0.8:0.2:0.6); xlim([-0.8 0.6]); xtickangle(0); ylim([-0.1 1.01]); set(gca,'YDir','normal');
title('Aerial image of Central Observatory 2, 22 July','FontSize',7,'FontWeight','normal');

ax(3) = nexttile; % density contourplot
range = [700 800 880 900 920 950 1000]; % Density graph accuracy
[~,~] = contourf(X,Y,Z,range,'-','ShowText','off','LabelSpacing',800,'LineColor','w'); hold on % clabel(C,h,'FontSize',6,'Color','k');
p = plot(td([1:3 5:6 8:23]),drho_fy_bot_lvl,'ko-'); p.MarkerSize = 1.5; set(p,'markerfacecolor',get(p,'color')); % bottom density core depth
p = plot(td([1:3 5:6 8:23]),drho_fy_top_lvl,'ko-'); p.MarkerSize = 1.5; set(p,'markerfacecolor',get(p,'color')); % top density core depth
p = text(td(19)-10,-0.06,'880'); set(p,'Color','k','HorizontalAlignment','right','FontSize',6);
p = text(td(19)-10,-0.12,'900'); set(p,'Color','k','HorizontalAlignment','right','FontSize',6);
p = text(td(19)-10,-1.41,'920'); set(p,'Color','k','HorizontalAlignment','right','FontSize',6);
p = text(td(20)-3,-0.16,'800'); set(p,'Color','k','HorizontalAlignment','right','FontSize',6);
p = text(td(20)-3,-0.28,'880'); set(p,'Color','k','HorizontalAlignment','right','FontSize',6);
p = text(td(20)-3,-0.34,'900'); set(p,'Color','k','HorizontalAlignment','right','FontSize',6);
p = text(td(20)-3,-0.43,'920'); set(p,'Color','k','HorizontalAlignment','right','FontSize',6);
p = text(td(20)-3,-0.57,'920'); set(p,'Color','k','HorizontalAlignment','right','FontSize',6);
p = text(td(20)-3,-0.74,'920'); set(p,'Color','k','HorizontalAlignment','right','FontSize',6);
p = text(td(20)-3,-1.01,'920'); set(p,'Color','k','HorizontalAlignment','right','FontSize',6);
p = text(td(20)-3,-1.38,'900'); set(p,'Color','k','HorizontalAlignment','right','FontSize',6);
p = text(td(22)-5,-0.11,'900'); set(p,'Color','k','HorizontalAlignment','right','FontSize',6);
p = text(td(22)-5,-1.17,'880'); set(p,'Color','k','HorizontalAlignment','right','FontSize',6);
p = text(td(22)-3,-0.55,'900'); set(p,'Color','k','HorizontalAlignment','right','FontSize',6);
p = text(td(23)-5,+0.05,'900'); set(p,'Color','k','HorizontalAlignment','right','FontSize',6);
p = text(td(23)-5,-0.26,'880'); set(p,'Color','k','HorizontalAlignment','right','FontSize',6);
p = text(td(23)-5,-0.48,'880'); set(p,'Color','k','HorizontalAlignment','right','FontSize',6);
p = text(td(23)-5,-0.66,'880'); set(p,'Color','k','HorizontalAlignment','right','FontSize',6);
p = text(td(23)-5,-0.87,'880'); set(p,'Color','k','HorizontalAlignment','right','FontSize',6);
p = text(td(23)-1,-0.10,'880'); set(p,'Color','k','HorizontalAlignment','right','FontSize',6);
p = text(td(23)-1,+0.07,'880'); set(p,'Color','k','HorizontalAlignment','right','FontSize',6);
p = text(td(23)-1,-0.40,'880'); set(p,'Color','k','HorizontalAlignment','right','FontSize',6);
hYLabel = ylabel('Depth (m)'); set([hYLabel gca],'FontSize',8,'FontWeight','normal');
colormap(ax(3),flipud(summer));
clim([749 951]);
yticks(-2.2:0.2:0.4); ylim([-1.7 0.35]); % depth limits
t_start = datenum(datetime('09-Jun-2020 00:00')); t_end = datenum(datetime('30-Jul-2020 00:00')); xlim([t_start t_end]);
ax = gca; ax.XTick = datenum(datetime(['08-Jun-2020';'20-Jun-2020';'30-Jun-2020';'10-Jul-2020';'20-Jul-2020';'30-Jul-2020'])); datetick('x','dd mmm','keepticks'); xtickangle(0); % time
xticklabels({'1 May','20 Jun','30 Jun','10 Jul','20 Jul','30 Jul'});
hBar1 = colorbar('West'); set(hBar1,'Position',[0.910 0.17 0.010 0.14]); zticks(750:50:950);
p = text(td(22)-0,-1.51,'Ice density'); set(p,'Color','k','HorizontalAlignment','right','FontSize',7);
p = text(td(22)-0,-1.60,'(kg m^-^3)'); set(p,'Color','k','HorizontalAlignment','right','FontSize',7);
p = text(td(18)+0,0.28,'Bulk ice density'); set(p,'Color',c{1},'HorizontalAlignment','left','FontSize',7);
p = fill([td([1:3 5:6 8:23])' fliplr(td([1:3 5:6 8:23])')],[drho_fy_top_lvl fliplr(drho_fy_top_lvl+hs_fy([1:3 5:6 8:23]))],1,'FaceColor',[0.7 0.7 0.7],'edgecolor','none'); set(p,'facealpha',.3); % shaded snow
p = plot([td(18)-10 td(23)],[0 0],'k:'); set(p,'linewidth',2.5); % waterline
p = text(td(20)-4,-0.05,'Waterline'); set(p,'Color','k','HorizontalAlignment','center','FontSize',7);
fill([td(18)+4 td(18)+8 td(18)+8 td(18)+4],[-2 -2 0.2 0.2],1,'FaceColor','w','edgecolor','none'); % shaded observation gap
title('First-year ice density','FontSize',7,'FontWeight','normal');
for i = 19:22; p = text(td(i),drho_fy_top_lvl(i-2)+0.08,sprintf('%.0f',rho_si_pr_bulk_fy(i) )); set(p,'Color',c{1},'HorizontalAlignment','center','FontSize',7); end % bulk density, text
for i = 23; p = text(td(i)-1.8,drho_fy_top_lvl(i-2)+0.08,sprintf('%.0f',rho_si_pr_bulk_fy(i) )); set(p,'Color',c{1},'HorizontalAlignment','center','FontSize',7); end % bulk density, text
for i = 18; p = text(td(i)+1.6,drho_fy_top_lvl(i-2)+0.08,sprintf('%.0f',rho_si_pr_bulk_fy(i) )); set(p,'Color',c{1},'HorizontalAlignment','center','FontSize',7); end % bulk density, text

axes('Position',[.7335 .133 .0182 .03]); % axis break (c)
px=[0 0]; py=[0 1]; w=2; px2=px+w; fill([px flip(px2)],[py flip(py)],'w','EdgeColor','none'); hold all;
plot(px,py,'k',px2,py,'k','linewidth',0.1); box off; axis off;

annotation('textbox',[0 .51 0.02 .51],'String','(a)','FontSize',9,'EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.225 .51 0.225 .51],'String','(b)','FontSize',9,'EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.44 .51 0.44 .51],'String','(c)','FontSize',9,'EdgeColor','none','HorizontalAlignment','center');

%% Figure 3: Illustration of different freeboard and draft for pre-melt sea ice, and for sea ice with undrained and drained melt ponds
amp = 0.20; hmp = 0.25; hi = 1.0; rhosi = 900; rhow = 1000; rhosn = 300; hsn = 0.15; % chosen physical properties for illustration
fb_in = hi - hi*rhosi/rhow; % freeboard before ponding
fb_ponded = hi*(1-(1-amp)*rhosi/rhow)+(hmp-hi)*amp*rhosi/rhow-amp*hmp; % freeboard, ponded undrained
fb_drained_unponded = hi-rhosi/rhow*(hi-amp*hmp); % drained, unponded part only
fb_drained_eff = fb_drained_unponded*(1-amp) + (fb_drained_unponded-hmp)*amp; % drained, both ponded and unponded parts
fb_deep = (rhow-rhosi)*(hmp*amp-hi)/((amp-1)*rhow); % freeboard, drained, deep, unponded part only
fb_deep_eff = fb_deep*(1-amp);

fb_in_sn = hi-(rhosn*hsn+rhosi*hi)/rhow; % freeboard, pre-melt, with snow
fb_ponded_sn = ((hi-hmp*amp)*(rhow-rhosi)-rhosn*(1-amp)*hsn)/rhow; % freeboard, ice with undrained ponds, with snow
fb_drained_sn = hi-(rhosi*(hi-hmp*amp)+rhosn*(hsn-hsn*amp))/rhow; % freeboard, unponded ice with drained ponds, with snow
fb_drained_sn_eff = fb_drained_sn*(1-amp)+(fb_drained_sn-hmp)*amp; % freeboard, both ponded and unponded ice with drained ponds, with snow
fb_deep_sn = (rhow-rhosi)*(hi-hmp*amp)/((1-amp)*rhow)-rhosn*hsn/rhow; % freeboard, drained, deep, unponded part only
fb_deep_sn_eff = fb_deep_sn*(1-amp); % freeboard, drained, deep, both unponded and unpoded ice
fb_drained_true = fb_deep; if fb_drained_unponded > hmp; fb_drained_true = fb_drained_unponded; end % condition for deep or shallow pond assumption
fb_drained_sn_true = fb_deep_sn; if fb_drained_sn > hmp; fb_drained_sn_true = fb_drained_sn; end % condition for deep or shallow pond assumption, with snow

figure
tile = tiledlayout(1,3); tile.TileSpacing = 'compact'; tile.Padding = 'none';
nexttile
w = 2.2; fb = fb_in_sn;
plot([-w/4 0],[0 0],'color',c{1}); hold on; % water
p = text(-w/4,-0.05,'sea'); set(p,'Color',c{1},'HorizontalAlignment','left','FontSize',11);
p = text(-w/4,-0.12,'level'); set(p,'Color',c{1},'HorizontalAlignment','left','FontSize',11);
plot([0 w w 0 0],[fb fb fb+hsn fb+hsn fb],'k'); % snow left
p = text(+0.05,fb+hsn-0.05,'snow'); set(p,'Color','k','HorizontalAlignment','left','FontSize',11);
plot([0 w w 0 0],[fb fb fb-hi fb-hi fb],'k'); % ice
p = text(+0.05,fb-0.05,'ice'); set(p,'Color','k','HorizontalAlignment','left','FontSize',11);
fill([0 w w 0 0],[fb fb fb-hi fb-hi fb],1,'FaceColor',c{1},'EdgeColor','none'); alpha(.1);
xlim([-w/4-0.05 w+0.05]); ylim([-1 0.35]); set(gca,'xcolor','w','ycolor','w','xtick',[],'ytick',[]);

nexttile
w = 2; fb = fb_ponded;
plot([-w/4 0],[0 0],'color',c{1}); hold on; % water
p = text(-w/4,-0.05,'sea'); set(p,'Color',c{1},'HorizontalAlignment','left','FontSize',11);
p = text(-w/4,-0.12,'level'); set(p,'Color',c{1},'HorizontalAlignment','left','FontSize',11);
plot([0 w w 0 0],[fb fb fb-hi fb-hi fb],'k'); % ice
fill([0 w w 0 0],[fb fb fb-hi fb-hi fb],1,'FaceColor',c{1},'EdgeColor','none'); alpha(.1); % ice fill
p = text(+0.05,fb-0.05,'ice'); set(p,'Color','k','HorizontalAlignment','left','FontSize',11); % ice text
p = text(+0.05,fb+0.15,'undrained'); set(p,'Color','k','HorizontalAlignment','left','FontSize',11); % text
p = text(+0.05,fb+0.07,'unponded'); set(p,'Color','k','HorizontalAlignment','left','FontSize',11); % text
p = text(w-0.05,fb+0.15,'undrained'); set(p,'Color','k','HorizontalAlignment','right','FontSize',11); % ice text
p = text(w-0.05,fb+0.07,'ponded'); set(p,'Color','k','HorizontalAlignment','right','FontSize',11); % text
plot([w/2 w w w/2 w/2],[fb fb fb-hmp fb-hmp fb],'k'); % pond
p = fill([w/2 w w w/2 w/2],[fb fb fb-hmp fb-hmp fb],1,'FaceColor',c{1},'EdgeColor','none'); alpha(p,.5); % pond fill
p = text(w-0.05,fb-0.05,'melt pond'); set(p,'Color','k','HorizontalAlignment','right','FontSize',11); % pond text
xlim([-w/4-0.05 w+0.05]); ylim([-1 0.35]); set(gca,'xcolor','w','ycolor','w','xtick',[],'ytick',[]);

nexttile
w = 2; fb = fb_drained_true;
plot([-w/4 0],[0 0],'color',c{1}); hold on; % water
p = text(-w/4,-0.05,'sea'); set(p,'Color',c{1},'HorizontalAlignment','left','FontSize',11);
p = text(-w/4,-0.12,'level'); set(p,'Color',c{1},'HorizontalAlignment','left','FontSize',11);
plot([0 w/2 w/2 w w 0 0],[fb fb fb-hmp fb-hmp fb-hi fb-hi fb],'k'); % ice
fill([0 w/2 w/2 w w 0 0],[fb fb fb-hmp fb-hmp fb-hi fb-hi fb],1,'FaceColor',c{1},'EdgeColor','none'); alpha(.1); % ice fill
p = text(+0.05,fb-0.05,'ice'); set(p,'Color','k','HorizontalAlignment','left','FontSize',11); % ice text
p = text(+0.05,fb+0.15,'drained'); set(p,'Color','k','HorizontalAlignment','left','FontSize',11); % text
p = text(+0.05,fb+0.07,'unponded'); set(p,'Color','k','HorizontalAlignment','left','FontSize',11); % text
p = text(w-0.05,fb+0.15,'drained'); set(p,'Color','k','HorizontalAlignment','right','FontSize',11); % ice text
p = text(w-0.05,fb+0.07,'ponded'); set(p,'Color','k','HorizontalAlignment','right','FontSize',11); % text
plot([w/2 w w w/2 w/2],[0 0 fb-hmp fb-hmp 0],'k'); % pond
p = fill([w/2 w w w/2 w/2],[0 0 fb-hmp fb-hmp 0],1,'FaceColor',c{1},'EdgeColor','none'); alpha(p,.5); % pond fill
p = text(w-0.05,0-0.04,'melt pond'); set(p,'Color','k','HorizontalAlignment','right','FontSize',11); % pond text
xlim([-w/4-0.05 w+0.05]); ylim([-1 0.35]); set(gca,'xcolor','w','ycolor','w','xtick',[],'ytick',[]);

annotation('textbox',[0.01 .51 0.01 .51],'String','(a)','FontSize',10,'EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.24 .51 0.25 .51],'String','(b)','FontSize',10,'EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.47 .51 0.48 .51],'String','(c)','FontSize',10,'EdgeColor','none','HorizontalAlignment','center');
clearvars amp fb fb_deep fb_deep_eff fb_deep_sn fb_deep_sn_eff fb_drained_sn fb_drained_eff fb_drained_true fb_drained_sn_eff fb_drained_sn_true fb_drained_unponded fb_in fb_in_sn fb_ponded fb_ponded_sn hi hmp hsn rhosi rhosn rhow w

%% Figure 4: Summer gas and brine volume, ice density, and freeboard
a_mp_cor = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 15 18 16 20]/100; % melt pond fraction from transect, from Webster et al. (2022)
hs_fy(19:23) = [7.5 5 5 5 5]/100.*[1-a_mp_cor(19) 1-a_mp_cor(20) 1-a_mp_cor(21) 1-a_mp_cor(22) 1-a_mp_cor(23)];
CT = T_fy_b; SA(1:20) = 33.7; SA(19:23) = Sb_fy_ct(3:7); rho_sw = gsw_rho(SA,CT,0); % seawater density
rho_fyi_bal = (h_fy_avg/100.*rho_sw - hs_fy.*rho_s - fb_fy.*rho_sw)./(h_fy_avg/100); % FYI density from coring freeboard and draft
vg_tape = 1 - vb_bulk_fy/1000 .* (1 - rho_sw/917.1) - rho_fyi_bal/917.1; % FYI gas volume from coring freeboard, draft, and brine volume
melt = 19:23; premelt = 15:18;

figure
tile = tiledlayout(2,3); tile.TileSpacing = 'compact'; tile.Padding = 'compact';
nexttile % freeboard
t_ALS = datetime(['21-Mar-2020';'08-Apr-2020';'23-Apr-2020';'10-May-2020';'30-Jun-2020';'04-Jul-2020';'17-Jul-2020';'22-Jul-2020']); fb_sn_FYI_ALS = [0.2376 0.2643 0.2589 0.2805 0.206 0.111 0.177 0.192];
sn_FYI_ALS = [0.124 0.124 0.147 0.157 0.08 0.048 0.048 0.048]; fb_FYI_ALS = fb_sn_FYI_ALS - sn_FYI_ALS;
hs_web = hs_fy;
msz = 2.2;
fill([t_fy(19)-15 t_fy(19)-15 t_fy(19)-5 t_fy(19)-5],[1 0 0 1],1,'FaceColor',c{1},'EdgeColor','none') ;alpha(.1); hold on % shade
fill([t_fy(21)-3 t_fy(21)-3 t_fy(21)+1 t_fy(21)+1],[1 0 0 1],1,'FaceColor',c{2},'EdgeColor','none') ;alpha(.1); hold on % shade, melt pond drainage
errorbar(t_fy(19)-10,mean(h_fy_avg(premelt)-d_fy_avg(premelt))/100,std(h_fy_avg(premelt)-d_fy_avg(premelt))/100,'color',[0.72 0.82 0.92]);
errorbar(t_fy(19)-12,mean((h_fy_avg(premelt)-d_fy_avg(premelt))/100+hs_web(premelt)),std((h_fy_avg(premelt)-d_fy_avg(premelt))/100+hs_web(premelt)),'color',[0.72 0.82 0.92]);
errorbar(t_fy(melt),(h_fy_avg(melt)-d_fy_avg(melt))/100,[3.8 3.7 4.5 4.8 2.8]/100,'color',[0.72 0.82 0.92]);
errorbar(t_fy(19)-8,mean(fb_FYI_ALS(1:4)),std(fb_FYI_ALS(1:4)),'color',c{2});
errorbar(t_fy(19)-8,mean(fb_sn_FYI_ALS(1:4)),std(fb_sn_FYI_ALS(1:4)),'color',c{2});
p = plot([t_fy(19)-12 t_fy(melt)],[mean((h_fy_avg(premelt)-d_fy_avg(premelt))/100+hs_web(premelt)) (h_fy_avg(melt)-d_fy_avg(melt))/100+hs_web(melt)],'o:','Color',c{1},'LineWidth',2.5); p.MarkerSize = msz; set(p,'markerfacecolor','w'); hold on
p = plot([t_fy(19)-10 t_fy(19:23)],[mean(h_fy_avg(premelt)-d_fy_avg(premelt))/100 (h_fy_avg(19:23)-d_fy_avg(19:23))/100],'o-','Color',c{1}); p.MarkerSize = msz; set(p,'markerfacecolor',get(p,'color')); hold on
prt = [5 7:8]; p = plot([t_fy(19)-8; t_ALS(prt)],[mean(fb_sn_FYI_ALS(1:4)) fb_sn_FYI_ALS(prt)],':o','Color',c{2},'LineWidth',2.5); p.MarkerSize = msz; set(p,'markerfacecolor','w'); hold on
prt = [5 7:8]; p = plot([t_fy(19)-8; t_ALS(prt)],[mean(fb_FYI_ALS(1:4)) fb_FYI_ALS(prt)],'o-','Color',c{2}); p.MarkerSize = msz; set(p,'markerfacecolor',get(p,'color'));
prt = [5 7:8]; p = plot([t_fy(19)-8; t_ALS(prt)],[mean(fb_sn_FYI_ALS(1:4)) fb_sn_FYI_ALS(prt)],'o','Color',c{2},'LineWidth',1); p.MarkerSize = msz; set(p,'markerfacecolor','w');
p = plot([t_fy(19)-12 t_fy(melt)],[mean((h_fy_avg(premelt)-d_fy_avg(premelt))/100+hs_web(premelt)) (h_fy_avg(melt)-d_fy_avg(melt))/100+hs_web(melt)],'o','Color',c{1},'LineWidth',1); p.MarkerSize = msz; set(p,'markerfacecolor','w'); hold on
% prt = 6; p = plot([t_fy(19)-8; t_ALS(prt)],[mean(fb_sn_FYI_ALS(1:4)) fb_sn_FYI_ALS(prt)],'o','Color',c{2},'LineWidth',1); p.MarkerSize = msz; set(p,'markerfacecolor','w');
% prt = 6; p = plot([t_fy(19)-8; t_ALS(prt)],[mean(fb_FYI_ALS(1:4)) fb_FYI_ALS(prt)],'o','Color',c{2}); p.MarkerSize = msz; set(p,'markerfacecolor',get(p,'color'));
fill([t_fy(19)-5 t_fy(19)-5 t_fy(19)-5+2 t_fy(19)-5+2],[1 0 0 1],1,'FaceColor','w','EdgeColor','none'); hold on % white break in axis
hYLabel = ylabel('FYI freeboard (m)'); set([hYLabel gca],'FontSize',8,'FontWeight','normal'); ylim([0.0 0.40])
leg = legend('','','','','','','','Snow, coring','Ice, coring','Snow, ROV','Ice, ROV','box','off','NumColumns',2); set(leg,'FontSize',7,'Location','best'); leg.ItemTokenSize = [30*0.44,18*0.44];
t_start = datetime('09-Jun-2020'); t_end = datetime('31-Jul-2020'); xlim([t_start t_end]);
ax = gca; ax.XTick = datetime(['08-Jun-2020';'20-Jun-2020';'30-Jun-2020';'10-Jul-2020';'20-Jul-2020';'31-Jul-2020']); datetick('x','dd mmm','keepticks'); xtickangle(0); % time
xticklabels({'1 May', '20 Jun','30 Jun','10 Jul','20 Jul','31 Jul'}); set(gca,'FontSize',7,'FontWeight','normal');
title('Freeboard: coring, ROV sites','FontSize',7,'FontWeight','normal');

nexttile([2 1]) % density
premelt = [5:9 11:19];
rho_als_may = (1020*(1.755+0.157-0.2805)-300*0.157)/(1.755); % 1.755 m - ice thickness assuming 0.235 m melt, 0.157 m level ice snow depth from Itkin, 0.2805 m ALS freeboard, 300 snow density, May 10
rho_als_apr = (1020*(1.73+0.147-0.2589)-300*0.147)/(1.73); % 1.732 m - ice thickness assuming 0.21 m melt, 0.147 m level ice snow depth from Itkin, 0.2589 m ALS freeboard, 300 snow density, Apr 23
rho_als_8_apr = (1020*(1.67+0.124-0.2343)-300*0.124)/(1.67); % 1.67 m - ice thickness assuming 0.15 m melt, 0.124 m level ice snow depth from Itkin, 0.2343 m ALS freeboard, 300 snow density, Apr 8
rho_als_21_mar = (1020*(1.587+0.124-0.2375)-300*0.124)/(1.587); % 1.67 m - ice thickness assuming 0.07 m melt, 0.124 m level ice snow depth from Itkin, 0.2376 m ALS freeboard, 300 snow density, Mar 21
rho_als_premelt = [rho_als_21_mar rho_als_8_apr rho_als_apr rho_als_may]; % pre-melt densities from ALS for March-May 
fill([t_fy(19)-15 t_fy(19)-15 t_fy(19)-5 t_fy(19)-5],[1000 0 0 1000],1,'FaceColor',c{1},'EdgeColor','none') ;alpha(.1); hold on % shade
fill([t_fy(21)-3 t_fy(21)-3 t_fy(21)+1 t_fy(21)+1],[1000 0 0 1000],1,'FaceColor',c{2},'EdgeColor','none') ;alpha(.1); % shade, melt pond drainage
p = errorbar(t_fy(melt),rho_fyi_bal(melt),[22.3 21.1 27.5 33.7 25.2],'o-','Color',[0.72 0.82 0.92]); p.MarkerSize = 2; set(p,'markerfacecolor',get(p,'color'));
errorbar(t_fy(19)-11,nanmean(rho_fyi_bal(premelt)),nanstd(rho_fyi_bal(premelt)),'color',[0.72 0.82 0.92]);
p = plot([t_fy(19)-11 t_fy(melt)],[nanmean(rho_fyi_bal(premelt)) rho_fyi_bal(melt)],'o-','Color',c{1},'linewidth',2); p.MarkerSize = 2; set(p,'markerfacecolor',get(p,'color')); % coring freeboard
p = plot([t_fy(19)-9 t_fy(melt)],[nanmean(rho_si_pr_bulk_fy(premelt)) rho_si_pr_bulk_fy(melt)],'o-','Color',c{3},'linewidth',2); p.MarkerSize = 2; set(p,'markerfacecolor',get(p,'color')); % disconnected pockets
p = plot([t_fy(19)-9 td_int(1) td_int(409) td_int(end)],[mean(rho_als_premelt) rho_si_rov_als(1) rho_si_rov_als(409) rho_si_rov_als(end)],'o-','Color',c{2},'linewidth',2); p.MarkerSize = 2; set(p,'markerfacecolor',get(p,'color'));
p = plot([t_ALS(5) t_ALS(7) t_ALS(8)],rho_si_rov_als_mp,'o:','Color',c{2},'linewidth',2); p.MarkerSize = 2; set(p,'markerfacecolor',get(p,'color')); % filled ALS freeboard + melt ponds from Fuchs, markers
p = plot([t_ALS(5) t_ALS(7) t_ALS(8)],rho_si_rov_als_mp,':','Color',c{2},'linewidth',3); p.MarkerSize = 2; set(p,'markerfacecolor',get(p,'color')); % filled ALS freeboard + melt ponds from Fuchs
p = plot([t_fy(21:23)],[858.7 888.9 882.9],'o:','Color',c{1}); p.MarkerSize = 2.5; set(p,'markerfacecolor','w'); % correction for coring at unponded ice only
errorbar(t_fy(19)-9,mean(rho_als_premelt),std(rho_als_premelt),'color',c{2});
p = text(t_fy(19)-10,870,'pre-melt'); set(p,'Color',c{1},'HorizontalAlignment','right','FontSize',8,'Rotation',90); % pre-melt text
p = text(t_fy(21)-1.9,840,'drainage'); set(p,'Color',c{2},'HorizontalAlignment','right','FontSize',7,'Rotation',90); % drainage text
errorbar(t_fy(19)-9,nanmean(rho_si_pr_bulk_fy(premelt)),nanstd(rho_si_pr_bulk_fy(premelt)),'color',c{3});
fill([t_fy(19)-5 t_fy(19)-5 t_fy(19)-5+2 t_fy(19)-5+2],[1000 800 800 1000],1,'FaceColor','w','EdgeColor','none'); hold on % white break in axis
leg = legend('','','','','Coring','Weighing','ROV site','ROV site, w/ melt ponds','box','off','NumColumns',1); set(leg,'FontSize',7,'Location','best'); leg.ItemTokenSize = [30*0.66,18*0.66];
hYLabel = ylabel('FYI density (kg  m^-^3)'); set([hYLabel gca],'FontSize',8,'FontWeight','normal');
t_start = datetime('09-Jun-2020'); t_end = datetime('31-Jul-2020'); xlim([t_start t_end]);
ax = gca; ax.XTick = datetime(['08-Jun-2020';'20-Jun-2020';'30-Jun-2020';'10-Jul-2020';'20-Jul-2020';'31-Jul-2020']); datetick('x','dd mmm','keepticks'); xtickangle(0); % time
xticklabels({'1 May', '20 Jun','30 Jun','10 Jul','20 Jul','31 Jul'}); set(gca,'FontSize',7,'FontWeight','normal');
yticks(820:20:980); ylim([815 960]);
title('FYI density','FontSize',7,'FontWeight','normal');

nexttile([2 1]) % brine and air volume
melt = 19:23; premelt = [5:9 11:19];
fill([t_fy(19)-15 t_fy(19)-15 t_fy(19)-5 t_fy(19)-5],[100 -10 -10 100],1,'FaceColor',c{1},'EdgeColor','none') ;alpha(.1); hold on % shade, pre-melt
fill([t_fy(21)-3 t_fy(21)-3 t_fy(21)+1 t_fy(21)+1],[100 -10 -10 100],1,'FaceColor',c{2},'EdgeColor','none') ;alpha(.1); % shade, melt pond drainage
errorbar(t_fy(melt),vg_tape(melt)*100,[2.43 2.30 3.00 3.67 2.75],'color',[0.72 0.82 0.92]);
errorbar(t_fy(19)-9,nanmean(vg_tape(premelt)*100),nanstd(vg_tape(premelt)*100),'color',[0.72 0.82 0.92]);
p = plot([t_fy(19)-11 t_fy(melt)],[nanmean(vb_bulk_fy(premelt)/10) vb_bulk_fy(melt)/10],'o:','color',c{3},'linewidth',2.0); p.MarkerSize = 2.5; set(p,'markerfacecolor','w'); hold on
p = plot([t_fy(19)-9 t_fy(melt)],[nanmean(vg_tape(premelt)*100) max(vg_tape(melt)*100,0)],'o-','Color',c{1}); p.MarkerSize = 2; set(p,'markerfacecolor',get(p,'color'));
p = plot([t_fy(19)-11 t_fy(melt)],[nanmean(vg_pr_bulk_fy(premelt)/10) vg_pr_bulk_fy(melt)/10],'o-','color',c{3}); p.MarkerSize = 2; set(p,'markerfacecolor',get(p,'color'));
p = plot(t_fy(melt),max(vg_tape(melt)*100,0),'o','Color',c{1}); p.MarkerSize = 2; set(p,'markerfacecolor',get(p,'color'));
errorbar(t_fy(19)-11,nanmean(vg_pr_bulk_fy(premelt)/10),nanstd(vg_pr_bulk_fy(premelt)/10),'color',c{3});
errorbar(t_fy(19)-11,nanmean(vb_bulk_fy(premelt)/10),nanstd(vb_bulk_fy(premelt)/10),'color',c{3});
p = text(t_fy(19)-10,20,'pre-melt'); set(p,'Color',c{1},'HorizontalAlignment','right','FontSize',8,'Rotation',90); % pre-melt text
p = text(t_fy(21)-1.3,2.7,'drainage'); set(p,'Color',c{2},'HorizontalAlignment','right','FontSize',7,'Rotation',90); % pre-melt text
fill([t_fy(19)-5 t_fy(19)-5 t_fy(19)-5+2 t_fy(19)-5+2],[35 -5 -5 35],1,'FaceColor','w','EdgeColor','none'); hold on % white break in axis
leg = legend('','','','','Brine','Air, coring','Air, weighing','box','off','NumColumns',1); set(leg,'FontSize',7,'Location','best'); leg.ItemTokenSize = [30*0.66,18*0.66];
hYLabel = ylabel('FYI volume (%)'); set([hYLabel gca],'FontSize',8,'FontWeight','normal'); ylim([-4 31])
t_start = datetime('09-Jun-2020'); t_end = datetime('31-Jul-2020'); xlim([t_start t_end]);
ax = gca; ax.XTick = datetime(['08-Jun-2020';'20-Jun-2020';'30-Jun-2020';'10-Jul-2020';'20-Jul-2020';'31-Jul-2020']); datetick('x','dd mmm','keepticks'); xtickangle(0); % time
xticklabels({'1 May', '20 Jun','30 Jun','10 Jul','20 Jul','31 Jul'}); set(gca,'FontSize',7,'FontWeight','normal');
title('Brine and air volume','FontSize',7,'FontWeight','normal');

nexttile % ice draft
time = ncread('IMB_CMW_MOSAiC_density_premelt.nc','time'); hi_max_imb = ncread('IMB_CMW_MOSAiC_density_premelt.nc','hi_max'); bot = ncread('IMB_CMW_MOSAiC_density_premelt.nc','bot'); % Import of IMB data
hsn_max_imb = ncread('IMB_CMW_MOSAiC_density_premelt.nc','hsn_max'); sno = ncread('IMB_CMW_MOSAiC_density_premelt.nc','sno'); T_bot = ncread('IMB_CMW_MOSAiC_density_premelt.nc','T_bot'); ice_type = ncread('IMB_CMW_MOSAiC_density_premelt.nc','imb_ice_type');
sur = ncread('IMB_CMW_MOSAiC_density_premelt.nc','sur'); t_0 = (datetime('1979-01-01 00:00:00')); t_IMB = t_0 + days(time); % Import of IMB data
ice_thick_imb = -ncread('IMB_CMW_MOSAiC_density_premelt.nc','ice_thick'); snow_thick_imb = -ncread('IMB_CMW_MOSAiC_density_premelt.nc','snow_thick');
rho_fyi_imb = interp1(convertTo(t_fy,"datenum"),rho_si_pr_bulk_fy,convertTo(t_IMB,"datenum"),'linear','extrap'); % interpolation of FYI bulk density for IMB
rho_syi_imb = interp1(convertTo(t_sy,"datenum"),rho_si_pr_bulk_sy,convertTo(t_IMB,"datenum"),'linear','extrap'); % interpolation of SYI bulk density for IMB
t_sp = datetime(['03-May-2020';'22-Jun-2020';'06-Jul-2020';'13-Jul-2020';'20-Jul-2020';'27-Jul-2020']); % snow pits, June-July 2020
rho_sn_sp =     [         295         416.7         420.5         289.5         427.5          389.5];
rho_sn_imb = interp1(convertTo(t_sp,"datenum"),rho_sn_sp,convertTo(t_IMB,"datenum"),'linear','extrap'); % interpolation of snow bulk density for IMB
good = [1:6 9]; Tw_imb = mean(T_bot(good,:)); Sw_imb = gsw_SA_freezing_from_t(Tw_imb,0,0); rho_sw_imb = gsw_rho(Sw_imb,Tw_imb,0);
good = [1:6 9]; hi_imb = mean(ice_thick_imb(good,:),1); % IMB ice thickness
good = [1:6 9]; hs_imb = mean(snow_thick_imb(good,:),1); % IMB snow thickness
d_imb = -(rho_fyi_imb.*hi_imb'+rho_sn_imb.*hs_imb')./rho_sw_imb';
good = [4 5 9]; Tw_imb = mean(T_bot(good,:)); Sw_imb = gsw_SA_freezing_from_t(Tw_imb,0,0); rho_fyi_sw_imb = gsw_rho(Sw_imb,Tw_imb,0);
good = [4 5 9]; hi_fyi_imb = mean(ice_thick_imb(good,:),1); % IMB ice thickness
good = [4 5 9]; hs_fyi_imb = mean(snow_thick_imb(good,:),1); % IMB snow thickness
d_fyi_imb = -(rho_fyi_imb.*hi_fyi_imb'+rho_sn_imb.*hs_fyi_imb')./rho_fyi_sw_imb';
good = [1:3 6]; Tw_imb = mean(T_bot(good,:)); Sw_imb = gsw_SA_freezing_from_t(Tw_imb,0,0); rho_syi_sw_imb = gsw_rho(Sw_imb,Tw_imb,0);
good = [1:3 6]; hi_syi_imb = mean(ice_thick_imb(good,:),1); % IMB ice thickness
good = [1:3 6]; hs_syi_imb = mean(snow_thick_imb(good,:),1); % IMB snow thickness
d_syi_imb = -(rho_syi_imb.*hi_syi_imb'+rho_sn_imb.*hs_syi_imb')./rho_syi_sw_imb';
fill([t_fy(19)-15 t_fy(19)-15 t_fy(19)-5 t_fy(19)-5],[1 -2 -2 1],1,'FaceColor',c{1},'EdgeColor','none') ;alpha(.1); hold on % shade
fill([t_fy(21)-3 t_fy(21)-3 t_fy(21)+1 t_fy(21)+1],[1 -2 -2 1],1,'FaceColor',c{2},'EdgeColor','none') ;alpha(.1); hold on % shade, melt pond drainage
errorbar(t_fy(19)-11,nanmean(-d_fy_avg(15:18)/100),nanstd(-d_fy_avg(15:18)/100),'color',[0.72 0.82 0.92]);
p = plot(t_fy(19)-11,nanmean(-d_fy_avg(15:18)/100),'o','color',c{1}); p.MarkerSize = 2; set(p,'markerfacecolor',get(p,'color'));
p = plot(t_fy(15:23),-d_fy_avg(15:23)/100,'o','color',c{1}); p.MarkerSize = 2; set(p,'markerfacecolor',get(p,'color'));
p = plot(t_ROV,d_FYI_ROV,'o','color',c{2}); p.MarkerSize = 2; set(p,'markerfacecolor',get(p,'color'));
plot(t_IMB(1753+50:end),d_fyi_imb(1753+50:end),':','color',c{1},'LineWidth',2.5); % plot(t_IMB,d_imb,'--','color',c{3},'LineWidth',1.0);
p = text(t_fy(21)-1.3,-1.37,'drainage'); set(p,'Color',c{2},'HorizontalAlignment','right','FontSize',6,'Rotation',90); % drainage text
fill([t_fy(19)-5 t_fy(19)-5 t_fy(19)-5+2 t_fy(19)-5+2],[0 -2 -2 0],1,'FaceColor','w','EdgeColor','none'); hold on % white break in axis
leg = legend('','','','','Coring site','ROV site','IMB','box','off','NumColumns',1); set(leg,'FontSize',7,'Location','northwest'); leg.ItemTokenSize = [30*0.5,18*0.5];
hYLabel = ylabel('FYI draft (m)'); set([hYLabel gca],'FontSize',8,'FontWeight','normal'); ylim([-1.75 -0.6])
t_start = datetime('09-Jun-2020'); t_end = datetime('31-Jul-2020'); xlim([t_start t_end]);
ax = gca; ax.XTick = datetime(['08-Jun-2020';'20-Jun-2020';'30-Jun-2020';'10-Jul-2020';'20-Jul-2020';'31-Jul-2020']); datetick('x','dd mmm','keepticks'); xtickangle(0); % time
xticklabels({'1 May', '20 Jun','30 Jun','10 Jul','20 Jul','31 Jul'}); set(gca,'FontSize',7,'FontWeight','normal');
title('Ice draft','FontSize',7,'FontWeight','normal');

axes('Position',[.11 .549 .0065 .017]); % axis break (a)
px=[0 0]; py=[0 1]; w=2; px2=px+w; fill([px flip(px2)],[py flip(py)],'w','EdgeColor','none'); hold all;
plot(px,py,'k',px2,py,'k','linewidth',0.1); box off; axis off;
axes('Position',[.11 .061 .0065 .02]); % axis break (d)
px=[0 0]; py=[0 1]; w=2; px2=px+w; fill([px flip(px2)],[py flip(py)],'w','EdgeColor','none'); hold all;
plot(px,py,'k',px2,py,'k','linewidth',0.1); box off; axis off;
axes('Position',[.43 .061 .0065 .02]); % axis break (b)
px=[0 0]; py=[0 1]; w=2; px2=px+w; fill([px flip(px2)],[py flip(py)],'w','EdgeColor','none'); hold all;
plot(px,py,'k',px2,py,'k','linewidth',0.1); box off; axis off;
axes('Position',[.75 .061 .0065 .02]); % axis break (c)
px=[0 0]; py=[0 1]; w=2; px2=px+w; fill([px flip(px2)],[py flip(py)],'w','EdgeColor','none'); hold all;
plot(px,py,'k',px2,py,'k','linewidth',0.1); box off; axis off;

annotation('textbox',[0.01 .51 0.01 .51],'String','(a)','FontSize',8,'EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.22 .51 0.22 .51],'String','(b)','FontSize',8,'EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.43 .51 0.44 .51],'String','(c)','FontSize',8,'EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.01 .28 0.01 .28],'String','(d)','FontSize',8,'EdgeColor','none','HorizontalAlignment','center');

%% Figure 5: freeboard, draft, thickness at different scales (6.2 x 2.5 in)
a_mp_cor = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 15 18 16 20]/100; % melt pond area from Webster et al. during coring events
hs_fy(19:23) = [7.5 5 5 5 5]/100.*[1-a_mp_cor(19) 1-a_mp_cor(20) 1-a_mp_cor(21) 1-a_mp_cor(22) 1-a_mp_cor(23)]; % snow thickness from Webster et al. during coring events
melt = 19:23; premelt = [5:9 11:19];

figure
tile = tiledlayout(2,4); tile.TileSpacing = 'compact'; tile.Padding = 'compact';
nexttile % freeboard
t_ALS = datetime(['21-Mar-2020';'08-Apr-2020';'23-Apr-2020';'10-May-2020';'30-Jun-2020';'04-Jul-2020';'17-Jul-2020';'22-Jul-2020']); fb_sn_FYI_ALS = [0.2376 0.2643 0.2589 0.2805 0.206 0.111 0.177 0.192];
sn_FYI_ALS = [0.124 0.124 0.147 0.157 0.08 0.048 0.048 0.048]; fb_FYI_ALS = fb_sn_FYI_ALS - sn_FYI_ALS;
hs_web = hs_fy;
msz = 1.9;
fill([t_fy(15)-15 t_fy(15)-15 t_fy(18)+6 t_fy(18)+6],[1 -2 -2 1],1,'FaceColor',c{1},'EdgeColor','none'); alpha(.1); hold on % shade
fill([t_fy(21)-3 t_fy(21)-3 t_fy(21)+1 t_fy(21)+1],[1 0 0 1],1,'FaceColor',c{2},'EdgeColor','none') ;alpha(.1); % shade, melt pond drainage
prt = [5 7:8]; p = plot([t_ALS(1:4); t_ALS(prt)],[fb_sn_FYI_ALS(1:4) fb_sn_FYI_ALS(prt)],'o:','Color',c{2},'LineWidth',2.5); p.MarkerSize = msz; set(p,'markerfacecolor','w'); % ALS ROV snow freeboard
prt = [5 7:8]; p = plot([t_ALS(1:4); t_ALS(prt)],[fb_FYI_ALS(1:4) fb_FYI_ALS(prt)],'o-','Color',c{2}); p.MarkerSize = msz; set(p,'markerfacecolor',get(p,'color')); % ALS ROV ice freeboard
p = plot([t_fy(15:18) t_fy(melt)],[(h_fy_avg(15:18)-d_fy_avg(15:18))/100 (h_fy_avg(19:23)-d_fy_avg(19:23))/100],'o-','Color',c{1},'LineWidth',0.5); p.MarkerSize = msz; set(p,'markerfacecolor',get(p,'color')); % Ice, FYI coring
prt = [5 7:8]; p = plot([t_ALS(1:4); t_ALS(prt)],[fb_sn_FYI_ALS(1:4) fb_sn_FYI_ALS(prt)],'o','Color',c{2},'LineWidth',1); p.MarkerSize = msz; set(p,'markerfacecolor','w');
p = plot([t_fy(15:18) t_fy(melt)],[((h_fy_avg(15:18)-d_fy_avg(15:18))/100) (h_fy_avg(melt)-d_fy_avg(melt))/100],'o','Color',c{1},'LineWidth',.1); p.MarkerSize = msz; set(p,'markerfacecolor',get(p,'color'));
p = text(t_fy(15)+16,0.47,'pre-melt'); set(p,'Color',c{1},'HorizontalAlignment','right','FontSize',7,'Rotation',90); % text
hYLabel = ylabel('Freeboard (m)'); set([hYLabel gca],'FontSize',8,'FontWeight','normal'); ylim([0.0 0.50]); yticks(0:0.1:0.6); 
leg = legend('','','Snow, ROV site','Ice, ROV site','Ice, coring site','box','off','NumColumns',1); set(leg,'FontSize',6,'Location','best'); leg.ItemTokenSize = [30*0.44,18*0.44];
ax = gca; ax.XTick = (datetime(['01-Apr-2020';'01-May-2020';'01-Jun-2020';'01-Jul-2020';'01-Aug-2020'])); datetick('x','mmm','keepticks'); xtickangle(0); % time
datetick('x','mmm','keepticks','keeplimits');
title('Freeboard: coring, ROV sites','FontSize',7,'FontWeight','normal'); set(gca,'FontSize',7,'FontWeight','normal');

t_als(1) = datetime('08-Apr-2020'); t_als(2) = datetime('23-Apr-2020'); t_als(3) = datetime('10-May-2020');  t_als(4) = datetime('16-Jun-2020'); t_als(5) = datetime('30-Jun-2020'); t_als(6) = datetime('4-Jul-2020');
t_als(7) = datetime('7-Jul-2020'); t_als(8) = datetime('17-Jul-2020'); t_als(9) = datetime('22-Jul-2020');
h_sn_web = [0.284 0.300 0.311 0.182 0.10 0.068 0.052 0.064 0.061]; % snow or SSL depth from Webster (2022), all ice
h_sn_level_web = [0.263 0.263 0.279 0.182 0.10 0.068 0.052 0.064 0.061]; % snow or SSL depth from Webster (2022), hi < 2.5 m
fb_all_mean = [0.391 0.422 0.446 0.3823 0.3102 0.3044 0.3152 0.2970 0.2638];
fb_level_mean = [0.322 0.352 0.371 0.3168 0.2534 0.2462 0.2576 0.2385 0.2279];
fb_co_mean = [0.5471 0.5403 0.6106 NaN  0.3629 0.2934 0.3311 0.3168 0.3176];
fb_co_level_mean = [0.456 0.4581 0.4660 NaN 0.2746 0.2284 0.2535 0.2515 0.2508];
fb_co_fyi_level_mean = [0.3573 0.3516 0.3903 NaN 0.2225 0.1785 0.1946 0.2135 0.2071];
fb_co_syi_level_mean = [0.5060 0.5116 0.4895 NaN 0.3015 0.2593 0.2850 0.2714 0.2727];

nexttile([2 1]) % ALS CO2
fill([t_als(1)-7 t_als(1)-7 t_als(3) t_als(3)],[1 0 0 1],1,'FaceColor',c{1},'EdgeColor','none') ;alpha(.1); hold on % shade
fill([t_als(7)+3 t_als(7)+3 t_als(7)+7 t_als(7)+7],[1 0 0 1],1,'FaceColor',c{2},'EdgeColor','none') ;alpha(.1); % shade, melt pond drainage
good = [1 2 3 5:9]; msz = 2.0;
p = plot(t_als(good),fb_co_mean(good),':','color',c{4},'linewidth',2); p.MarkerSize = msz; set(p,'markerfacecolor','w');
p = plot(t_als(good),fb_co_mean(good),'o:','color',c{4},'linewidth',1); p.MarkerSize = msz; set(p,'markerfacecolor','w');
p = plot(t_als(good),fb_co_level_mean(good),':','color',c{3},'linewidth',2); p.MarkerSize = msz; set(p,'markerfacecolor',get(p,'color'));
p = plot(t_als(good),fb_co_level_mean(good),'o:','color',c{3},'linewidth',1); p.MarkerSize = msz; set(p,'markerfacecolor','w');
p = plot(t_als(good),fb_co_syi_level_mean(good)-h_sn_web(good),'go-','color',c{2}); p.MarkerSize = msz; set(p,'markerfacecolor',get(p,'color'));
p = plot(t_als(good),fb_co_level_mean(good)-h_sn_web(good),'o-','color',c{3}); p.MarkerSize = msz; set(p,'markerfacecolor',get(p,'color'));
p = plot(t_als(good),fb_co_fyi_level_mean(good)-h_sn_web(good),'mo-','color',c{1}); p.MarkerSize = msz; set(p,'markerfacecolor',get(p,'color'));
p = text(t_als(3)-19,0.40,'pre-melt'); set(p,'Color',c{1},'HorizontalAlignment','right','FontSize',8,'Rotation',90); % text
p = text(t_als(7)+5,0.13,'drainage'); set(p,'Color',c{2},'HorizontalAlignment','right','FontSize',7,'Rotation',90); % drainage text
hYLabel = ylabel('Freeboard (m)'); set([hYLabel gca],'FontSize',7,'FontWeight','normal');
leg = legend('','','','Snow','','Snow, level','Ice, level SYI','Ice, level','Ice, level FYI','box','off','NumColumns',1); set(leg,'FontSize',7,'Location','best'); leg.ItemTokenSize = [30*0.45,18*0.45];
ax = gca; ax.XTick = (datetime(['01-Apr-2020';'01-May-2020';'01-Jun-2020';'01-Jul-2020';'01-Aug-2020'])); datetick('x','mmm','keepticks'); xtickangle(0); % time
datetick('x','mmm','keepticks','keeplimits');
title('Freeboard: CO2','FontSize',8,'FontWeight','normal'); ylim([0 .25]); ylim([0 .63]); set(gca,'FontSize',7,'FontWeight','normal');

nexttile([2 1]) % ALS full
fill([t_als(1)-7 t_als(1)-7 t_als(3) t_als(3)],[1 0 0 1],1,'FaceColor',c{1},'EdgeColor','none') ;alpha(.1); hold on % shade
fill([t_als(7)+3 t_als(7)+3 t_als(7)+7 t_als(7)+7],[1 0 0 1],1,'FaceColor',c{2},'EdgeColor','none') ;alpha(.1); % shade, melt pond drainage
p = plot(t_als,fb_all_mean,':','color',c{4},'linewidth',2); p.MarkerSize = msz; set(p,'markerfacecolor',get(p,'color'));
p = plot(t_als,fb_all_mean,'o:','color',c{4},'linewidth',1); p.MarkerSize = msz; set(p,'markerfacecolor','w');
p = plot(t_als,fb_level_mean,':','color',c{3},'linewidth',2); p.MarkerSize = msz; set(p,'markerfacecolor',get(p,'color'));
p = plot(t_als,fb_level_mean,'o:','color',c{3},'linewidth',1); p.MarkerSize = msz; set(p,'markerfacecolor','w');
p = plot(t_als,fb_all_mean - h_sn_web,'o-','color',c{4}); p.MarkerSize = msz; set(p,'markerfacecolor',get(p,'color'));
p = plot(t_als,fb_level_mean - h_sn_level_web,'o-','color',c{3}); p.MarkerSize = msz; set(p,'markerfacecolor',get(p,'color'));
p = text(t_als(3)-19,0.30,'pre-melt'); set(p,'Color',c{1},'HorizontalAlignment','right','FontSize',8,'Rotation',90); % text
p = text(t_als(7)+5,0.13,'drainage'); set(p,'Color',c{2},'HorizontalAlignment','right','FontSize',7,'Rotation',90); % drainage text
hYLabel = ylabel('Freeboard (m)'); set([hYLabel gca],'FontSize',7,'FontWeight','normal');
leg = legend('','','','Snow','','Snow, level','Ice','Ice, level','box','off','NumColumns',1); set(leg,'FontSize',7,'Location','best'); leg.ItemTokenSize = [30*0.45,18*0.45];
ax = gca; ax.XTick = (datetime(['01-Apr-2020';'01-May-2020';'01-Jun-2020';'01-Jul-2020';'01-Aug-2020'])); datetick('x','mmm','keepticks'); xtickangle(0); % time
datetick('x','mmm','keepticks','keeplimits');
title('Freeboard: ALS full','FontSize',8,'FontWeight','normal'); ylim([0 .3]); ylim([0 .63]); set(gca,'FontSize',7,'FontWeight','normal');

nexttile([2 1]) % Accumulated ice melt
time = ncread('IMB_CMW_MOSAiC_density_new.nc','time'); bot = ncread('IMB_CMW_MOSAiC_density_new.nc','bot'); sur = ncread('IMB_CMW_MOSAiC_density_new.nc','sur'); t_0 = (datetime('1979-01-01 00:00:00')); t_IMB = t_0 + days(time); % Import of IMB data
t_ROV = datetime(['24-Jun-2020 12:00';'28-Jun-2020 12:00';'01-Jul-2020 12:00';'07-Jul-2020 12:00';'14-Jul-2020 12:00';'21-Jul-2020 12:00';'28-Jul-2020 12:00']);
h_FYI_ROV = [1.5653 1.5351 1.5213 1.5162 1.2766 1.1496 1.0008]; h_SYI_ROV = [2.9080 2.9008 2.7690 2.9315 2.2971 2.3560 1.9792]; h_JR_ROV  = [4.2998 4.1892 4.1315 3.9847 3.6286 3.3473 2.5271];
fill([t_ROV(4)+3 t_ROV(4)+3 t_ROV(4)+7 t_ROV(4)+7],[-1 2 2 -1],1,'FaceColor',c{2},'EdgeColor','none') ;alpha(.1); hold on % shade, melt pond drainage
p = plot(t_ROV,h_JR_ROV(1)-h_JR_ROV,'o-','color',c{5}); p.MarkerSize = 2; set(p,'markerfacecolor',get(p,'color')); % total melt, Jaridge
p = plot(t_ROV,h_SYI_ROV(1)-h_SYI_ROV,'o-','color',c{2}); p.MarkerSize = 2; set(p,'markerfacecolor',get(p,'color')); % Total melt, SYI (ROV)
good = [1:6 9]; plot(t_IMB(181:end),mean(bot(good,181:end),1)/0.9+mean(sur(good,181:end),1)/0.9-mean(bot(good,397),1)/0.9-mean(sur(good,397),1)/0.9,'color',c{3},'LineWidth',1.0); hold on % FYI total melt, IMB
p = plot(t_ROV,h_FYI_ROV(1)-h_FYI_ROV,'o--','color',c{1}); p.MarkerSize = 2; set(p,'markerfacecolor','w'); % Total melt, FYI (ROV)
p = plot(t_fy(19:23),h_fy_avg(19)/100-h_fy_avg(19:23)/100,'o-','color',c{1}); p.MarkerSize = 2; set(p,'markerfacecolor',get(p,'color')); % FYI total melt, coring
p = text(t_ROV(4)+4.6,1.0,'drainage'); set(p,'Color',c{2},'HorizontalAlignment','right','FontSize',7,'Rotation',90); % drainage text
leg = legend('','Ridge, ROV','SYI, ROV','FYI & SYI, IMB','FYI, ROV','FYI, coring','box','off'); set(leg,'FontSize',7,'Location','northwest'); leg.ItemTokenSize = [30*0.55,18*0.55];
hYLabel = ylabel('Accumulated ice melt (m)'); set([hYLabel gca],'FontSize',7,'FontWeight','normal'); grid off;
t_start = datetime('15-Jun-2020'); t_end = datetime('31-Jul-2020'); xlim([t_start t_end]); xData = linspace(t_start,t_end,4); ax = gca; ax.XTick = xData; xtickangle(0);
datetick('x','dd mmm','keepticks','keeplimits'); title('Accumulated ice melt','FontSize',8,'FontWeight','normal'); set(gca,'FontSize',7,'FontWeight','normal');
yticks(-0.2:0.2:2.2); ylim([-0.1 1.81]);

nexttile % ice draft
time = ncread('IMB_CMW_MOSAiC_density_premelt.nc','time'); hi_max_imb = ncread('IMB_CMW_MOSAiC_density_premelt.nc','hi_max'); bot = ncread('IMB_CMW_MOSAiC_density_premelt.nc','bot'); % Import of IMB data
hsn_max_imb = ncread('IMB_CMW_MOSAiC_density_premelt.nc','hsn_max'); sno = ncread('IMB_CMW_MOSAiC_density_premelt.nc','sno'); T_bot = ncread('IMB_CMW_MOSAiC_density_premelt.nc','T_bot'); ice_type = ncread('IMB_CMW_MOSAiC_density_premelt.nc','imb_ice_type');
sur = ncread('IMB_CMW_MOSAiC_density_premelt.nc','sur'); t_0 = (datetime('1979-01-01 00:00:00')); t_IMB = t_0 + days(time); % Import of IMB data
ice_thick_imb = -ncread('IMB_CMW_MOSAiC_density_premelt.nc','ice_thick'); snow_thick_imb = -ncread('IMB_CMW_MOSAiC_density_premelt.nc','snow_thick');
rho_fyi_imb = interp1(convertTo(t_fy,"datenum"),rho_si_pr_bulk_fy,convertTo(t_IMB,"datenum"),'linear','extrap'); % interpolation of FYI bulk density for IMB
rho_syi_imb = interp1(convertTo(t_sy,"datenum"),rho_si_pr_bulk_sy,convertTo(t_IMB,"datenum"),'linear','extrap'); % interpolation of SYI bulk density for IMB
t_sp = datetime(['03-May-2020';'22-Jun-2020';'06-Jul-2020';'13-Jul-2020';'20-Jul-2020';'27-Jul-2020']); % snow pits, June-July 2020
rho_sn_sp =     [         295         416.7         420.5         289.5         427.5          389.5];
rho_sn_imb = interp1(convertTo(t_sp,"datenum"),rho_sn_sp,convertTo(t_IMB,"datenum"),'linear','extrap'); % interpolation of snow bulk density for IMB
good = [1:6 9]; Tw_imb = mean(T_bot(good,:)); Sw_imb = gsw_SA_freezing_from_t(Tw_imb,0,0); rho_sw_imb = gsw_rho(Sw_imb,Tw_imb,0);
good = [1:6 9]; hi_imb = mean(ice_thick_imb(good,:),1); % IMB ice thickness
good = [1:6 9]; hs_imb = mean(snow_thick_imb(good,:),1); % IMB snow thickness
d_imb = -(rho_fyi_imb.*hi_imb'+rho_sn_imb.*hs_imb')./rho_sw_imb';
good = [4 5 9]; Tw_imb = mean(T_bot(good,:)); Sw_imb = gsw_SA_freezing_from_t(Tw_imb,0,0); rho_fyi_sw_imb = gsw_rho(Sw_imb,Tw_imb,0);
good = [4 5 9]; hi_fyi_imb = mean(ice_thick_imb(good,:),1); % IMB ice thickness
good = [4 5 9]; hs_fyi_imb = mean(snow_thick_imb(good,:),1); % IMB snow thickness
d_fyi_imb = -(rho_fyi_imb.*hi_fyi_imb'+rho_sn_imb.*hs_fyi_imb')./rho_fyi_sw_imb';
good = [1:3 6]; Tw_imb = mean(T_bot(good,:)); Sw_imb = gsw_SA_freezing_from_t(Tw_imb,0,0); rho_syi_sw_imb = gsw_rho(Sw_imb,Tw_imb,0);
good = [1:3 6]; hi_syi_imb = mean(ice_thick_imb(good,:),1); % IMB ice thickness
good = [1:3 6]; hs_syi_imb = mean(snow_thick_imb(good,:),1); % IMB snow thickness
d_syi_imb = -(rho_syi_imb.*hi_syi_imb'+rho_sn_imb.*hs_syi_imb')./rho_syi_sw_imb';
fill([t_fy(15)-15 t_fy(15)-15 t_fy(18)+6 t_fy(18)+6],[1 -2 -2 1],1,'FaceColor',c{1},'EdgeColor','none') ;alpha(.1); hold on % shade
fill([t_fy(21)-3 t_fy(21)-3 t_fy(21)+1 t_fy(21)+1],[1 -2 -2 1],1,'FaceColor',c{2},'EdgeColor','none') ;alpha(.1); hold on % shade, melt pond drainage
p = plot(t_fy(15:23),-d_fy_avg(15:23)/100,'o','color',c{1}); p.MarkerSize = 2; set(p,'markerfacecolor',get(p,'color'));
p = plot(t_ROV,d_FYI_ROV,'o','color',c{2}); p.MarkerSize = 2; set(p,'markerfacecolor',get(p,'color'));
plot(t_IMB,d_fyi_imb,':','color',c{1},'LineWidth',2.5); % plot(t_IMB,d_imb,'--','color',c{3},'LineWidth',1.0);
plot(t_IMB,d_syi_imb,':','color',c{2},'LineWidth',2.5);
p = text(t_fy(21)-1.3,-1.40,'drainage'); set(p,'Color',c{2},'HorizontalAlignment','right','FontSize',6,'Rotation',90); % drainage text
leg = legend('','','FYI, coring site','FYI, ROV site','FYI, IMB','SYI, IMB','box','off','NumColumns',1); set(leg,'FontSize',7,'Location','best'); leg.ItemTokenSize = [30*0.5,18*0.5];
hYLabel = ylabel('Draft (m)'); set([hYLabel gca],'FontSize',8,'FontWeight','normal'); ylim([-1.8 -0.6]); set(gca,'FontSize',7,'FontWeight','normal');
t_start = datetime('01-Apr-2020'); t_end = datetime('01-Aug-2020'); xlim([t_start t_end]); datetick('x','mmm','keepticks'); xtickangle(0);
set(gca,'FontSize',7,'FontWeight','normal');
title('Ice draft','FontSize',8,'FontWeight','normal'); set(gca,'FontSize',7,'FontWeight','normal');

annotation('textbox',[0.01 .51 0.01 .51],'String','(a)','FontSize',7,'EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.17 .51 0.17 .51],'String','(b)','FontSize',7,'EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.33 .51 0.33 .51],'String','(c)','FontSize',7,'EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.49 .51 0.50 .51],'String','(d)','FontSize',7,'EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.01 .26 0.01 .26],'String','(e)','FontSize',7,'EdgeColor','none','HorizontalAlignment','center');

%% Figure 7: Bulk parameters for all ice types
figure
tile = tiledlayout(4,2); tile.TileSpacing = 'compact'; tile.Padding = 'compact';
nexttile % snow thickness
p = plot(t_fy,hs_fy,'o','LineWidth',0.5,'color',c{1}); p.MarkerSize = 1.5; set(p,'markerfacecolor',get(p,'color')); hold on  % FYI coring
p = plot(t_sy,hs_sy/100,'o','Color',c{2}); p.MarkerSize = 1.5; set(p,'markerfacecolor',get(p,'color')); % SYI coring
p = plot([t_tr_old; t_tr(1); t_tr],[h_sn_tr_old NaN h_sn_tr],'-','color',c{3},'LineWidth',1); p.MarkerSize = 2; set(p,'markerfacecolor',get(p,'color')); % Level ice, transect
p = plot(t_fy,hs_fy,'o','LineWidth',0.5,'color',c{1}); p.MarkerSize = 1.5; set(p,'markerfacecolor',get(p,'color')); hold on  % FYI coring
leg = legend('FYI','SYI','Transect, level ice','box','off','NumColumns',3); set(leg,'FontSize',7,'Location','best'); leg.ItemTokenSize = [30*0.33,18*0.33];
hYLabel = ylabel('Snow thickness (m)'); set(gca,'FontSize',9,'FontWeight','normal'); set(hYLabel,'FontSize',8,'FontWeight','normal'); ylim([0 0.30]); yticks(0:0.10:0.3);
ax = gca; ax.XTick = datetime(['01-Oct-2019';'01-Dec-2019';'01-Feb-2020';'01-Apr-2020';'01-Jun-2020';'01-Aug-2020']); datetick('x','mmm','keepticks'); xtickangle(0); % time

nexttile % salinity
p = plot(t_fy,S_fy_bulk,'o-','LineWidth',0.5,'color',c{1}); p.MarkerSize = 1.5; set(p,'markerfacecolor',get(p,'color')); hold on
p = plot(t_sy,S_sy_bulk,'o-','LineWidth',0.5,'color',c{2}); p.MarkerSize = 1.5; set(p,'markerfacecolor',get(p,'color'));
leg = legend('FYI','SYI','box','off','NumColumns',1); set(leg,'FontSize',7,'Location','best'); leg.ItemTokenSize = [30*0.5,18*0.5];
hYLabel = ylabel('Salinity'); set(gca,'FontSize',9,'FontWeight','normal'); set(hYLabel,'FontSize',8,'FontWeight','normal'); yticks(0:2:8); ylim([0 8]);
ax = gca; ax.XTick = datetime(['01-Oct-2019';'01-Dec-2019';'01-Feb-2020';'01-Apr-2020';'01-Jun-2020';'01-Aug-2020']); datetick('x','mmm','keepticks'); xtickangle(0); % time

nexttile % freeboard
p = plot(t_fy,fb_fy,'o-','LineWidth',0.5,'color',c{1}); p.MarkerSize = 1.5; set(p,'markerfacecolor',get(p,'color')); hold on
p = plot(t_sy,fb_sy,'o-','LineWidth',0.5,'color',c{2}); p.MarkerSize = 1.5; set(p,'markerfacecolor',get(p,'color')); hold on
leg = legend('FYI','SYI','box','off','NumColumns',1); set(leg,'FontSize',7,'Location','best'); leg.ItemTokenSize = [30*0.33,18*0.33];
hYLabel = ylabel('Freeboard (m)'); set(gca,'FontSize',9,'FontWeight','normal'); set(hYLabel,'FontSize',8,'FontWeight','normal'); ylim([0 0.30]); yticks(0:0.10:0.3);
ax = gca; ax.XTick = datetime(['01-Oct-2019';'01-Dec-2019';'01-Feb-2020';'01-Apr-2020';'01-Jun-2020';'01-Aug-2020']); datetick('x','mmm','keepticks'); xtickangle(0); % time

nexttile % temperature
p = plot(t_fy,T_fy_bulk,'o','LineWidth',0.5,'color',c{1}); p.MarkerSize = 1.5; set(p,'markerfacecolor',get(p,'color')); hold on
p = plot(t_sy,T_sy_bulk,'o','LineWidth',0.5,'color',c{2}); p.MarkerSize = 1.5; set(p,'markerfacecolor',get(p,'color'));
plot(t_T66(1:end),mean(frz_T66(1:end,:),2,'omitnan'),'color',[0.92 0.82 0.82],'linewidth',0.1); hold on % FYI bulk IMB
plot(t_T62(1:end),mean(frz_T62(1:end,:),2,'omitnan'),'color',[0.72 0.82 0.92],'linewidth',0.1); % SYI bulk IMB
p = plot(t_fy,T_fy_bulk,'o','LineWidth',0.5,'color',c{1}); p.MarkerSize = 1.5; set(p,'markerfacecolor',get(p,'color'));
p = plot(t_sy,T_sy_bulk,'o','LineWidth',0.5,'color',c{2}); p.MarkerSize = 1.5; set(p,'markerfacecolor',get(p,'color'));
leg = legend('FYI','SYI','FYI (IMB)','SYI (IMB)','box','off','NumColumns',2); set(leg,'FontSize',7,'Location','north'); leg.ItemTokenSize = [30*0.33,18*0.33];
hYLabel = ylabel('Temperature (°C)'); set(gca,'FontSize',9,'FontWeight','normal'); set(hYLabel,'FontSize',8,'FontWeight','normal'); ylim([-12 1])
ax = gca; ax.XTick = datetime(['01-Oct-2019';'01-Dec-2019';'01-Feb-2020';'01-Apr-2020';'01-Jun-2020';'01-Aug-2020']); datetick('x','mmm','keepticks'); xtickangle(0); % time

nexttile % ice thickness
CT = T_fy_b; SA(1:20) = 33.7; SA(19:23) = Sb_fy_ct(3:7); rho_sw = gsw_rho(SA,CT,0); % seawater density
a_mp_cor = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 15 18 16 20]/100;
hs_fy(22:23) = [5.2 4.4]/100;
fb_bal = h_fy_avg/100 - (hs_fy.*rho_s.*(1-a_mp_cor) + rho_si_pr_bulk_fy.*h_fy_avg/100)./rho_sw; % FYI freeboard assuming hydrostatic balance
p = plot(t_fy,h_fy_avg/100,'o-','LineWidth',0.5,'color',c{1}); p.MarkerSize = 1.5; set(p,'markerfacecolor',get(p,'color')); hold on
p = plot(t_sy,h_sy_avg/100,'o-','LineWidth',0.5,'color',c{2}); p.MarkerSize = 1.5; set(p,'markerfacecolor',get(p,'color')); hold on
leg = legend('FYI','SYI','box','off','NumColumns',1); set(leg,'FontSize',7,'Location','best'); leg.ItemTokenSize = [30*0.33,18*0.33];
hYLabel = ylabel('Thickness (m)'); set(gca,'FontSize',9,'FontWeight','normal'); set(hYLabel,'FontSize',8,'FontWeight','normal'); ylim([0 3]); set(gca,'YDir','reverse'); % reverse y-axis
ax = gca; ax.XTick = datetime(['01-Oct-2019';'01-Dec-2019';'01-Feb-2020';'01-Apr-2020';'01-Jun-2020';'01-Aug-2020']); datetick('x','mmm','keepticks'); xtickangle(0); % time

nexttile % brine volume
p = plot(t_fy(~isnan(vb_bulk_fy)),vb_bulk_fy(~isnan(vb_bulk_fy))/10,'o','color',c{1}); p.MarkerSize = 1.5; set(p,'markerfacecolor',get(p,'color')); hold on;
p = plot(t_sy(~isnan(vb_bulk_sy)),vb_bulk_sy(~isnan(vb_bulk_sy))/10,'o','color',c{2}); p.MarkerSize = 1.5; set(p,'markerfacecolor',get(p,'color'));
p = plot(t_T66,smooth(vb_bulk_fy_int/10,50,'sgolay'),'-','color',[0.72 0.82 0.92]); set(p,'LineWidth',0.5);
p = plot(t_T62(1:end-47),smooth(vb_bulk_sy_int/10,50,'sgolay'),'b-','color',[0.92 0.82 0.82]); set(p,'LineWidth',0.5);
p = plot(t_fy(~isnan(vb_bulk_fy)),vb_bulk_fy(~isnan(vb_bulk_fy))/10,'o','color',c{1}); p.MarkerSize = 1.5; set(p,'markerfacecolor',get(p,'color'));
p = plot(t_sy(~isnan(vb_bulk_sy)),vb_bulk_sy(~isnan(vb_bulk_sy))/10,'o','color',c{2}); p.MarkerSize = 1.5; set(p,'markerfacecolor',get(p,'color'));
leg = legend('FYI','SYI','FYI, IMB temp.','SYI, IMB temp.','box','off','NumColumns',2); set(leg,'FontSize',7,'Location','best'); leg.ItemTokenSize = [30*0.5,18*0.5];
hYLabel = ylabel('Brine volume (%)'); set(gca,'FontSize',9,'FontWeight','normal'); set(hYLabel,'FontSize',8,'FontWeight','normal'); 
ax = gca; ax.XTick = datetime(['01-Oct-2019';'01-Dec-2019';'01-Feb-2020';'01-Apr-2020';'01-Jun-2020';'01-Aug-2020']); datetick('x','mmm','keepticks'); xtickangle(0); % time
ylim([0 30]); yticks(0:5:30);

nexttile % density
t_cor_ridge = [datetime('10-Jul-2020 12:00') datetime('15-Jul-2020 12:00') datetime('27-Jul-2020 12:00')]; vg_ridge = [3.23 3.23 3.48]; rho_ridge = [896.6 894.1 889.2];
t_cor_bhj = [datetime('02-Jul-2020 12:00')]; vg_bhj = 51/10;
p = plot(t_fy([1:3 5:6 8:23]),rho_si_pr_bulk_fy([1:3 5:6 8:23]),'o-','Color',c{1}); p.MarkerSize = 1.5; set(p,'markerfacecolor',get(p,'color')); hold on
p = plot(t_sy([1:3 5:14]),rho_si_pr_bulk_sy([1:3 5:14]),'o-','Color',c{2}); p.MarkerSize = 1.5; set(p,'markerfacecolor',get(p,'color'));
p = plot(t_sy(15:18),rho_si_pr_bulk_sy(15:18),'o-','Color',c{2}); p.MarkerSize = 1.5; set(p,'markerfacecolor',get(p,'color'));
p = plot(t_cor_ridge,rho_ridge,'o-','Color',c{4}); p.MarkerSize = 1.5; set(p,'markerfacecolor',get(p,'color'));
p = plot([datetime('22-Apr-2020') datetime('05-May-2020')],[907.4 903.0],'o','Color',c{4}); p.MarkerSize = 1.5; set(p,'markerfacecolor',get(p,'color')); % RMRF
leg = legend('FYI','SYI','','ridge','box','off','NumColumns',3); set(leg,'FontSize',7,'Location','south'); leg.ItemTokenSize = [30*0.5,18*0.5];
hYLabel = ylabel('Density (kg / m^3)'); set(gca,'FontSize',9,'FontWeight','normal'); set(hYLabel,'FontSize',8,'FontWeight','normal'); ylim([878 921]); yticks(880:10:920);
ax = gca; ax.XTick = datetime(['01-Oct-2019';'01-Dec-2019';'01-Feb-2020';'01-Apr-2020';'01-Jun-2020';'01-Aug-2020']); datetick('x','mmm','keepticks'); xtickangle(0); % time

nexttile % air volume
p = plot(t_fy(~isnan(vg_pr_bulk_fy)),vg_pr_bulk_fy(~isnan(vg_pr_bulk_fy))/10,'o-'); p.MarkerSize = 1.5; set(p,'markerfacecolor',get(p,'color')); hold on
p = plot(t_sy([1:3 5:14]),vg_pr_bulk_sy([1:3 5:14])/10,'o-','Color',c{2}); p.MarkerSize = 1.5; set(p,'markerfacecolor',get(p,'color'));
p = plot(t_sy(15:18),vg_pr_bulk_sy(15:18)/10,'o-','Color',c{2}); p.MarkerSize = 1.5; set(p,'markerfacecolor',get(p,'color'));
p = plot(t_cor_ridge,vg_ridge,'o-','Color',c{4}); p.MarkerSize = 1.5; set(p,'markerfacecolor',get(p,'color')); % leg 4 ridges
p = plot([datetime('22-Apr-2020') datetime('05-May-2020')],[2.81 2.13],'o','Color',c{4}); p.MarkerSize = 1.5; set(p,'markerfacecolor',get(p,'color')); % RMRF, leg 3
leg = legend('FYI','SYI','','Ridge','','box','off','NumColumns',3); set(leg,'FontSize',7,'Location','north'); leg.ItemTokenSize = [30*0.5,18*0.5];
hYLabel = ylabel('Air volume (%)'); set(gca,'FontSize',9,'FontWeight','normal'); set(hYLabel,'FontSize',8,'FontWeight','normal'); ylim([0 8]); yticks(0:2:8);
ax = gca; ax.XTick = datetime(['01-Oct-2019';'01-Dec-2019';'01-Feb-2020';'01-Apr-2020';'01-Jun-2020';'01-Aug-2020']); datetick('x','mmm','keepticks'); xtickangle(0); % time

annotation('textbox',[0 .505 0.02 .51],'String','(a)','EdgeColor','none','HorizontalAlignment','center','FontSize',8);
annotation('textbox',[0.33 .505 0.325 .51],'String','(b)','EdgeColor','none','HorizontalAlignment','center','FontSize',8);
annotation('textbox',[0 .38 0.02 .38],'String','(c)','EdgeColor','none','HorizontalAlignment','center','FontSize',8);
annotation('textbox',[0.33 .38 0.325 .38],'String','(d)','EdgeColor','none','HorizontalAlignment','center','FontSize',8);
annotation('textbox',[0 .51/2 0.02 .51/2],'String','(e)','EdgeColor','none','HorizontalAlignment','center','FontSize',8);
annotation('textbox',[0.33 .51/2 0.325 .51/2],'String','(f)','EdgeColor','none','HorizontalAlignment','center','FontSize',8);
annotation('textbox',[0 .14 0.02 .14],'String','(g)','EdgeColor','none','HorizontalAlignment','center','FontSize',8);
annotation('textbox',[0.33 .14 0.325 .14],'String','(h)','EdgeColor','none','HorizontalAlignment','center','FontSize',8);

%% Figure 8: summer air, salinity, temperature profiles (6.2 x 3.5 in)
figure
tile = tiledlayout(1,3); tile.TileSpacing = 'compact'; tile.Padding = 'compact';
nexttile
i = 19; p = plot(vg_pr_fy{1,i}(1:end-2)/10,-drho_fy{i,1}(1:end-2)+drho_fy{i,1}(end-2)+1.625+bot_T66_cor(18)-bot_T66_cor(21),'o-','color',[0.8 0.8 0.8]); p.MarkerSize = 2; set(p,'markerfacecolor',get(p,'color')); hold on
i = 21; p = plot(vg_pr_fy{1,i}(1:end-2)/10,-drho_fy{i,1}(1:end-2)+drho_fy{i,1}(end-2)+1.625+bot_T66_cor(18)-bot_T66_cor(21),'o-','color',c{4}); p.MarkerSize = 2; set(p,'markerfacecolor',get(p,'color')); hold on
i = 22; p = plot(vg_pr_fy{1,i}(1:end-2)/10,-drho_fy{i,1}(1:end-2)+drho_fy{i,1}(end-2)+1.625+bot_T66_cor(18)-bot_T66_cor(22),'o-','color',c{5}); p.MarkerSize = 2; set(p,'markerfacecolor',get(p,'color'));
i = 23; p = plot(vg_pr_fy{1,i}(1:end)/10,-drho_fy{i,1}(1:end)+drho_fy{i,1}(end)+1.625+bot_T66_cor(18)-bot_T66_cor(23),'o-','color',c{6}); p.MarkerSize = 2; set(p,'markerfacecolor',get(p,'color'));
leg = legend('22 Jun','13 Jul','21 Jul','29 Jul','box','off'); set(leg,'FontSize',6,'Location','east'); leg.ItemTokenSize = [30*0.33,18*0.33];
hXLabel = xlabel('FYI air volume (%)'); hYLabel = ylabel('Initial ice thickness (m)'); set([hXLabel hYLabel gca],'FontSize',8,'FontWeight','normal'); set(gca,'YDir','reverse');
yticks(-1.8:0.2:1.6); ylim([-0.05 1.55]); xticks(0:2:16); xlim([-1 14]); xtickangle(0);

nexttile
i = 19; p = plot(rho_si_pr_fy{1,i}(1:end-2),-drho_fy{i,1}(1:end-2)+drho_fy{i,1}(end-2)+1.625+bot_T66_cor(18)-bot_T66_cor(21),'o-','color',[0.8 0.8 0.8]); p.MarkerSize = 2; set(p,'markerfacecolor',get(p,'color')); hold on
i = 21; p = plot(rho_si_pr_fy{1,i}(1:end-2),-drho_fy{i,1}(1:end-2)+drho_fy{i,1}(end-2)+1.625+bot_T66_cor(18)-bot_T66_cor(21),'o-','color',c{4}); p.MarkerSize = 2; set(p,'markerfacecolor',get(p,'color')); hold on
i = 22; p = plot(rho_si_pr_fy{1,i}(1:end-2),-drho_fy{i,1}(1:end-2)+drho_fy{i,1}(end-2)+1.625+bot_T66_cor(18)-bot_T66_cor(22),'o-','color',c{5}); p.MarkerSize = 2; set(p,'markerfacecolor',get(p,'color'));
i = 23; p = plot(rho_si_pr_fy{1,i}(1:end),-drho_fy{i,1}(1:end)+drho_fy{i,1}(end)+1.625+bot_T66_cor(18)-bot_T66_cor(23),'o-','color',c{6}); p.MarkerSize = 2; set(p,'markerfacecolor',get(p,'color'));
leg = legend('22 Jun','13 Jul','21 Jul','29 Jul','box','off'); set(leg,'FontSize',6,'Location','west'); leg.ItemTokenSize = [30*0.33,18*0.33];
hXLabel = xlabel('FYI density (kg m^-^3)'); hYLabel = ylabel('Initial ice thickness (m)'); set([hXLabel hYLabel gca],'FontSize',8,'FontWeight','normal'); set(gca,'YDir','reverse');
yticks(-1.8:0.2:1.6); ylim([-0.05 1.55]); xticks(800:20:940); xlim([810 940]); xtickangle(0);

% CryoSat-2 freeboard radar penetration, IMB data:
time = ncread('IMB_CMW_MOSAiC_density_premelt.nc','time'); t_0 = (datetime('1979-01-01 00:00:00')); t_IMB = t_0 + days(time); % Import of IMB data
ice_thick_imb = -ncread('IMB_CMW_MOSAiC_density_premelt.nc','ice_thick'); snow_thick_imb = -ncread('IMB_CMW_MOSAiC_density_premelt.nc','snow_thick');
hi_imb = mean(ice_thick_imb([1:6 9],:),1); % IMB ice thickness
hs_imb = mean(snow_thick_imb([1:6 9],:),1); % IMB snow thickness

% CS2 ice thickness estimate from Landy and Dawson (2022) doi:10.5285/d8c66670-57ad-44fc-8fef-942a46734ecb
time = ncread('uit_cryosat2_seaicethickness_nh_80km_v1p7.nc','Time');
t_0 = (datetime('01-Jan-0000 00:00')); t_alt = t_0 + days(time-1);
lat = ncread('uit_cryosat2_seaicethickness_nh_80km_v1p7.nc','Latitude');
lon = ncread('uit_cryosat2_seaicethickness_nh_80km_v1p7.nc','Longitude');
hi = ncread('uit_cryosat2_seaicethickness_nh_80km_v1p7.nc','Sea_Ice_Thickness');
range = 217+12:237; t_alt = t_alt(range); hi = hi(:,:,range);
for i = 1:length(t_alt)
    [~,idx(i)] = min(((lat(:)-lat_CO(i)).^2+(lon(:)-lon_CO(i)).^2)); [idx_2D(1,i),idx_2D(2,i)] = ind2sub(size(lat),idx(i)); % GPS index of MOSAiC CO2
    hi_alt(i) = nanmean(nanmean(-hi(idx_2D(1,i)-1:idx_2D(1,i)+1,idx_2D(2,i)-1:idx_2D(2,i)+1,i))); % Landy, 3x3 grid
end

rho_fyi_imb = interp1(datenum(t_fy),rho_si_pr_bulk_fy,datenum(t_IMB'),'linear','extrap'); % interpolation of FYI bulk density for IMB
t_sp = datetime(['03-May-2020';'22-Jun-2020';'06-Jul-2020';'20-Jul-2020';'27-Jul-2020']); % snow pits, June-July 2020
rho_sn_sp =     [         295         416.7         420.5         427.5          389.5];
rho_sn_imb = interp1(datenum(t_sp),rho_sn_sp,datenum(t_IMB'),'linear','extrap'); % interpolation of snow bulk density for IMB
rho_sw_imb = 1024;
d_imb = -(rho_fyi_imb.*hi_imb+rho_sn_imb.*hs_imb)./rho_sw_imb; fb_imb = d_imb+hi_imb; % raw observations
rho_sw_c1 = 1024; rho_fyi_c1 = 917; rho_sn_c1 = 300; pen = 0.1; fb_imb_c1 = fb_imb; hs_imb_c1 = hs_imb; % retrieval from IMB freeboard estimate
hi_imb_c1 = (rho_sw_c1.*fb_imb_c1 + rho_sn_c1.*hs_imb_c1 + pen*hs_imb.*rho_sw_c1)./(rho_sw_c1-rho_fyi_c1);
rho_sn = 295; hs = max(hs_imb); hs = 0.30; d_imb = -(rho_fyi_imb.*hi_imb+rho_sn.*hs)./rho_sw_imb; fb_imb_1 = d_imb+hi_imb; % constant snow load
rho_sw_c2 = 1024; rho_fyi_c2 = 917; rho_sn_c2 = 300; pen = 0.0; fb_imb_c2 = fb_imb_1; hs_imb_c2 = hs_imb*30/16; % retrieval for constant snow load
hi_imb_c2 = (rho_sw_c2.*fb_imb_c2 + rho_sn_c2.*hs_imb_c2 + pen*hs_imb_c2.*rho_sw_c2)./(rho_sw_c2-rho_fyi_c2); % retrieval for constant snow load
pen = 0.27; hi_imb_c3 = (rho_sw_c2.*fb_imb_c2 + rho_sn_c2.*hs_imb_c2 + pen*hs_imb_c2.*rho_sw_c2)./(rho_sw_c2-rho_fyi_c2); % retrieval for constant snow load, small radar penetration

nexttile
p = plot(t_alt,hi_alt-min(hi_alt),'o','color','k'); p.MarkerSize = 2.5; set(p,'markerfacecolor',get(p,'color')); hold on % Landy
plot(t_IMB,max(hi_imb_c3)-hi_imb_c3,':','color',c{4},'linewidth',2.5); % IMB fb + constant snow load, 70 % penetration
plot(t_IMB,max(hi_imb_c2)-hi_imb_c2,'color',c{3}); % IMB fb + constant snow load, 100 % penetration
plot(t_IMB,max(hi_imb)-hi_imb,'color',c{2}); hold on % observations
p = plot(t_alt,hi_alt-min(hi_alt),'o','color','k'); p.MarkerSize = 2.5; set(p,'markerfacecolor',get(p,'color')); hold on % Landy
leg = legend('CryoSat-2','IMB, 70 % pen., const. snow load','IMB, const. snow load','IMB, observations','box','off','NumColumns',1); set(leg,'FontSize',7,'Location','best'); leg.ItemTokenSize = [30*0.5,18*0.5];
hYLabel = ylabel('Accumulated ice melt (m)'); set([hYLabel gca],'FontSize',8,'FontWeight','normal'); set(gca,'FontSize',7,'FontWeight','normal'); ylim([-0.1 3.5]);
t_start = datetime('01-May-2020'); t_end = datetime('01-Aug-2020'); xlim([t_start t_end]); datetick('x','mmm','keepticks'); xtickangle(0);
set(gca,'FontSize',8,'FontWeight','normal');

% nexttile
% p = plot(t_alt,-min(hi_alt)+hi_alt,'o','color',c{1}); p.MarkerSize = 2.5; set(p,'markerfacecolor',get(p,'color')); hold on % Landy
% plot(t_IMB,-hi_imb-9*hs_imb-min(-hi_imb-9*hs_imb),'-','color',c{2},'LineWidth',1.0); % Total melt, IMB
% plot(t_IMB,-hi_imb+max(hi_imb),'-','color',c{3},'LineWidth',1.0); % Ice melt, IMB
% p = plot(t_alt,-min(hi_alt)+hi_alt,'o','color',c{1}); p.MarkerSize = 2.5; set(p,'markerfacecolor',get(p,'color')); % Landy
% leg = legend('CryoSat-2','IMB, 0 % snow penetration','IMB, 100 % snow penetration','box','off','NumColumns',1); set(leg,'FontSize',7,'Location','best'); leg.ItemTokenSize = [30*0.66,18*0.66];
% hYLabel = ylabel('Accumulated ice melt (m)'); set([hYLabel gca],'FontSize',8,'FontWeight','normal'); set(gca,'FontSize',7,'FontWeight','normal'); ylim([-0.1 3.0]);
% t_start = datetime('01-May-2020'); t_end = datetime('01-Aug-2020'); xlim([t_start t_end]); datetick('x','mmm','keepticks'); xtickangle(0);
% set(gca,'FontSize',8,'FontWeight','normal');

annotation('textbox',[0.01 .51 0.01 .51],'String','(a)','FontSize',8,'EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.22 .51 0.22 .51],'String','(b)','FontSize',8,'EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.43 .51 0.44 .51],'String','(c)','FontSize',8,'EdgeColor','none','HorizontalAlignment','center');

%% Figure 10: geochemistry of sea-ice
vg_fy_cold = vertcat(vg_pr_fy{1:23})/10;
vb_fy_cold = vertcat(vb_rho_fy{1:23});
d_fy_cold = vertcat(drho_fy{1:23});
perm = repmat({'Above water, perm.'},length(vg_fy_cold),1);
perm(vb_fy_cold >= 50 & d_fy_cold >  -0.03) = {'Above water, perm.'};
perm(vb_fy_cold <  50 & d_fy_cold >  -0.03) = {'Above water, imperm.'};
perm(vb_fy_cold >= 50 & d_fy_cold <= -0.03) = {'Underwater, perm.'};
perm(vb_fy_cold <  50 & d_fy_cold <= -0.03) = {'Underwater, imperm.'};
season = repmat({'winter'},length(vg_fy_cold),1);
season(349:end) = {'summer'};

figure
tile = tiledlayout(2,4); tile.TileSpacing = 'compact'; tile.Padding = 'compact';
nexttile
boxchart(categorical(perm),vg_fy_cold,'orientation','horizontal','MarkerStyle','none','GroupByColor',season,'LineWidth',0.5); colororder([c{2}; c{1}]);
leg = legend('box','off','NumColumns',1); set(leg,'FontSize',7,'Location','northeast'); leg.ItemTokenSize = [30*0.33,18*0.33];
set(gca(),'YTickLabel',{sprintf('Above water\\newlineimpermeable') sprintf('Above water\\newlinepermeable') sprintf('Underwater\\newlineimpermeable') sprintf('Underwater\\newlinepermeable') });
hXLabel = xlabel('FYI air volume (%)'); set([hXLabel gca],'FontSize',7,'FontWeight','normal'); xlim([-1 12]); xticks(-0:2:12); xtickangle(0);

nexttile % FYI vg vs vb
msz = 1.6;
i = 1; h = plot(0.1*vb_rho_fy{i}(drho_fy{i}>=-0.03),0.1*vg_pr_fy{i}(drho_fy{i}>=-0.03),'o','LineWidth',0.05,'color',c{1}); h.MarkerSize = msz*1.0; set(h,'markerfacecolor',get(h,'color')); hold on
i = 1; h = plot(0.1*vb_rho_fy{i}(drho_fy{i}<-0.03),0.1*vg_pr_fy{i}(drho_fy{i}<-0.03),'o','LineWidth',0.05,'color',c{1}); h.MarkerSize = msz; set(h,'markerfacecolor','w'); hold on
i = 19; h = plot(0.1*vb_rho_fy{i}(drho_fy{i}>=-0.03),0.1*vg_pr_fy{i}(drho_fy{i}>=-0.03),'o','LineWidth',0.05,'color',c{2}); h.MarkerSize = msz; set(h,'markerfacecolor',get(h,'color')); hold on
i = 19; h = plot(0.1*vb_rho_fy{i}(drho_fy{i}<-0.03),0.1*vg_pr_fy{i}(drho_fy{i}<-0.03),'o','LineWidth',0.05,'color',c{2}); h.MarkerSize = msz; set(h,'markerfacecolor','w'); hold on
for i = [1:3 5:6 8:18] % cold season, below water
    h = plot(0.1*vb_rho_fy{i}(drho_fy{i}<-0.03),0.1*vg_pr_fy{i}(drho_fy{i}<-0.03),'o','LineWidth',0.05,'color',c{1}); h.MarkerSize = msz; set(h,'markerfacecolor','w'); hold on
end
for i = 19:23 % melt season, below water
    h = plot(0.1*vb_rho_fy{i}(drho_fy{i}<-0.03),0.1*vg_pr_fy{i}(drho_fy{i}<-0.03),'o','LineWidth',0.05,'color',c{2}); h.MarkerSize = msz; set(h,'markerfacecolor','w'); hold on
end
for i = [1:2 5:6 8:17] % cold season, above water
    h = plot(0.1*vb_rho_fy{i}(drho_fy{i}>=-0.03),0.1*vg_pr_fy{i}(drho_fy{i}>=-0.03),'o','LineWidth',0.05,'color',c{1}); h.MarkerSize = msz*1.0; set(h,'markerfacecolor',get(h,'color')); hold on
end
for i = 19:23 % melt season, above water
    h = plot(0.1*vb_rho_fy{i}(drho_fy{i}>=-0.03),0.1*vg_pr_fy{i}(drho_fy{i}>=-0.03),'o','LineWidth',0.05,'color',c{2}); h.MarkerSize = msz*1.0; set(h,'markerfacecolor',get(h,'color')); hold on
end
plot([0 50],[0 5],'--','LineWidth',0.05,'color','k');
p = text(38,3.4,'gas = 10% brine'); set(p,'Color','k','HorizontalAlignment','center','FontSize',7); set(p,'Rotation',28);
hXLabel = xlabel('FYI brine volume (%)'); hYLabel = ylabel('FYI air volume (%)'); set([hXLabel hYLabel gca],'FontSize',7,'FontWeight','normal'); xlim([0 50]); ylim([0 14]); xtickangle(0);
leg = legend('winter, above water','winter, below water','summer, above water','summer, below water','box','on','NumColumns',1); set(leg,'FontSize',6,'Location','northeast');
leg.ItemTokenSize = [30*0.2,18*0.2]; leg.BoxFace.ColorType='truecoloralpha'; leg.BoxFace.ColorData=uint8(255*[1 1 1 0.7]'); % transparent legend

nexttile % vg vs T
i = 1; h = plot(T_rho_fy{i}(drho_fy{i}>=-0.03),0.1*vg_pr_fy{i}(drho_fy{i}>=-0.03),'o','LineWidth',0.05,'color',c{1}); h.MarkerSize = msz*1.0; set(h,'markerfacecolor',get(h,'color')); hold on
i = 1; h = plot(T_rho_fy{i}(drho_fy{i}<-0.03),0.1*vg_pr_fy{i}(drho_fy{i}<-0.03),'o','LineWidth',0.05,'color',c{1}); h.MarkerSize = msz; set(h,'markerfacecolor','w'); hold on
i = 19; h = plot(T_rho_fy{i}(drho_fy{i}>=-0.03),0.1*vg_pr_fy{i}(drho_fy{i}>=-0.03),'o','LineWidth',0.05,'color',c{2}); h.MarkerSize = msz*1.0; set(h,'markerfacecolor',get(h,'color')); hold on
i = 19; h = plot(T_rho_fy{i}(drho_fy{i}<-0.03),0.1*vg_pr_fy{i}(drho_fy{i}<-0.03),'o','LineWidth',0.05,'color',c{2}); h.MarkerSize = msz; set(h,'markerfacecolor','w'); hold on
for i = [1:3 5:6 8:18] % cold season, below water
    h = plot(T_rho_fy{i}(drho_fy{i}<-0.03),0.1*vg_pr_fy{i}(drho_fy{i}<-0.03),'o','LineWidth',0.05,'color',c{1}); h.MarkerSize = msz; set(h,'markerfacecolor','w'); hold on
end
for i = 19:23 % melt season, below water
    h = plot(T_rho_fy{i}(drho_fy{i}<-0.03),0.1*vg_pr_fy{i}(drho_fy{i}<-0.03),'o','LineWidth',0.05,'color',c{2}); h.MarkerSize = msz; set(h,'markerfacecolor','w'); hold on
end
for i = [1:2 5:6 8:17] % cold season, above water
    h = plot(T_rho_fy{i}(drho_fy{i}>=-0.03),0.1*vg_pr_fy{i}(drho_fy{i}>=-0.03),'o','LineWidth',0.05,'color',c{1}); h.MarkerSize = msz*1.0; set(h,'markerfacecolor',get(h,'color')); hold on
end
for i = 19:23 % melt season, above water
    h = plot(T_rho_fy{i}(drho_fy{i}>=-0.03),0.1*vg_pr_fy{i}(drho_fy{i}>=-0.03),'o','LineWidth',0.05,'color',c{2}); h.MarkerSize = msz*1.0; set(h,'markerfacecolor',get(h,'color')); hold on
end
hXLabel = xlabel('FYI temperature (°C)'); hYLabel = ylabel('FYI air volume (%)'); set([hXLabel hYLabel gca],'FontSize',7,'FontWeight','normal'); xlim([-20 0]); ylim([0 14]); xtickangle(0);
% set(gca,'XScale','log'); xtickformat('$%,.0f')
leg = legend('winter, above water','winter, below water','summer, above water','summer, below water','box','on','NumColumns',1); set(leg,'FontSize',6,'Location','northeast');
leg.ItemTokenSize = [30*0.2,18*0.2]; leg.BoxFace.ColorType='truecoloralpha'; leg.BoxFace.ColorData=uint8(255*[1 1 1 0.7]'); % transparent legend

nexttile % vg vs T
i = 1; h = plot(T_rho_fy{i}(drho_fy{i}>=-0.03),0.1*vg_pr_fy{i}(drho_fy{i}>=-0.03),'o','LineWidth',0.05,'color',c{1}); h.MarkerSize = msz*1.0; set(h,'markerfacecolor',get(h,'color')); hold on
i = 1; h = plot(T_rho_fy{i}(drho_fy{i}<-0.03),0.1*vg_pr_fy{i}(drho_fy{i}<-0.03),'o','LineWidth',0.05,'color',c{1}); h.MarkerSize = msz; set(h,'markerfacecolor','w'); hold on
i = 19; h = plot(T_rho_fy{i}(drho_fy{i}>=-0.03),0.1*vg_pr_fy{i}(drho_fy{i}>=-0.03),'o','LineWidth',0.05,'color',c{2}); h.MarkerSize = msz*1.0; set(h,'markerfacecolor',get(h,'color')); hold on
i = 19; h = plot(T_rho_fy{i}(drho_fy{i}<-0.03),0.1*vg_pr_fy{i}(drho_fy{i}<-0.03),'o','LineWidth',0.05,'color',c{2}); h.MarkerSize = msz; set(h,'markerfacecolor','w'); hold on
for i = [1:3 5:6 8:18] % cold season, below water
    h = plot(T_rho_fy{i}(drho_fy{i}<-0.03),0.1*vg_pr_fy{i}(drho_fy{i}<-0.03),'o','LineWidth',0.05,'color',c{1}); h.MarkerSize = msz; set(h,'markerfacecolor','w'); hold on
end
for i = 19:23 % melt season, below water
    h = plot(T_rho_fy{i}(drho_fy{i}<-0.03),0.1*vg_pr_fy{i}(drho_fy{i}<-0.03),'o','LineWidth',0.05,'color',c{2}); h.MarkerSize = msz; set(h,'markerfacecolor','w'); hold on
end
for i = [1:2 5:6 8:17] % cold season, above water
    h = plot(T_rho_fy{i}(drho_fy{i}>=-0.03),0.1*vg_pr_fy{i}(drho_fy{i}>=-0.03),'o','LineWidth',0.05,'color',c{1}); h.MarkerSize = msz*1.0; set(h,'markerfacecolor',get(h,'color')); hold on
end
for i = 19:23 % melt season, above water
    h = plot(T_rho_fy{i}(drho_fy{i}>=-0.03),0.1*vg_pr_fy{i}(drho_fy{i}>=-0.03),'o','LineWidth',0.05,'color',c{2}); h.MarkerSize = msz*1.0; set(h,'markerfacecolor',get(h,'color')); hold on
end
plot([-5 -1.127 -0.1],[1.2 1.2 7.7413],'-','LineWidth',2,'color','k');
p = text(-2.85,8.5,'y = max (6.4x + 8.4, 1.2)'); set(p,'Color','k','HorizontalAlignment','left','FontSize',8,'FontWeight','bold');
hXLabel = xlabel('FYI temperature (°C)'); hYLabel = ylabel('FYI air volume (%)'); set([hXLabel hYLabel gca],'FontSize',7,'FontWeight','normal'); xlim([-3 0]); ylim([0 14]);
leg = legend('winter, above water','winter, below water','summer, above water','summer, below water','box','on','NumColumns',1); set(leg,'FontSize',6,'Location','northeast');
leg.ItemTokenSize = [30*0.2,18*0.2]; leg.BoxFace.ColorType='truecoloralpha'; leg.BoxFace.ColorData=uint8(255*[1 1 1 0.7]'); % transparent legend

vg_sy_cold = vertcat(vg_pr_sy{1:18})/10;
vb_sy_cold = vertcat(vb_rho_sy{1:18});
d_sy_cold = vertcat(drho_sy{1:18});
perm = repmat({'Above water, perm.'},length(vg_sy_cold),1);
perm(vb_sy_cold >= 50 & d_sy_cold >  -0.03) = {'Above water, perm.'};
perm(vb_sy_cold <  50 & d_sy_cold >  -0.03) = {'Above water, imperm.'};
perm(vb_sy_cold >= 50 & d_sy_cold <= -0.03) = {'Underwater, perm.'};
perm(vb_sy_cold <  50 & d_sy_cold <= -0.03) = {'Underwater, imperm.'};
season = repmat({'winter'},length(vg_sy_cold),1);
season(266:end) = {'summer'};

nexttile
boxchart(categorical(perm),vg_sy_cold,'orientation','horizontal','MarkerStyle','none','GroupByColor',season,'LineWidth',0.05); colororder([c{2}; c{1}]);
leg = legend('box','off','NumColumns',1); set(leg,'FontSize',7,'Location','northeast'); leg.ItemTokenSize = [30*0.33,18*0.33];
set(gca(),'YTickLabel',{sprintf('Above water\\newlineimpermeable') sprintf('Above water\\newlinepermeable') sprintf('Underwater\\newlineimpermeable') sprintf('Underwater\\newlinepermeable') });
hXLabel = xlabel('SYI air volume (%)'); set([hXLabel gca],'FontSize',7,'FontWeight','normal'); xlim([-1 12]); xticks(-0:2:12); xtickangle(0);

nexttile % SYI vg vs vb
msz = 1.6;
i = 1; h = plot(0.1*vb_rho_sy{i}(drho_sy{i}>=-0.03),0.1*vg_pr_sy{i}(drho_sy{i}>=-0.03),'o','LineWidth',0.05,'color',c{1}); h.MarkerSize = msz*1.0; set(h,'markerfacecolor',get(h,'color')); hold on
i = 1; h = plot(0.1*vb_rho_sy{i}(drho_sy{i}<-0.03),0.1*vg_pr_sy{i}(drho_sy{i}<-0.03),'o','LineWidth',0.05,'color',c{1}); h.MarkerSize = msz; set(h,'markerfacecolor','w'); hold on
i = 15; h = plot(0.1*vb_rho_sy{i}(drho_sy{i}>=-0.03),0.1*vg_pr_sy{i}(drho_sy{i}>=-0.03),'o','LineWidth',0.05,'color',c{2}); h.MarkerSize = msz; set(h,'markerfacecolor',get(h,'color')); hold on
i = 15; h = plot(0.1*vb_rho_sy{i}(drho_sy{i}<-0.03),0.1*vg_pr_sy{i}(drho_sy{i}<-0.03),'o','LineWidth',0.05,'color',c{2}); h.MarkerSize = msz; set(h,'markerfacecolor','w'); hold on
for i = [1:3 5:14] % cold season, below water
    h = plot(0.1*vb_rho_sy{i}(drho_sy{i}<-0.03),0.1*vg_pr_sy{i}(drho_sy{i}<-0.03),'o','LineWidth',0.05,'color',c{1}); h.MarkerSize = msz; set(h,'markerfacecolor','w'); hold on
end
for i = 15:18 % melt season, below water
    h = plot(0.1*vb_rho_sy{i}(drho_sy{i}<-0.03),0.1*vg_pr_sy{i}(drho_sy{i}<-0.03),'o','LineWidth',0.05,'color',c{2}); h.MarkerSize = msz; set(h,'markerfacecolor','w'); hold on
end
for i = [1:3 5:14] % cold season, above water
    h = plot(0.1*vb_rho_sy{i}(drho_sy{i}>=-0.03),0.1*vg_pr_sy{i}(drho_sy{i}>=-0.03),'o','LineWidth',0.05,'color',c{1}); h.MarkerSize = msz*1.0; set(h,'markerfacecolor',get(h,'color')); hold on
end
for i = 15:18 % melt season, above water
    h = plot(0.1*vb_rho_sy{i}(drho_sy{i}>=-0.03),0.1*vg_pr_sy{i}(drho_sy{i}>=-0.03),'o','LineWidth',0.05,'color',c{2}); h.MarkerSize = msz*1.0; set(h,'markerfacecolor',get(h,'color')); hold on
end
plot([0 50],[0 5],'--','LineWidth',0.05,'color','k');
% p = text(10,2.0,'gas = 10% brine'); set(p,'Color','k','HorizontalAlignment','center','FontSize',7); set(p,'Rotation',31);
hXLabel = xlabel('SYI brine volume (%)'); hYLabel = ylabel('SYI air volume (%)'); set([hXLabel hYLabel gca],'FontSize',7,'FontWeight','normal'); xlim([0 15]); ylim([0 14]); xtickangle(0);
leg = legend('winter, above water','winter, below water','summer, above water','summer, below water','box','on','NumColumns',1); set(leg,'FontSize',6,'Location','northeast');
leg.ItemTokenSize = [30*0.2,18*0.2]; leg.BoxFace.ColorType='truecoloralpha'; leg.BoxFace.ColorData=uint8(255*[1 1 1 0.7]'); % transparent legend

nexttile % vg vs T
i = 1; h = plot(T_rho_sy{i}(drho_sy{i}>=-0.03),0.1*vg_pr_sy{i}(drho_sy{i}>=-0.03),'o','LineWidth',0.05,'color',c{1}); h.MarkerSize = msz*1.0; set(h,'markerfacecolor',get(h,'color')); hold on
i = 1; h = plot(T_rho_sy{i}(drho_sy{i}<-0.03),0.1*vg_pr_sy{i}(drho_sy{i}<-0.03),'o','LineWidth',0.05,'color',c{1}); h.MarkerSize = msz; set(h,'markerfacecolor','w'); hold on
i = 15; h = plot(T_rho_sy{i}(drho_sy{i}>=-0.03),0.1*vg_pr_sy{i}(drho_sy{i}>=-0.03),'o','LineWidth',0.05,'color',c{2}); h.MarkerSize = msz*1.0; set(h,'markerfacecolor',get(h,'color')); hold on
i = 15; h = plot(T_rho_sy{i}(drho_sy{i}<-0.03),0.1*vg_pr_sy{i}(drho_sy{i}<-0.03),'o','LineWidth',0.05,'color',c{2}); h.MarkerSize = msz; set(h,'markerfacecolor','w'); hold on
for i = [1:3 5:14] % cold season, below water
    h = plot(T_rho_sy{i}(drho_sy{i}<-0.03),0.1*vg_pr_sy{i}(drho_sy{i}<-0.03),'o','LineWidth',0.05,'color',c{1}); h.MarkerSize = msz; set(h,'markerfacecolor','w'); hold on
end
for i = 15:18 % melt season, below water
    h = plot(T_rho_sy{i}(drho_sy{i}<-0.03),0.1*vg_pr_sy{i}(drho_sy{i}<-0.03),'o','LineWidth',0.05,'color',c{2}); h.MarkerSize = msz; set(h,'markerfacecolor','w'); hold on
end
for i = [1:3 5:14] % cold season, above water
    h = plot(T_rho_sy{i}(drho_sy{i}>=-0.03),0.1*vg_pr_sy{i}(drho_sy{i}>=-0.03),'o','LineWidth',0.05,'color',c{1}); h.MarkerSize = msz*1.0; set(h,'markerfacecolor',get(h,'color')); hold on
end
for i = 15:18 % melt season, above water
    h = plot(T_rho_sy{i}(drho_sy{i}>=-0.03),0.1*vg_pr_sy{i}(drho_sy{i}>=-0.03),'o','LineWidth',0.05,'color',c{2}); h.MarkerSize = msz*1.0; set(h,'markerfacecolor',get(h,'color')); hold on
end
hXLabel = xlabel('SYI temperature (°C)'); hYLabel = ylabel('SYI air volume (%)'); set([hXLabel hYLabel gca],'FontSize',7,'FontWeight','normal'); xlim([-20 0]); ylim([0 14]); xtickangle(0);
leg = legend('winter, above water','winter, below water','summer, above water','summer, below water','box','on','NumColumns',1); set(leg,'FontSize',6,'Location','northeast');
leg.ItemTokenSize = [30*0.2,18*0.2]; leg.BoxFace.ColorType='truecoloralpha'; leg.BoxFace.ColorData=uint8(255*[1 1 1 0.7]'); % transparent legend

nexttile % vg vs T
i = 1; h = plot(T_rho_sy{i}(drho_sy{i}>=-0.03),0.1*vg_pr_sy{i}(drho_sy{i}>=-0.03),'o','LineWidth',0.05,'color',c{1}); h.MarkerSize = msz*1.0; set(h,'markerfacecolor',get(h,'color')); hold on
i = 1; h = plot(T_rho_sy{i}(drho_sy{i}<-0.03),0.1*vg_pr_sy{i}(drho_sy{i}<-0.03),'o','LineWidth',0.05,'color',c{1}); h.MarkerSize = msz; set(h,'markerfacecolor','w'); hold on
i = 15; h = plot(T_rho_sy{i}(drho_sy{i}>=-0.03),0.1*vg_pr_sy{i}(drho_sy{i}>=-0.03),'o','LineWidth',0.05,'color',c{2}); h.MarkerSize = msz*1.0; set(h,'markerfacecolor',get(h,'color')); hold on
i = 15; h = plot(T_rho_sy{i}(drho_sy{i}<-0.03),0.1*vg_pr_sy{i}(drho_sy{i}<-0.03),'o','LineWidth',0.05,'color',c{2}); h.MarkerSize = msz; set(h,'markerfacecolor','w'); hold on
for i = [1:3 5:14] % cold season, below water
    h = plot(T_rho_sy{i}(drho_sy{i}<-0.03),0.1*vg_pr_sy{i}(drho_sy{i}<-0.03),'o','LineWidth',0.05,'color',c{1}); h.MarkerSize = msz; set(h,'markerfacecolor','w'); hold on
end
for i = 15:18 % melt season, below water
    h = plot(T_rho_sy{i}(drho_sy{i}<-0.03 & T_rho_sy{i}<-0.15),0.1*vg_pr_sy{i}(drho_sy{i}<-0.03 & T_rho_sy{i}<-0.15),'o','LineWidth',0.05,'color',c{2}); h.MarkerSize = msz; set(h,'markerfacecolor','w'); hold on
end
for i = [1:3 5:14] % cold season, above water
    h = plot(T_rho_sy{i}(drho_sy{i}>=-0.03),0.1*vg_pr_sy{i}(drho_sy{i}>=-0.03),'o','LineWidth',0.05,'color',c{1}); h.MarkerSize = msz*1.0; set(h,'markerfacecolor',get(h,'color')); hold on
end
for i = 15:18 % melt season, above water
    h = plot(T_rho_sy{i}(drho_sy{i}>=-0.03),0.1*vg_pr_sy{i}(drho_sy{i}>=-0.03),'o','LineWidth',0.05,'color',c{2}); h.MarkerSize = msz*1.0; set(h,'markerfacecolor',get(h,'color')); hold on
end
plot([-5 -1.6 -0.1],[1.02 1.02 3.30],'-','LineWidth',2,'color','k');
p = text(-2.85,4,'y = max (1.7x + 3.5, 1.0)'); set(p,'Color','k','HorizontalAlignment','left','FontSize',8,'FontWeight','bold');
hXLabel = xlabel('SYI temperature (°C)'); hYLabel = ylabel('SYI air volume (%)'); set([hXLabel hYLabel gca],'FontSize',7,'FontWeight','normal'); xlim([-3 0]); ylim([0 14]);
leg = legend('winter, above water','winter, below water','summer, above water','summer, below water','box','on','NumColumns',1); set(leg,'FontSize',6,'Location','northeast');
leg.ItemTokenSize = [30*0.2,18*0.2]; leg.BoxFace.ColorType='truecoloralpha'; leg.BoxFace.ColorData=uint8(255*[1 1 1 0.7]'); % transparent legend

annotation('textbox',[0 .50 0.02 .51],'String','(a)','FontSize',8,'EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.185 .50 0.185 .51],'String','(b)','FontSize',8,'EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.34 .50 0.34 .51],'String','(c)','FontSize',8,'EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.495 .50 0.495 .51],'String','(d)','FontSize',8,'EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0 .51/2 0.02 .51/2],'String','(e)','FontSize',8,'EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.185 .51/2 0.185 .51/2],'String','(f)','FontSize',8,'EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.34 .51/2 0.34 .51/2],'String','(g)','FontSize',8,'EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.495 .51/2 0.495 .51/2],'String','(i)','FontSize',8,'EdgeColor','none','HorizontalAlignment','center');

%% Figure 6: ALS freeboard for FYI ROV site: 10 May, 30 June, 17 July, 22 July 
clear; clc;
tile = tiledlayout(2,4); tile.TileSpacing = 'compact'; tile.Padding = 'compact';

nexttile % 10 May
project = 'awi-mosaic-l4-als-vq580-stere-0p50m-20200510T135817-20200510T154321-fv2p0.nc'; t = datetime('10-May-2020'); % from Hutter et al., 10.1594/PANGAEA.950896
freeboard = ncread(project,'freeboard'); xc = ncread(project,'xc'); yc = ncread(project,'yc');
xc = -xc - 430.0; yc = -yc - 600 -1000; % corrections for May 10
[~,idx_x_min] = min(abs(-xc--300)); [~,idx_x_max] = min(abs(-xc-0)); % ROV
[~,idx_y_min] = min(abs(-yc--100)); [~,idx_y_max] = min(abs(-yc-80)); % ROV
xv = -xc(idx_x_min:idx_x_max); % selected area X
yv = -yc(idx_y_min:idx_y_max)'; % selected area Y
fb = freeboard(idx_x_min:idx_x_max,idx_y_min:idx_y_max)'; % selected area Z

load('density_data.mat',"als_fyi","als_syi","area_6");
theta = 41; R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
for i = 1:length(als_fyi); FYI_area(i,:) = R*als_fyi(i,:)'; end % rotation
FYI_area(:,1) = FYI_area(:,1)+180; FYI_area(:,2) = FYI_area(:,2)+105; % motion
for i = 1:length(als_syi); SYI_area(i,:) = R*als_syi(i,:)'; end % rotation
SYI_area(:,1) = SYI_area(:,1)+180; SYI_area(:,2) = SYI_area(:,2)+105; % motion

theta = -12; R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)]; for j = 1:length(area_6); FYI_ROV(j,:) = R*area_6(j,:)'; end
FYI_ROV(:,1) = FYI_ROV(:,1)-460; FYI_ROV(:,2) = FYI_ROV(:,2)-29; clearvars j R % motion
theta = 41; R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
for i = 1:length(FYI_ROV); FYI_ROV(i,:) = R*FYI_ROV(i,:)'; end % rotation
FYI_ROV(:,1) = FYI_ROV(:,1)+180; FYI_ROV(:,2) = FYI_ROV(:,2)+105; % motion

[XX,YY] = meshgrid(xv, yv);
[in,~] = inpolygon(XX(:),YY(:),FYI_ROV(:,1),FYI_ROV(:,2)); inside = reshape(in,length(xv),length(yv));
fb_rov = freeboard(idx_x_min:idx_x_max,idx_y_min:idx_y_max)'; fb_rov(inside==0) = NaN;
% rotation
dx = -205; dy = 38;
XY = [XX(:) YY(:)];
theta=-41.5; R=[cosd(theta) -sind(theta); sind(theta) cosd(theta)]; rotXY=XY*R'; % rotation
Xqr = reshape(rotXY(:,1),size(XX,1),[]); Yqr = reshape(rotXY(:,2),size(YY,1),[]);
xv_rot = Xqr+dx; yv_rot = Yqr+dy;
for i = 1:length(FYI_ROV); FYI_ROV_rot(i,:) = R*FYI_ROV(i,:)'; end % rotation
FYI_ROV_rot(:,1) = FYI_ROV_rot(:,1)+dx; FYI_ROV_rot(:,2) = FYI_ROV_rot(:,2)+dy; % motion

range = -0.1:0.025:0.4;
contourf(xv_rot,yv_rot,fb_rov,range,'LabelSpacing',10,'edgecolor','none'); hold on % Rotated & Shifted
x_new{1} = xv_rot; y_new{1} = yv_rot; fb_new{1} = fb_rov;
load('density_data.mat',"batlow"); colormap(batlow);
p = text(-405,50,sprintf('fb = %.2f m',nanmean(fb_rov(:)) )); set(p,'Color','k','HorizontalAlignment','left','FontSize',7);
p = text(-405,60,'sn = 0.16 m'); set(p,'Color','k','HorizontalAlignment','left','FontSize',7);
hXLabel = xlabel('x (m)'); hYLabel = ylabel('y (m)'); set([hXLabel hYLabel gca],'FontSize',7,'FontWeight','normal'); xtickangle(0);
xticks(-410:40:-200); xlim([-410 -210]); yticks(40:20:220); ylim([40 220]);
title(sprintf(['FYI freeboard, ', datestr(t,'dd mmmm')]),'FontSize',8,'FontWeight','normal');

nexttile % 30 June
project = "awi-mosaic-l4-als-vq580-stere-0p50m-20200630T074032-20200630T090829-fv2p0.nc"; t = datetime('30-Jun-2020'); % from Hutter et al., 10.1594/PANGAEA.950896
freeboard = ncread(project,'freeboard'); xc = ncread(project,'xc'); yc = ncread(project,'yc');
xc = xc - 0.0; yc = yc - 9.5; % corrections for June 30 scan
   
load('density_data.mat',"area_6");
theta = -12; R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)]; for i = 1:length(area_6); FYI_area(i,:) = R*area_6(i,:)'; end % rotation
FYI_area(:,1) = FYI_area(:,1)-460; FYI_area(:,2) = FYI_area(:,2)-29; % motion
[~,idx_x_min] = min(abs(-xc--200)); [~,idx_x_max] = min(abs(-xc--450)); % x from -500 to -100
[~,idx_y_min] = min(abs(-yc-250)); [~,idx_y_max] = min(abs(-yc-30)); % y from 0 to 250
xv = -xc(idx_x_min:idx_x_max); % selected area X
yv = -yc(idx_y_min:idx_y_max)'; % selected area Y
[~,idx_x_min_tent] = min(abs(-xc--240)); [~,idx_x_max_tent] = min(abs(-xc--265)); % x from -500 to -100
[~,idx_y_min_tent] = min(abs(-yc-65)); [~,idx_y_max_tent] = min(abs(-yc-50)); % y from 0 to 250
freeboard_tent = freeboard(idx_x_min_tent:idx_x_max_tent,idx_y_min_tent:idx_y_max_tent); % ROV tent fix
freeboard_tent(freeboard_tent > 0.8) = 0.2;
freeboard(idx_x_min_tent:idx_x_max_tent,idx_y_min_tent:idx_y_max_tent) = freeboard_tent;
fb = freeboard(idx_x_min:idx_x_max,idx_y_min:idx_y_max)'; % selected area Z
[XX,YY] = meshgrid(xv, yv);
[in,~] = inpolygon(XX(:),YY(:),FYI_area(:,1),FYI_area(:,2)); inside = reshape(in,length(xv),length(yv));
fb_fyi = freeboard(idx_x_min:idx_x_max,idx_y_min:idx_y_max)'; fb_fyi(inside==0) = NaN;
fb_fyi_filled = fillmissing2(fb_fyi,"linear"); fb_fyi_filled(inside==0) = NaN;
range = -0.1:0.025:0.4; % freeboard range
contourf(xv,yv,fb_fyi_filled,range,'LabelSpacing',10,'edgecolor','none'); hold on;
x_new{2} = xv; y_new{2} = yv; fb_new{2} = fb_fyi;
colormap(batlow);
p = text(-405,50,sprintf('fb = %.2f m',nanmean(fb_fyi(:)) )); set(p,'Color','k','HorizontalAlignment','left','FontSize',7);
p = text(-405,60,'sn = 0.08 m'); set(p,'Color','k','HorizontalAlignment','left','FontSize',7);
hXLabel = xlabel('x (m)'); hYLabel = ylabel('y (m)'); set([hXLabel hYLabel gca],'FontSize',7,'FontWeight','normal'); xtickangle(0);
xticks(-410:40:-200); xlim([-410 -210]); yticks(40:20:220); ylim([40 220]);
title(sprintf(['FYI freeboard, ', datestr(t,'dd mmmm')]),'FontSize',8,'FontWeight','normal');

nexttile % 17 July
clearvars -except x_new y_new fb_new
project = 'awi-mosaic-l4-als-vq580-stere-0p50m-20200717T152302-20200717T171735-fv2p0.nc'; t = datetime('17-Jul-2020'); % from Hutter et al., 10.1594/PANGAEA.950896
freeboard = ncread(project,'freeboard'); xc = ncread(project,'xc'); yc = ncread(project,'yc');
    
load('density_data.mat',"area_6");
theta = -12; R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)]; for i = 1:length(area_6); FYI_area(i,:) = R*area_6(i,:)'; end % rotation
FYI_area(:,1) = FYI_area(:,1)-460; FYI_area(:,2) = FYI_area(:,2)-29; % motion
[~,idx_x_min] = min(abs(-xc--100)); [~,idx_x_max] = min(abs(-xc--500)); % x from -500 to -100
[~,idx_y_min] = min(abs(-yc-250)); [~,idx_y_max] = min(abs(-yc-0)); % y from 0 to 250
xv = -xc(idx_x_min:idx_x_max); % selected area X
yv = -yc(idx_y_min:idx_y_max)'; % selected area Y
[~,idx_x_min_tent] = min(abs(-xc--240)); [~,idx_x_max_tent] = min(abs(-xc--265)); % x from -500 to -100
[~,idx_y_min_tent] = min(abs(-yc-65)); [~,idx_y_max_tent] = min(abs(-yc-50)); % y from 0 to 250
freeboard_tent = freeboard(idx_x_min_tent:idx_x_max_tent,idx_y_min_tent:idx_y_max_tent); % ROV tent fix
freeboard_tent(freeboard_tent > 0.8) = 0.2;
freeboard(idx_x_min_tent:idx_x_max_tent,idx_y_min_tent:idx_y_max_tent) = freeboard_tent;
fb = freeboard(idx_x_min:idx_x_max,idx_y_min:idx_y_max)'; % selected area Z
[XX,YY] = meshgrid(xv, yv);
[in,~] = inpolygon(XX(:),YY(:),FYI_area(:,1),FYI_area(:,2)); inside = reshape(in,length(xv),length(yv));
fb_fyi = freeboard(idx_x_min:idx_x_max,idx_y_min:idx_y_max)'; fb_fyi(inside==0) = NaN;
fb_fyi_filled = fillmissing2(fb_fyi,"linear"); fb_fyi_filled(inside==0) = NaN;

range = -0.1:0.025:0.4; % freeboard range
contourf(xv,yv,fb_fyi_filled,range,'LabelSpacing',10,'edgecolor','none'); hold on;
x_new{3} = xv; y_new{3} = yv; fb_new{3} = fb_fyi;
% p = plot(FYI_area(:,1),FYI_area(:,2),':','color','k','LineWidth',3.0);
p = text(-405,50,sprintf('fb = %.2f m',nanmean(fb_fyi(:)) )); set(p,'Color','k','HorizontalAlignment','left','FontSize',7);
p = text(-405,60,'sn = 0.05 m'); set(p,'Color','k','HorizontalAlignment','left','FontSize',7);
load('density_data.mat',"batlow"); colormap(batlow);
hXLabel = xlabel('x (m)'); hYLabel = ylabel('y (m)'); set([hXLabel hYLabel gca],'FontSize',7,'FontWeight','normal'); xtickangle(0);
xticks(-410:40:-200); xlim([-410 -210]); yticks(40:20:220); ylim([40 220]);
title(sprintf(['FYI freeboard, ', datestr(t,'dd mmmm')]),'FontSize',8,'FontWeight','normal');

nexttile % 22 July
clearvars -except x_new y_new fb_new
project = 'awi-mosaic-l4-als-vq580-stere-0p50m-20200722T151614-20200722T172748-fv2p0.nc'; t = datetime('22-Jul-2020'); % from Hutter et al., 10.1594/PANGAEA.950896
freeboard = ncread(project,'freeboard'); xc = ncread(project,'xc'); yc = ncread(project,'yc');
dx_new = 0; dy_new = 25; % moving both als scan and ROV FYI patch (+y = up)
xc = xc - 25.5; yc = yc - 51 - dy_new; % corrections for July 22 scan

load('density_data.mat',"area_6");
theta = 0; R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)]; for j = 1:length(area_6); FYI_area(j,:) = R*area_6(j,:)'; end % rotation
FYI_area(:,1) = FYI_area(:,1)-419.7; FYI_area(:,2) = FYI_area(:,2)-60; clearvars j R % motion
theta = 0; R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
for i = 1:length(FYI_area); FYI_area(i,:) = R*FYI_area(i,:)'; end % rotation
FYI_area(:,1) = FYI_area(:,1)-3; FYI_area(:,2) = FYI_area(:,2)-17 + dy_new; % motion

[~,idx_x_min] = min(abs(-xc--100)); [~,idx_x_max] = min(abs(-xc--500)); % x from -500 to -100
[~,idx_y_min] = min(abs(-yc-250)); [~,idx_y_max] = min(abs(-yc-0)); % y from 0 to 250
xv = -xc(idx_x_min:idx_x_max); % selected area X
yv = -yc(idx_y_min:idx_y_max)'; % selected area Y
[~,idx_x_min_tent] = min(abs(-xc--225)); [~,idx_x_max_tent] = min(abs(-xc--240)); % x from -500 to -100
[~,idx_y_min_tent] = min(abs(-yc-80)); [~,idx_y_max_tent] = min(abs(-yc-70)); % y from 0 to 250
freeboard_tent = freeboard(idx_x_min_tent:idx_x_max_tent,idx_y_min_tent:idx_y_max_tent); % ROV tent fix
freeboard_tent(freeboard_tent > 0.8) = 0.2;
freeboard(idx_x_min_tent:idx_x_max_tent,idx_y_min_tent:idx_y_max_tent) = freeboard_tent;
fb = freeboard(idx_x_min:idx_x_max,idx_y_min:idx_y_max)'; % selected area Z
[XX,YY] = meshgrid(xv, yv);
[in,~] = inpolygon(XX(:),YY(:),FYI_area(:,1),FYI_area(:,2)); inside = reshape(in,length(xv),length(yv));
fb_fyi = freeboard(idx_x_min:idx_x_max,idx_y_min:idx_y_max)'; fb_fyi(inside==0) = NaN;
fb_fyi_filled = fillmissing2(fb_fyi,"linear"); fb_fyi_filled(inside==0) = NaN;

% rotation
dx = -33; dy = -61;
XY = [XX(:) YY(:)];
theta=-11; R=[cosd(theta) -sind(theta); sind(theta) cosd(theta)]; rotXY=XY*R'; % rotation
Xqr = reshape(rotXY(:,1),size(XX,1),[]); Yqr = reshape(rotXY(:,2),size(YY,1),[]);
xv_rot = Xqr+dx; yv_rot = Yqr+dy;
for i = 1:length(FYI_area); FYI_area_rot(i,:) = R*FYI_area(i,:)'; end % rotation
FYI_area_rot(:,1) = FYI_area_rot(:,1)+dx; FYI_area_rot(:,2) = FYI_area_rot(:,2)+dy; % motion

range = -0.1:0.025:0.4;
contourf(xv_rot,yv_rot,fb_fyi_filled,range,'LabelSpacing',10,'edgecolor','none'); hold on;
x_new{4} = xv_rot; y_new{4} = yv_rot; fb_new{4} = fb_fyi;
p = text(-405,50,sprintf('fb = %.2f m',nanmean(fb_fyi(:)) )); set(p,'Color','k','HorizontalAlignment','left','FontSize',7);
p = text(-405,60,'sn = 0.05 m'); set(p,'Color','k','HorizontalAlignment','left','FontSize',7);
load('density_data.mat',"batlow"); colormap(batlow);
hBar1 = colorbar; ylabel(hBar1,'Freeboard (m)','FontSize',7); set(hBar1,'Position',[0.893 0.76 0.007 0.15]); clim([0 0.4]); hBar1.Color = "k";
hXLabel = xlabel('x (m)'); hYLabel = ylabel('y (m)'); set([hXLabel hYLabel gca],'FontSize',7,'FontWeight','normal'); xtickangle(0);
xticks(-410:40:-200); xlim([-410 -210]); yticks(40:20:220); ylim([40 220]);
title(sprintf(['FYI freeboard, ', datestr(t,'dd mmmm')]),'FontSize',8,'FontWeight','normal');

clear; clc;
load('density_data.mat',"x_edge","y_edge","area_6","t");
load('density_data.mat',"XX","YY","ZZ"); % ROV sonar data from Katlein et al., doi:10.1594/PANGAEA.945846
[in,~] = inpolygon(XX{1}(:),YY{1}(:),x_edge-4,y_edge);
inside = reshape(in,601,801);
good = ~isnan(ZZ{1}).*~isnan(ZZ{2}).*~isnan(ZZ{3}).*~isnan(ZZ{4}).*~isnan(ZZ{5}).*~isnan(ZZ{6}).*inside;
for i=1:7
    ZZ_clean{i} = ZZ{i};
    ZZ_clean{i}(good==0) = NaN;
end

[in,~] = inpolygon(XX{1}(:),YY{1}(:),area_6(:,1)',area_6(:,2)'); inside = reshape(in,601,801); % draft in FYI
for i=1:7
    ZZ_fyi{i} = ZZ_clean{i};
    ZZ_fyi{i}(inside==0) = NaN;
end

for i = [1 3 6 7]
    nexttile
    x_off = -419.7; y_off = -60;
    XY = [XX{i}(:) YY{i}(:)]; dx = -40; dy = 31;
    theta=-12; R=[cosd(theta) -sind(theta); sind(theta) cosd(theta)]; rotXY=XY*R';
    XX_rot = reshape(rotXY(:,1),size(XX{i},1),[]); YY_rot = reshape(rotXY(:,2),size(YY{i},1),[]);
    XX_rot = XX_rot+dx; YY_rot = YY_rot+dy;
    for j = 1:length(area_6); FYI_area_rot(j,:) = R*area_6(j,:)'; end % rotation
    FYI_area_rot(:,1) = FYI_area_rot(:,1)+dx; FYI_area_rot(:,2) = FYI_area_rot(:,2)+dy; % motion
    range = -1:0.1:2; % Graph accuracy
    contourf(XX_rot+x_off,YY_rot+y_off,-ZZ_fyi{i},range,'-','ShowText','on','LabelSpacing',400,'LineColor','none'); hold on
    % p = plot(FYI_area_rot(:,1)+x_off,FYI_area_rot(:,2)+y_off,':','color','k','LineWidth',3.0);
    p = text(-405,50,sprintf('d = %.2f m',nanmean(-ZZ_fyi{i}(:)) )); set(p,'Color','k','HorizontalAlignment','left','FontSize',7);
    set(gca,'YDir','normal'); % reverse y-axis for image
    clim([0.0 2]);
    hXLabel = xlabel('x (m)'); hYLabel = ylabel('y (m)'); set([hXLabel hYLabel gca],'FontSize',7,'FontWeight','normal'); xtickangle(0);
    load('density_data.mat',"batlow"); colormap(batlow);
        if i == 7
            hBar1 = colorbar; ylabel(hBar1,'Draft (m)','FontSize',7); set(hBar1,'Position',[0.893 0.28 0.007 0.15]); clim([0 2]); hBar1.Color = "k";
        end
    title(sprintf(['FYI draft, ', datestr(t(i),'dd mmmm')]),'FontSize',8,'FontWeight','normal');
    xticks(-410:40:-210);xlim([-410 -210]); yticks(40:20:220); ylim([40 220]);
end
annotation('textbox',[0 .50 0.02 .51],'String','(a)','FontSize',8,'EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.185 .50 0.185 .51],'String','(b)','FontSize',8,'EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.34 .50 0.34 .51],'String','(c)','FontSize',8,'EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.495 .50 0.495 .51],'String','(d)','FontSize',8,'EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0 .51/2 0.02 .51/2],'String','(e)','FontSize',8,'EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.185 .51/2 0.185 .51/2],'String','(f)','FontSize',8,'EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.34 .51/2 0.34 .51/2],'String','(g)','FontSize',8,'EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.495 .51/2 0.495 .51/2],'String','(i)','FontSize',8,'EdgeColor','none','HorizontalAlignment','center');
clearvars batlow dx dy FYI_area_rot good hXLabel hYLabel hBar1 i in inside j p R range theta x_edge x_off y_off rotXY area_6

%% Figure 9: effective ice density (6.2 x 4.0 in)
% ROV data repgridding for ALS scan from 30 June
clear; clc; load('density_data.mat',"X","Y","Z","c");
% water level correction for ROV
Z_cor = Z;
Z_cor{1} = Z{1} + 0.72; Z_cor{2} = Z{2} + 0.82; Z_cor{3} = Z{3} + 0.68; Z_cor{4} = Z{4} + 0.69; Z_cor{5} = Z{5} + 0.75; Z_cor{6} = Z{6} + 0.69; Z_cor{7} = Z{7} + 0.62;
x0 = min(X{1}); y0 = min(Y{1});
x_cor   = [0 +1.5 +4 2.5 +3.0 +1.0-1 -2.0];
y_cor   = [1 -1.5 +0 1.0 -1.0 +0.5-0 +1.5];
x_cor_2 = [5.5 4.0 0 4.5 0 5.0+1.5 0];
for i = 1:7; X{i} = X{i} - x0 - x_cor(i) - x_cor_2(i); Y{i} = Y{i} - y0 - y_cor(i); end
clearvars i x0 y0 Z x_cor y_cor x_cor_2
i = 3; % ROV scan
x = X{i}; y = Y{i}; z = Z_cor{i};
xv = linspace(-50, 350, 2*(350+50)+1); % all X
yv = linspace(+50, 350, 2*(350-50)+1); % all Y
[XX{i},YY{i}] = meshgrid(xv, yv);
ZZ{i} = griddata(x,y,z,XX{i},YY{i});
% rotation
XY = [XX{i}(:) YY{i}(:)];
theta = -7;
R=[cosd(theta) -sind(theta); sind(theta) cosd(theta)]; rotXY=XY*R'; 
XX_rot{i} = reshape(rotXY(:,1), size(XX{i},1), []);
YY_rot{i} = reshape(rotXY(:,2), size(YY{i},1), []);
% shifting
dX = -435.2; dY = -51; XX_rot{i} = XX_rot{i} + dX; YY_rot{i} = YY_rot{i} + dY;
clearvars theta R x y z Z_cor xv yv rotXY XY X Y 
% FYI area
load('density_data.mat',"area_6");
theta = -12; R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)]; for j = 1:length(area_6); FYI_area(j,:) = R*area_6(j,:)'; end % rotation
FYI_area(:,1) = FYI_area(:,1)-460; FYI_area(:,2) = FYI_area(:,2)-29; clearvars j R % motion
[in,~] = inpolygon(XX_rot{i}(:),YY_rot{i}(:),FYI_area(:,1),FYI_area(:,2)); inside = reshape(in,601,801);
ZZ_fyi{i} = ZZ{i};

% FYI freeboard for selected ALS scans
project(1) = "awi-mosaic-l4-als-vq580-stere-0p50m-20200630T074032-20200630T090829-fv2p0.nc"; % from Hutter et al., 10.1594/PANGAEA.950896
t(1) = datetime('30-Jun-2020'); % ncdisp(project); 30 June
project = project(1);
freeboard = ncread(project,'freeboard'); xc = ncread(project,'xc'); yc = ncread(project,'yc');
xc = xc - 0.0; yc = yc - 9.5; % corrections for June 30 scan
   
load('density_data.mat',"area_6");
theta = -12; R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)]; for i = 1:length(area_6); FYI_area(i,:) = R*area_6(i,:)'; end % rotation
FYI_area(:,1) = FYI_area(:,1)-460; FYI_area(:,2) = FYI_area(:,2)-29; % motion
[~,idx_x_min] = min(abs(-xc--200)); [~,idx_x_max] = min(abs(-xc--450)); % x from -500 to -100
[~,idx_y_min] = min(abs(-yc-250)); [~,idx_y_max] = min(abs(-yc-30)); % y from 0 to 250
xv = -xc(idx_x_min:idx_x_max); % selected area X
yv = -yc(idx_y_min:idx_y_max)'; % selected area Y
fb = freeboard(idx_x_min:idx_x_max,idx_y_min:idx_y_max)'; % selected area Z
freeboard(freeboard > 1.0) = 0.2; % ROV tent fix
[XX,YY] = meshgrid(xv, yv);
[in,~] = inpolygon(XX(:),YY(:),FYI_area(:,1),FYI_area(:,2)); inside = reshape(in,length(xv),length(yv));
fb_fyi = freeboard(idx_x_min:idx_x_max,idx_y_min:idx_y_max)'; fb_fyi(inside==0) = NaN;
fb_filled = fillmissing2(freeboard(idx_x_min:idx_x_max,idx_y_min:idx_y_max)',"linear");
fb_filled_cut = fillmissing2(freeboard(idx_x_min:idx_x_max,idx_y_min:idx_y_max)',"linear"); fb_filled_cut(inside==0) = NaN;
fb_fyi_regrid = griddata(xv,yv,fb_fyi,XX_rot{3},YY_rot{3});
fb_filled_regrid = griddata(xv,yv,fb_filled,XX_rot{3},YY_rot{3});
fb_filled_cut_regrid = griddata(xv,yv,fb_filled_cut,XX_rot{3},YY_rot{3}); fb_30jun = fb_filled_cut_regrid; x_30jun = XX_rot{3}; y_30jun = YY_rot{3};

% Figure: freeboard and derived density, 30 June
step = 20; B = ones(step)/step^2;
fb_filled_regrid_sm = conv2(fb_filled_regrid,B,'same');
fb_filled_regrid_sm(isnan(fb_filled_cut_regrid)) = NaN;
ZZ_fyi_sm = conv2(ZZ_fyi{3},B,'same');
di = -ZZ_fyi_sm; fb_sn = fb_filled_regrid_sm; rho_sn = 400;
h_sn = 0.08;
fb_i = fb_sn - (fb_sn-min(fb_sn(:)))/nanmean((fb_sn(:)-min(fb_sn(:))))*h_sn;
h_sn_cor = fb_sn - fb_i;
hi = di + fb_i; rhow = 1017.5;
rhosi = (rhow*di - rho_sn*h_sn_cor) ./ hi;

figure
tile = tiledlayout(2,3); tile.TileSpacing = 'compact'; tile.Padding = 'compact';
ax(1) = nexttile;
range = 0.0:0.05:0.4; % freeboard range
i = 3; contourf(XX_rot{i},YY_rot{i},fb_filled_cut_regrid,range,'-','ShowText','on','LabelSpacing',400,'LineColor','none'); hold on
p = text(-405,50,sprintf('%.2f ± %.2f m',nanmean(fb_filled_cut_regrid(:)),nanstd(fb_filled_cut_regrid(:)) )); set(p,'Color','k','HorizontalAlignment','left','FontSize',8);
load('density_data.mat',"batlow"); colormap(ax(1),batlow);
hBar1 = colorbar; ylabel(hBar1,'Freeboard (m)','FontSize',7); clim([0 0.4]); set(hBar1,'Position',[0.24 0.76 0.007 0.15]); hBar1.Color = "k";
title('FYI freeboard, 30 June','FontSize',8,'FontWeight','Normal');
hXLabel = xlabel('x (m)'); hYLabel = ylabel('y (m)'); set([hXLabel hYLabel gca],'FontSize',8,'FontWeight','normal'); xtickangle(0);
xticks(-410:40:-200); xlim([-410 -210]); yticks(40:20:220); ylim([40 220]); % xlim([-410 -220]); ylim([40 220]);

ax(2) = nexttile; % density
range = 800:2:1000; % Graph accuracy
i = 3; contourf(XX_rot{i},YY_rot{i},rhosi,range,'-','ShowText','on','LabelSpacing',400,'LineColor','none'); hold on
p = text(-405,50,sprintf('%.0f ± %.0f kg m^-^3',nanmean(rhosi(:)),nanstd(rhosi(:)) )); set(p,'Color','k','HorizontalAlignment','left','FontSize',8);
colormap(ax(2),batlow);
hBar1 = colorbar; ylabel(hBar1,'FYI density (kg m^-^3)','FontSize',7); clim([800 1000]); set(hBar1,'Position',[0.55 0.76 0.007 0.15]); hBar1.Color = "k";
title('FYI density, 30 June','FontSize',8,'FontWeight','Normal');
hXLabel = xlabel('x (m)'); hYLabel = ylabel('y (m)'); set([hXLabel hYLabel gca],'FontSize',8,'FontWeight','normal'); xtickangle(0);
xticks(-410:40:-200); xlim([-410 -210]); yticks(40:20:220); ylim([40 220]); % xlim([-410 -220]); ylim([40 220]);

ax(3) = nexttile;
histogram(rhosi,40,'BinLimits',[800 1000],'Normalization','probability','FaceColor',[0.7 0.7 0.7]); hold on
p = plot(910.0,0.0185,'ko','color',c{1}); p.MarkerSize = 3; set(p,'markerfacecolor',get(p,'color')); % weighing
p = plot(899.8,0.0185,'ko','color',c{1}); p.MarkerSize = 3; set(p,'markerfacecolor',get(p,'color')); % weighing
leg = legend('FYI ROV','Weighing','box','off','NumColumns',1); set(leg,'FontSize',7,'Location','west'); leg.ItemTokenSize = [30*0.15,18*0.15];
xlim([820 960]); xticks(820:20:960); title('FYI density, 30 June','FontSize',8,'FontWeight','Normal');
hXLabel = xlabel('FYI density (kg m^-^3)'); hYLabel = ylabel('pdf'); set([hXLabel hYLabel gca],'FontSize',8,'FontWeight','normal'); xtickangle(0);

annotation('textbox',[0 .505 0.02 .51],'String','(a)','FontSize',9,'EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.22 .505 0.23 .51],'String','(b)','FontSize',9,'EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.43 .505 0.44 .51],'String','(c)','FontSize',9,'EdgeColor','none','HorizontalAlignment','center');

% FYI freeboard for selected ALS scans (copy to make corrections) - 22 July co-location
clearvars -except fb_30jun x_30jun y_30jun ZZ_fyi batlow
project = 'awi-mosaic-l4-als-vq580-stere-0p50m-20200722T151614-20200722T172748-fv2p0.nc'; % from Hutter et al., 10.1594/PANGAEA.950896
freeboard = ncread(project,'freeboard'); xc = ncread(project,'xc'); yc = ncread(project,'yc');
dx_new = 0; dy_new = 25; % moving both als scan and ROV FYI patch (+y = up)
xc = xc - 25.5; yc = yc - 51 - dy_new; % corrections for July 22 scan
    
load('density_data.mat',"area_6");
theta = 0; R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)]; for j = 1:length(area_6); FYI_area(j,:) = R*area_6(j,:)'; end % rotation
FYI_area(:,1) = FYI_area(:,1)-419.7; FYI_area(:,2) = FYI_area(:,2)-60; clearvars j R % motion
theta = 0; R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
for i = 1:length(FYI_area); FYI_area(i,:) = R*FYI_area(i,:)'; end % rotation
FYI_area(:,1) = FYI_area(:,1)-3; FYI_area(:,2) = FYI_area(:,2)-17 + dy_new; % motion

[~,idx_x_min] = min(abs(-xc--100)); [~,idx_x_max] = min(abs(-xc--500)); % x from -500 to -100
[~,idx_y_min] = min(abs(-yc-250)); [~,idx_y_max] = min(abs(-yc-0)); % y from 0 to 250
xv = -xc(idx_x_min:idx_x_max); % selected area X
yv = -yc(idx_y_min:idx_y_max)'; % selected area Y
fb = freeboard(idx_x_min:idx_x_max,idx_y_min:idx_y_max)'; % selected area Z
[~,idx_x_min_tent] = min(abs(-xc-225)); [~,idx_x_max_tent] = min(abs(-xc--240)); % x from -500 to -100
[~,idx_y_min_tent] = min(abs(-yc-80)); [~,idx_y_max_tent] = min(abs(-yc-70)); % y from 0 to 250
freeboard_tent = freeboard(idx_x_min_tent:idx_x_max_tent,idx_y_min_tent:idx_y_max_tent); % ROV tent fix
freeboard_tent(freeboard_tent > 0.8) = 0.2;
freeboard(idx_x_min_tent:idx_x_max_tent,idx_y_min_tent:idx_y_max_tent) = freeboard_tent;
[XX,YY] = meshgrid(xv, yv);
[in,~] = inpolygon(XX(:),YY(:),FYI_area(:,1),FYI_area(:,2)); inside = reshape(in,length(xv),length(yv));
fb_fyi = freeboard(idx_x_min:idx_x_max,idx_y_min:idx_y_max)'; fb_fyi(inside==0) = NaN;
fb_filled = fillmissing2(freeboard(idx_x_min:idx_x_max,idx_y_min:idx_y_max)',"linear");
fb_filled_cut = fillmissing2(freeboard(idx_x_min:idx_x_max,idx_y_min:idx_y_max)',"linear"); fb_filled_cut(inside==0) = NaN;

load('density_data.mat',"X","Y","Z","c"); % ROV sonar data from Katlein et al., doi:10.1594/PANGAEA.945846
% water level correction for ROV
Z_cor = Z;
Z_cor{1} = Z{1} + 0.72; Z_cor{2} = Z{2} + 0.82; Z_cor{3} = Z{3} + 0.68; Z_cor{4} = Z{4} + 0.69;
Z_cor{5} = Z{5} + 0.75; Z_cor{6} = Z{6} + 0.75-0.06; Z_cor{7} = Z{7} + 0.62;
x0 = min(X{1}); y0 = min(Y{1});
x_cor   = [0 +1.5 +4 2.5 +3.0 +1.0-1 -2.0];
y_cor   = [1 -1.5 +0 1.0 -1.0 +0.5-0 +1.5];
x_cor_2 = [5.5 4.0 0 4.5 0 5.0+1.5 0];
for i = 1:7; X{i} = X{i} - x0 - x_cor(i) - x_cor_2(i); Y{i} = Y{i} - y0 - y_cor(i); end
clearvars i x0 y0 Z x_cor y_cor x_cor_2
i = 6; % ROV scan
x = X{i}; y = Y{i} + dy_new; z = Z_cor{i};
xv_rov = linspace(-50, 350, 2*(350+50)+1); % all X
yv_rov = linspace(+50, 350, 2*(350-50)+1); % all Y
[XX_rov{i},YY_rov{i}] = meshgrid(xv_rov, yv_rov);
ZZ{i} = griddata(x,y,z,XX_rov{i},YY_rov{i});
% rotation
XY = [XX_rov{i}(:) YY_rov{i}(:)];
theta = 0;
R=[cosd(theta) -sind(theta); sind(theta) cosd(theta)]; rotXY=XY*R'; 
XX_rot{i} = reshape(rotXY(:,1), size(XX_rov{i},1), []);
YY_rot{i} = reshape(rotXY(:,2), size(YY_rov{i},1), []);
% shifting
dX = -419.7-3; dY = -60-17; XX_rot{i} = XX_rot{i} + dX; YY_rot{i} = YY_rot{i} + dY;
clearvars theta R x y z Z_cor xv_rov yv_rov rotXY XY X Y 
i = 6; [in,~] = inpolygon(XX_rot{i}(:),YY_rot{i}(:),FYI_area(:,1),FYI_area(:,2)); inside = reshape(in,601,801);
ZZ_fyi{i} = ZZ{i};

ax(4) = nexttile; % regridded ALS freeboard
fb_filled_regrid = griddata(xv,yv,fb_filled,XX_rot{6},YY_rot{6});
fb_fyi_regrid = griddata(xv,yv,fb_fyi,XX_rot{6},YY_rot{6});
fb_filled_cut_regrid = griddata(xv,yv,fb_filled_cut,XX_rot{6},YY_rot{6});

% new rotation
dx = -33; dy = -61;
XY = [XX_rot{i}(:) YY_rot{i}(:)];
theta=-11; R=[cosd(theta) -sind(theta); sind(theta) cosd(theta)]; rotXY=XY*R'; % rotation
Xqr = reshape(rotXY(:,1),size(XX_rot{i},1),[]); Yqr = reshape(rotXY(:,2),size(YY_rot{i},1),[]);
xv_rot = Xqr+dx; yv_rot = Yqr+dy;
range = 0.0:0.05:0.4; % freeboard range
contourf(xv_rot,yv_rot,fb_filled_cut_regrid,range,'-','ShowText','on','LabelSpacing',400,'LineColor','none'); hold on
p = text(-405,50,sprintf('%.2f ± %.2f m',nanmean(fb_filled_cut_regrid(:)),nanstd(fb_filled_cut_regrid(:)) )); set(p,'Color','k','HorizontalAlignment','left','FontSize',8);
colormap(ax(4),batlow);
hBar1 = colorbar; ylabel(hBar1,'Freeboard (m)','FontSize',7); clim([0 0.4]); set(hBar1,'Position',[0.24 0.28 0.007 0.15]); hBar1.Color = "k";
title('FYI freeboard, 22 July','FontSize',8,'FontWeight','Normal');
hXLabel = xlabel('x (m)'); hYLabel = ylabel('y (m)'); set([hXLabel hYLabel gca],'FontSize',8,'FontWeight','normal'); xtickangle(0);
xticks(-410:40:-200); xlim([-410 -210]); yticks(40:20:220); ylim([40 220]);

ax(5) = nexttile; % density
step = 20; B = ones(step)/step^2;
fb_filled_regrid_sm = conv2(fb_filled_regrid,B,'same');
fb_filled_regrid_sm(isnan(fb_filled_cut_regrid)) = NaN;
ZZ_fyi_sm = conv2(ZZ_fyi{6},B,'same');
di = -ZZ_fyi_sm; fb_sn = fb_filled_regrid_sm; rho_sn = 400;
h_sn = 0.05;
fb_i = fb_sn - (fb_sn-min(fb_sn(:)))/nanmean((fb_sn(:)-min(fb_sn(:))))*h_sn;
h_sn_cor = fb_sn - fb_i;
di = di - 0.0214; % 0.9976 m - interpolated draft - correction
hi = di + fb_i; rhow = 1017.5;
rhosi = (rhow*di - rho_sn*h_sn_cor) ./ hi;
range = 800:2:1000; % Graph accuracy
contourf(xv_rot,yv_rot,rhosi,range,'-','ShowText','on','LabelSpacing',400,'LineColor','none'); hold on
p = text(-405,50,sprintf('%.0f ± %.0f kg m^-^3',nanmean(rhosi(:)),nanstd(rhosi(:)) )); set(p,'Color','k','HorizontalAlignment','left','FontSize',8);
colormap(ax(5),batlow);
hBar1 = colorbar; ylabel(hBar1,'FYI density (kg m^-^3)','FontSize',7); clim([800 1000]); set(hBar1,'Position',[0.55 0.28 0.007 0.15]); hBar1.Color = "k";
title('FYI density, 22 July','FontSize',8,'FontWeight','Normal');
hXLabel = xlabel('x (m)'); hYLabel = ylabel('y (m)'); set([hXLabel hYLabel gca],'FontSize',8,'FontWeight','normal'); xtickangle(0);
xticks(-410:40:-200); xlim([-410 -210]); yticks(40:20:220); ylim([40 220]);

nexttile % density histogram
histogram(rhosi,40,'BinLimits',[800 1000],'Normalization','probability','FaceColor',[0.7 0.7 0.7]); hold on
p = plot(883.6,0.0185,'ko','color',c{1}); p.MarkerSize = 3; set(p,'markerfacecolor',get(p,'color')); % weighing
leg = legend('FYI ROV','Weighing','box','off','NumColumns',1); set(leg,'FontSize',7,'Location','west'); leg.ItemTokenSize = [30*0.15,18*0.15];
xlim([820 960]); xticks(820:20:960); ylim([0 0.02]); title('FYI density, 22 July','FontSize',8,'FontWeight','Normal');
hXLabel = xlabel('FYI density (kg m^-^3)'); hYLabel = ylabel('pdf'); set([hXLabel hYLabel gca],'FontSize',8,'FontWeight','normal'); xtickangle(0);

annotation('textbox',[0.00 .51/2 0.02 .51/2],'String','(d)','FontSize',9,'EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.22 .51/2 0.23 .51/2],'String','(e)','FontSize',9,'EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.43 .51/2 0.44 .51/2],'String','(f)','FontSize',9,'EdgeColor','none','HorizontalAlignment','center');
clearvars batlow area_6 ax B dx dX dx_new dy dY dy_new i idx_x_max idx_x_max_tent idx_x_min idx_x_min_tent idx_y_max idx_y_max_tent idx_y_min idx_y_min_tent in inside p R project
clearvars hXLabel hYLabel hBar1 range rotXY step theta