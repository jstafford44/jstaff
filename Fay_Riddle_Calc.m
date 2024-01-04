classdef Fay_Riddle_Calc
   methods
       function heating = heat_calc(obj,trajectoryData)
            addpath('EquilFlowM')
            r_n = 0.23;
            Le = 1; h_d = 0;
            % stagnation point heating code
            gamma = 1.4;
            R = trajectoryData(:,4)./(trajectoryData(:,5).*trajectoryData(:,6));
            Mach = trajectoryData(:,3)./(sqrt(gamma.*R.*trajectoryData(:,6)));
            heating = zeros(length(trajectoryData(:,3)),1);
            T_w = 300;
            for i = 1:length(Mach)
                if(Mach(i)>1)
                    postNormalShock =clib.EquilFlowM.EquilFlow.NormalShock(trajectoryData(i,2)/1000,Mach(i));
                    h_e = postNormalShock(20)*1000; %J/kg/K
                    rho_e = postNormalShock(22);
                    mu_e = postNormalShock(24);
                    Pr = postNormalShock(27);
                    % Le = postNormalShock(28);%not 1?
                    P_e = postNormalShock(29);
                    vel = postNormalShock(30);
                    P_inf = postNormalShock(31);
                    rho_w = P_e/(287*T_w);
                    mu_w = 0.000001458*(T_w)^1.5/(T_w+110.4);
                    dU_dx = 1/r_n*sqrt(2*(P_e-P_inf)/rho_e);
                    h_w = 1.0045*300*1000;
                    h_oe = h_e + vel^2/2;
                    heating(i) = (0.76*Pr^-.6*(rho_e*mu_e)^.4*(rho_w*mu_w)^.1*sqrt(dU_dx)*(h_oe-h_w)*(1+(Le^.52-1)*h_d/h_oe))/10000;
                else
                    conditions = clib.EquilFlowM.EquilFlow.Composition(1,trajectoryData(i,6),trajectoryData(i,4));
                    P_e = trajectoryData(i,4);
                    P_w = P_e;
                    rho_e = trajectoryData(i,5);
                    mu_e = conditions(24);
                    Pr = conditions(27);
                    rho_w = P_w/(287*T_w);
                    dU_dx = 1/r_n*sqrt(2*(P_e-P_inf)/rho_e);
                    h_w = 1.0045*300*1000;
                    h_oe = h_e + vel^2/2;
                    heating(i) = (0.76*Pr^-.6*(rho_e*mu_e)^.4*(rho_w*mu_w)^.1*sqrt(dU_dx)*(h_oe-h_w)*(1+(Le^.52-1)*h_d/h_oe))/10000;
                end
                if(heating(i)<0)
                    heating(i) = 0;
                end
            end
       end
   end
end