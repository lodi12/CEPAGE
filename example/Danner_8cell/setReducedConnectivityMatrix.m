function [g_in,g_ex,g_el] = setReducedConnectivityMatrix(Ntot,gSynIn,gSynEx,da)




p1 = 115.98;
  p2 = -231.71;
  p3 = 25.54;
  p4 = 329.37;
  p5 = -407.13;
  p6 = 235.88;
  p7 = -76.053;
  p8 = 13.751;
  p9 = -1.1155;
  p10 = 0.11545;
  p11 = 0.16808;

  g_ex_fore0 = p1*da.^10 + p2*da.^9 + p3*da.^8 + p4*da.^7 + p5*da.^6 + p6*da.^5+...
            p7*da.^4+p8*da.^3+p9*da.^2+p10*da+p11;
  
% g_ex_fore0 = p1*da.^6 + p2*da.^5 + p3*da.^4 + p4*da.^3 + p5*da.^2 + p6*da + p7;


  
p1 = 3058.8;
  p2 = -13011;
  p3 = 23662;
  p4 = -23916;
  p5 = 14651;
  p6 = -5568.3;
  p7 = 1292;
  p8 = -172.9;
  p9 = 12.005;
  p10 = -0.25126;
  p11 = 0.1689;
  
g_ex_hind0 = p1*da.^10 + p2*da.^9 + p3*da.^8 + p4*da.^7 + p5*da.^6 + p6*da.^5+...
            p7*da.^4+p8*da.^3+p9*da.^2+p10*da+p11;


% flag to activate syn
V0D = 1;
V0V = 1; %V0V del = 0.1;
V3 = 1.3;
trans = 1;
% Create Network
    

% if da > 0.3
%     V3 = 1;
% end

    g_in = zeros(Ntot);
    g_ex = zeros(Ntot);
    g_el = zeros(Ntot);

    % Fore
    g_in(1,2) = 0.2853*0.08*gSynIn; %HCO
    g_in(2,1) = 0.2853*1*gSynIn; %HCO
    
    g_in(1,4) = 0.1841*0.012*gSynIn*trans; %trans
    g_in(1,3) = 0.2850*0.0266*gSynIn*V0D; % V0D
    g_in(1,3) = g_in(1,3)+0.1113*0.2*gSynIn*V0V; % V0V
    g_in(1,5) = 0.3633*0.035*gSynIn; % Homolat
    
    g_in(3,4) = 0.2853*0.08*gSynIn; %HCO
    g_in(4,3) = 0.2853*1*gSynIn; %HCO
    
    g_in(3,2) = 0.1841*0.012*gSynIn*trans; %trans
    g_in(3,1) = 0.2850*0.0266*gSynIn*V0D; % V0D
    g_in(3,1) = g_in(3,1)+0.1113*0.2*gSynIn*V0V; % V0V
    g_in(3,7) = 0.3633*0.035*gSynIn; % Homolat
    
    
    % Hind
    g_in(5,6) = 0.2853*0.08*gSynIn; %HCO
    g_in(6,5) = 0.2853*1*gSynIn; %HCO
    
    g_in(5,8) = 0.1841*0.017*gSynIn*trans; %trans
    g_in(5,7) = 0.2852*0.04*gSynIn*V0D; % V0D
    g_in(5,7) = g_in(5,7)+0.1115*0.3*gSynIn*V0V; % V0V
    g_in(5,1) = 0.3635*0.015*gSynIn; % Homolat
    
    g_in(7,8) = 0.2853*0.08*gSynIn; %HCO
    g_in(8,7) = 0.2853*1*gSynIn; %HCO
    
    g_in(7,6) = 0.1841*0.017*gSynIn*trans; %trans
    g_in(7,5) = 0.2852*0.04*gSynIn*V0D; % V0D
    g_in(7,5) = g_in(7,5)+0.1115*0.3*gSynIn*V0V; % V0V
    g_in(7,3) = 0.3635*0.015*gSynIn; % Homolat

   
  
    g_ex(1,3) = g_ex_fore0*0.02*gSynEx*V3;
    g_ex(3,1) = g_ex_fore0*0.02*gSynEx*V3;
    
    g_ex(5,7) = g_ex_hind0*0.03*gSynEx*V3;
    g_ex(7,5) = g_ex_hind0*0.03*gSynEx*V3;