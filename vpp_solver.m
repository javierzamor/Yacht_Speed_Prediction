

function vpp_solver()
    set(0,'DefaultFigureWindowStyle','docked')
    clc
    close('all');
    fclose all;
    clear;
    
   
    [param] = load_data();
    testing_data(param);
    
    V_tw = param.V_tw;
    alfa_tw = param.alfa_tw;
    tic
    x= zeros(5,length(V_tw)*length(alfa_tw));
    for i = 1:length(V_tw)       
        for j =1:length(alfa_tw)
            [V,phi,contador]=solver(param,V_tw(i),alfa_tw(j));
            x(:,(i-1)*length(V_tw)+j)=[V_tw(i);alfa_tw(j);V;phi;(contador/1000)];           
        end
    end
    toc

    postprocessing(x,param);
   
end


function [V_o,phi_o,contador]=solver(param,V_tw,alfa_tw)

    err=1;
    contador=1;
    V_o=1.25*sqrt(param.LWL);   
    phi_o=0;
   
    if V_tw <= V_o;
        V_o=V_tw;
    end
    
    while(err>0.0001)
        
        V=V_o;
        phi=phi_o;
        Fn=V/sqrt(param.g*param.LWL);
        if (Fn>0.60)           
            V=0.6*rand(1)*sqrt(param.g*param.LWL);
        elseif (Fn<0)
            V=0.6*rand(1)*sqrt(param.g*param.LWL);
        end
        if (phi<0)
            phi=35*rand(1);
        elseif (phi>35)
            phi=35*rand(1);
        end
      
        
        [B2,B4,C2,C4,m_heel,f_heel,dFheel_dphi,dFheel_dV,f_drive]=aero(V,phi,...
            V_tw,alfa_tw,param);
        [B3,C3,m_right] = M_right(phi,param);
        [B1,C1,r_tot] = Rtot(V,phi,param,f_heel,dFheel_dphi,dFheel_dV);
        
        
        D1=r_tot-V*B1-(phi*pi/180)*C1;
        D2=f_drive-V*B2-(phi*pi/180)*C2;
        D3=m_right-V*B3-(phi*pi/180)*C3;
        D4=m_heel-V*B4-(phi*pi/180)*C4;

%         AAA=[B1-B2 C1-C2;B3-B4 C3-C4];   
%         BBB=[D2-D1;D4-D3];
%         
%         spectral_radius=max(abs(eig(AAA/max(max(AAA)))));
      
        phi_o=((D2-D1)-((B1-B2)/B4)*(D3-D4))/((C1-C2)+((B1-B2)/B4)*(C3-C4));        
        V_o=((D4-D3)-(C3-C4)*phi_o)/(-B4);
        phi_o=phi_o*180/pi;
 
        err_1=abs(V_o-V);
        err_2=abs(phi_o-phi);
        err=sqrt(err_1^2+err_2^2);
        
        if contador==1000
            V_o=NaN;
            phi_o=NaN;
            return;          
        end
        contador=contador+1;
    end
end

function [B1,C1,r_tot] = Rtot(V,phi,param,f_heel,dFheel_dphi,dFheel_dV)

    V = abs(V); 
    phi = abs(phi);
    
    [dri_dV,dri_dphi,ri]= Ri(V,phi,param,f_heel,dFheel_dphi,dFheel_dV);
    [drrh_dV,drrh_dphi,rrh]= Rrh(V,param);
    [drrhH_dV,drrhH_dphi,rrhH]= RrhH(V,phi,param);
    [drrk_dV,drrk_dphi,rrk]= Rrk(V,param);
    [drrkH_dV,drrkH_dphi,rrkH]= RrkH(V,phi,param);
    [drvh_dV,drvh_dphi,rvh]= Rvh(V,param);
    [drvhH_dV,drvhH_dphi,rvhH]= RvhH(V,phi,param);
    [drvk_dV,drvk_dphi,rvk]=Rvk(V,param);
    [drvr_dV,drvr_dphi,rvr]=Rvr(V,param);

    r_tot = ri + rrh + rrhH + rrk + rrkH + rvh + rvhH + rvk + rvr;
    C1=dri_dphi+drrh_dphi+drrhH_dphi+drrk_dphi+drrkH_dphi+drvh_dphi+drvhH_dphi+...
        drvk_dphi+drvr_dphi;
    B1=dri_dV+drrh_dV+drrhH_dV+drrk_dV+drrkH_dV+drvh_dV+drvhH_dV+drvk_dV+drvr_dV;
   
end

function [drvk_dV,drvk_dphi,rvk] = Rvk (V,param)

    Rn = param.CHMEK* V/ param.ni_w;
    Cf = 0.075/ ( log10(Rn) - 2 )^2;
    Rfk = 1/2* param.rho_w* V.^2* param.SK* Cf;
    rvk = Rfk* param.KEELFF;
    
    drvk_dphi=0;
    
    dRn_dV=param.CHMEK/ param.ni_w;
    dlogRn_dV=(log10(exp(1))/Rn)*dRn_dV;
    dCf_dV=-2*(0.075/(log10(Rn)-2)^3)*dlogRn_dV;
    dV2Cf_dV=2*V*Cf+V^2*dCf_dV;
      
    dRfk_dV=0.5*param.rho_w*param.SK*dV2Cf_dV;
    drvk_dV=param.KEELFF*dRfk_dV;
end

function [drvr_dV,drvr_dphi,rvr] = Rvr (V,param)
    
    Rn = param.CHMER* V/ param.ni_w;
    Cf = 0.075/ ( log10(Rn) - 2 )^2;
    Rfr = 1/2* param.rho_w* V^2* param.SR* Cf;
    rvr = Rfr* param.RUDDFF;
    
    drvr_dphi=0;
    
    dRn_dV=param.CHMER/ param.ni_w;
    dlogRn_dV=(log10(exp(1))/Rn)*dRn_dV;
    dCf_dV=-2*(0.075/(log10(Rn)-2)^3)*dlogRn_dV;
    dV2Cf_dV=2*V*Cf+V^2*dCf_dV; 
    dRfr_dV=0.5*param.rho_w*param.SR*dV2Cf_dV;
    drvr_dV=param.RUDDFF*dRfr_dV;
    
    
end

function [drvhH_dV,drvhH_dphi,rvhH] = RvhH(V,phi,param)
   
    if (phi>=0 && phi<5)
        s0=-0.8224*phi;
        ds0_dphi=-0.8224;
        s1=0.0108*phi;
        ds1_dphi=0.0108;
        s2=-0.0054*phi;
        ds2_dphi=-0.0054;
        s3=1.2658*phi;
        ds3_dphi=1.2658;
    elseif (phi>=5 && phi<10)
        s0=-0.082*phi-3.702;
        ds0_dphi=-0.082;
        s1=-0.0372*phi+0.24;
        ds1_dphi=-0.0372;
        s2=-0.01*phi+0.023;
        ds2_dphi=-0.01;
        s3=0.4818*phi+3.92;
        ds3_dphi=0.4818;
    elseif (phi>=10 && phi<15)
        s0=0.2462*phi-6.984;
        ds0_dphi=0.2462;
        s1=-0.0514*phi-0.382;
        ds1_dphi=-0.0514;
        s2=-0.0082*phi+0.005;
        ds2_dphi=-0.0082;
        s3=0.0422*phi+8.316;
        ds3_dphi=0.0422;
    elseif (phi>=15 && phi<20)
        s0=1.0282*phi-18.714;
        ds0_dphi=1.0282;
        s1=-0.1622*phi+2.044;
        ds1_dphi=-0.1622;
        s2=0.0018*phi-0.145;
        ds2_dphi=0.0018;
        s3=-0.717*phi+19.704;
        ds3_dphi=-0.717;
    elseif (phi>=20 && phi <25)
        s0=0.932*phi-16.79;
        ds0_dphi=0.932;
        s1=-0.221*phi+3.22;
        ds1_dphi=-0.221;
        s2=0.0086*phi-0.281;
        ds2_dphi=0.0086;
        s3=-0.3842*phi+13.048;
        ds3_dphi=-0.3842;
    elseif (phi>=25 && phi<30)
        s0=1.1648*phi-22.61;
        ds0_dphi=1.1648;
        s1=-0.3212*phi+5.725;
        ds1_dphi=-0.3212;
        s2=0.018*phi-0.516;
        ds2_dphi=0.018;
        s3=-0.3352*phi+11.823;
        ds3_dphi=-0.3352;
    elseif (phi>=30 && phi <=35)
        s0=0.4628*phi-1.55;
        ds0_dphi=0.4628;
        s1=-0.2542*phi+3.715;
        ds1_dphi=-0.2542;
        s2=0.0156*phi-0.444;
        ds2_dphi=0.0156;
        s3=0.346*phi-8.613;
        ds3_dphi=0.346;
    end
  
   
    SC=(1.97+0.171*param.BWL/param.TCAN)*sqrt(param.DIVCAN*param.LWL)*...
        (0.65/param.CMS)^(1/3);
    
    Scphi_tmp=s0+s1*(param.BWL./param.TCAN)+s2*((param.BWL./param.TCAN).^2)+...
        s3*(param.CMS);
    Scphi=SC*(1+1/100*Scphi_tmp);

    Rn = param.LWL* 0.7* V/ param.ni_w;
    Cf = 0.075/ ( log10(Rn) - 2 )^2;
    RfhH = 1/2* param.rho_w* V.^2* Cf* (Scphi - SC);
    rvhH = RfhH*(1+ param.HULLFF);
    
    Scphi_tmp=ds0_dphi+ds1_dphi*(param.BWL/param.TCAN)+...
        ds2_dphi*((param.BWL/param.TCAN)^2)+ds3_dphi*(param.CMS);
    dScphi_dphi=SC*(1+1/100*Scphi_tmp);
    
    drfhH_dphi=0.5*param.rho_w*V^2*Cf*dScphi_dphi;
    drvhH_dphi=(1+param.HULLFF)*drfhH_dphi;
    
    
    dRn_dV=0.7*param.LWL/param.ni_w;
    dlogRn_dV=(log10(exp(1))/Rn)*dRn_dV;
    dCf_dV=-2*(0.075/(log10(Rn)-2)^3)*dlogRn_dV;
    dV2Cf_dV=2*V*Cf+V^2*dCf_dV;
    drfhH_V=(1/2)*param.rho_w*(Scphi - SC)*dV2Cf_dV;
    drvhH_dV=(1+param.HULLFF)*drfhH_V;

end

function [drvh_dV,drvh_dphi,rvh] = Rvh(V,param)
   


    SC=(1.97+0.171*param.BWL/param.TCAN)*sqrt(param.DIVCAN*param.LWL)*...
        (0.65/param.CMS)^(1/3);
    Rn = param.LWL* 0.7* V/ param.ni_w;
    Cf = 0.075/ ( log10(Rn) - 2 )^2;
    Rfh = 1/2* param.rho_w* V^2* SC* Cf;
    rvh = Rfh*(1+param.HULLFF);
    
    drvh_dphi=0;
    
    dRn_dV=0.7*param.LWL/param.ni_w;
    dlogRn_dV=(log10(exp(1))/Rn)*dRn_dV;
    dCf_dV=-2*(0.075/(log10(Rn)-2)^3)*dlogRn_dV;
    dV2Cf_dV=2*V*Cf+V^2*dCf_dV;
    dRfh_dV=0.5*param.rho_w*SC*dV2Cf_dV;
    drvh_dV=(1+param.HULLFF)*dRfh_dV;
            
end

function [drrkH_dV,drrkH_dphi,rrkH] = RrkH(V,phi,param)

    Fn = V./sqrt(param.g*param.LWL);
    H1=-3.5837;
    H2=-0.0518;
    H3=0.5958;
    H4=0.2055;
    rrkH_tmp_1=H1*(param.TCAN./param.T)+H2*(param.BWL./param.TCAN);
    rrkH_tmp_2=H3*((param.BWL.*param.TCAN)/(param.T.*param.TCAN))+...
        H4*(param.LWL./param.DIVCAN.^1/3);
    rrkH_tmp_3=rrkH_tmp_1+rrkH_tmp_2;
    rrkH = param.DVK.*param.rho_w.*param.g.*rrkH_tmp_3*Fn.^2.*(phi*pi/180);
    
    drrkH_dphi=param.DVK.*param.rho_w.*param.g.*rrkH_tmp_3*Fn.^2;
    
    dFn2_dV=2*V/(param.g*param.LWL);
    drrkH_dV=param.DVK.*param.rho_w.*param.g.*rrkH_tmp_3*(phi*pi/180)*dFn2_dV;

end

function [drrk_dV,drrk_dphi,rrk] = Rrk(V,param)

    Fn = V/sqrt(param.g*param.LWL);
    
    if (Fn>=0 && Fn < 0.15)
        A0=0;
        A1=0;
        A2=0;
        A3=0;
        dA0_dV=0;
        dA1_dV=0;
        dA2_dV=0;
        dA3_dV=0;
    elseif (Fn>=0.15 && Fn < 0.20)
        A0=-0.0208*Fn+0.00312;
        A1=0.0344*Fn-0.00516;
        A2=0.0234*Fn-0.00351;
        A3=-0.0016*Fn+0.00024;
        dA0_dV=(sqrt(param.g*param.LWL))*(-0.0208);
        dA1_dV=(sqrt(param.g*param.LWL))*(0.0344);
        dA2_dV=(sqrt(param.g*param.LWL))*(0.0234);
        dA3_dV=(sqrt(param.g*param.LWL))*(-0.0016);
    elseif (Fn>=0.20 && Fn<0.25)
        A0=-0.0892*Fn+0.0168;
        A1=0.085*Fn-0.01528;
        A2=0.0546*Fn-0.00975;
        A3=-0.0002*Fn-0.00004;
        dA0_dV=(sqrt(param.g*param.LWL))*(-0.0892);
        dA1_dV=(sqrt(param.g*param.LWL))*(0.085);
        dA2_dV=(sqrt(param.g*param.LWL))*(0.0546);
        dA3_dV=(sqrt(param.g*param.LWL))*(-0.0002);
    elseif (Fn>=0.25 && Fn<0.30)
        A0=-0.112*Fn-0.0225;
        A1=0.1648*Fn-0.03523;
        A2=-0.0642*Fn+0.01995;
        A3=0.006*Fn-0.00159;
        dA0_dV=(sqrt(param.g*param.LWL))*(-0.112);
        dA1_dV=(sqrt(param.g*param.LWL))*(0.1648);
        dA2_dV=(sqrt(param.g*param.LWL))*(-0.0642);
        dA3_dV=(sqrt(param.g*param.LWL))*(0.006);
    elseif (Fn>=0.30 && Fn < 0.35)
        A0=0.0794*Fn-0.03492;
        A1=0.2422*Fn-0.05845;
        A2=-0.0602*Fn+0.01875;
        A3=0.0036*Fn-0.00087;
        dA0_dV=(sqrt(param.g*param.LWL))*(0.0794);
        dA1_dV=(sqrt(param.g*param.LWL))*(0.2422);
        dA2_dV=(sqrt(param.g*param.LWL))*(-0.0602);
        dA3_dV=(sqrt(param.g*param.LWL))*(0.0036);
    elseif (Fn>=0.35 && Fn < 0.40)
        A0=-0.5736*Fn+0.19363;
        A1=1.2034*Fn-0.39487;
        A2=0.2462*Fn-0.08849;
        A3=-0.0044*Fn+0.00193;
        dA0_dV=(sqrt(param.g*param.LWL))*(-0.5736);
        dA1_dV=(sqrt(param.g*param.LWL))*(1.2034);
        dA2_dV=(sqrt(param.g*param.LWL))*(0.2462);
        dA3_dV=(sqrt(param.g*param.LWL))*(-0.0044);
    elseif (Fn>=0.40 && Fn <0.45)
        A0=0.6222*Fn-0.28469;
        A1=0.5886*Fn-0.14895;
        A2=-0.2126*Fn+0.09503;
        A3=0.0036*Fn-0.00127;
        dA0_dV=(sqrt(param.g*param.LWL))*(0.6222);
        dA1_dV=(sqrt(param.g*param.LWL))*(0.5886);
        dA2_dV=(sqrt(param.g*param.LWL))*(-0.2126);
        dA3_dV=(sqrt(param.g*param.LWL))*(0.0036);
    elseif (Fn>=0.45 && Fn< 0.50)
        A0=0.2046*Fn-0.09677;
        A1=-0.8442*Fn+0.49581;
        A2=1.211*Fn-0.54559;
        A3=-0.0298*Fn+0.01376;
        dA0_dV=(sqrt(param.g*param.LWL))*(0.2046);
        dA1_dV=(sqrt(param.g*param.LWL))*(-0.8442);
        dA2_dV=(sqrt(param.g*param.LWL))*(1.211);
        dA3_dV=(sqrt(param.g*param.LWL))*(-0.0298);
    elseif (Fn>=0.50 && Fn <0.55)
        A0=0.8538*Fn-0.42137;
        A1=-1.3422*Fn+0.74481;
        A2=0.2114*Fn-0.04579;
        A3=0.0158*Fn-0.00904;
        dA0_dV=(sqrt(param.g*param.LWL))*(0.8538);
        dA1_dV=(sqrt(param.g*param.LWL))*(-1.3422);
        dA2_dV=(sqrt(param.g*param.LWL))*(0.2114);
        dA3_dV=(sqrt(param.g*param.LWL))*(0.0158);
    elseif (Fn>=0.55 && Fn<=0.60)
        A0=-0.7602*Fn+0.46633;
        A1=2.7026*Fn-1.47983;
        A2=-0.1278*Fn+0.14077;
        A3=-0.0314*Fn+0.01692;
        dA0_dV=(sqrt(param.g*param.LWL))*(-0.7602);
        dA1_dV=(sqrt(param.g*param.LWL))*(2.7026);
        dA2_dV=(sqrt(param.g*param.LWL))*(-0.1278);
        dA3_dV=(sqrt(param.g*param.LWL))*(-0.0314);
    end
    
    Rrk_tmp=A0+A1*(param.T./param.BWL)+A2*((param.T-param.ZCBK).^3./param.DVK)+...
        A3*(param.DIVCAN./param.DVK);
    rrk=(param.DVK.*param.rho_w.*param.g).*Rrk_tmp;
    
    drrk_dphi=0;
    
    dRrk_dV_tmp_1=dA0_dV+dA1_dV*(param.T./param.BWL);
    dRrk_dV_tmp_2=dA2_dV*((param.T-param.ZCBK).^3./param.DVK)+...
        dA3_dV*(param.DIVCAN./param.DVK);
    dRrk_dV_tmp_3=dRrk_dV_tmp_1+dRrk_dV_tmp_2;
    drrk_dV=(param.DVK*param.rho_w*param.g)*dRrk_dV_tmp_3;
  
end

function [drrhH_dV,drrhH_dphi,rrhH] = RrhH(V,phi,param)
    
    Fn = V/sqrt(param.g*param.LWL);
    if (Fn>0.55)
        Fn=0.55;
    end
    
    
    if (Fn>=0 && Fn<0.20)
        u0=0;
        u1=0;
        u2=0;
        u3=0;
        u4=0;
        u5=0;
        du0_dV=0;
        du1_dV=0;
        du2_dV=0;
        du3_dV=0;
        du4_dV=0;
        du5_dV=0;
    elseif (Fn>=0.20 && Fn <0.25)
        u0=0.000536*Fn-0.0001072;
        u1=-0.000028*Fn+0.0000056;
        u2=-0.000114*Fn+0.0000228;
        u3=0.000032*Fn-0.0000064;
        u4=-0.00014*Fn+0.000028;
        u5=-0.000034*Fn+0.0000068;
        du0_dV=(sqrt(param.g*param.LWL))*(0.000536);
        du1_dV=(sqrt(param.g*param.LWL))*(-0.000028);
        du2_dV=(sqrt(param.g*param.LWL))*(-0.000114);
        du3_dV=(sqrt(param.g*param.LWL))*(0.000032);
        du4_dV=(sqrt(param.g*param.LWL))*(-0.00014);
        du5_dV=(sqrt(param.g*param.LWL))*(-0.000034);
    elseif (Fn>=0.25 && Fn<0.30)
        u0=0.0007896*Fn-0.0001706;
        u1=-0.001236*Fn+0.0003076;
        u2=-0.001284*Fn+0.0003153;
        u3=0.000106*Fn-0.0000249;
        u4=0.001058*Fn-0.0002715;
        u5=0.000026*Fn-0.0000082;
        du0_dV=(sqrt(param.g*param.LWL))*(0.0007896);
        du1_dV=(sqrt(param.g*param.LWL))*(-0.001236);
        du2_dV=(sqrt(param.g*param.LWL))*(-0.001284);
        du3_dV=(sqrt(param.g*param.LWL))*(0.000106);
        du4_dV=(sqrt(param.g*param.LWL))*(0.001058);
        du5_dV=(sqrt(param.g*param.LWL))*(0.000026);
    elseif (Fn>=0.30 && Fn<0.35)
        u0=0.0315404*Fn-0.00939584;
        u1=-0.003024*Fn+0.000844;
        u2=-0.001882*Fn+0.0004947;
        u3=0.00026*Fn-0.0000711;
        u4=-0.001998*Fn+0.0006453;
        u5=-0.000528*Fn+0.000158;
        du0_dV=(sqrt(param.g*param.LWL))*(0.0315404);
        du1_dV=(sqrt(param.g*param.LWL))*(-0.003024);
        du2_dV=(sqrt(param.g*param.LWL))*(-0.001882);
        du3_dV=(sqrt(param.g*param.LWL))*(0.00026);
        du4_dV=(sqrt(param.g*param.LWL))*(-0.001998);
        du5_dV=(sqrt(param.g*param.LWL))*(-0.000528);
    elseif (Fn>=0.35 && Fn < 0.40)
        u0=-0.050184*Fn + 0.0192077;
        u1=0.00358*Fn-0.0014674;
        u2=0.007732*Fn-0.0028702;
        u3=-0.000022*Fn+0.0000276;
        u4=-0.01052*Fn+0.003628;
        u5=-0.00173*Fn+0.0005787;
        du0_dV=(sqrt(param.g*param.LWL))*(-0.050184);
        du1_dV=(sqrt(param.g*param.LWL))*(0.00358);
        du2_dV=(sqrt(param.g*param.LWL))*(0.007732);
        du3_dV=(sqrt(param.g*param.LWL))*(-0.000022);
        du4_dV=(sqrt(param.g*param.LWL))*(-0.01052);
        du5_dV=(sqrt(param.g*param.LWL))*(-0.00173);
    elseif (Fn>=0.40 && Fn<0.45)
        u0=-0.048112*Fn+0.0183789;
        u1=0.003452*Fn-0.0014162;
        u2=0.006642*Fn-0.0024342;
        u3=0.00016*Fn-0.0000452;
        u4=-0.008528*Fn+0.0028312;
        u5=0.006318*Fn-0.0026405;
        du0_dV=(sqrt(param.g*param.LWL))*(-0.048112);
        du1_dV=(sqrt(param.g*param.LWL))*(0.003452);
        du2_dV=(sqrt(param.g*param.LWL))*(0.006642);
        du3_dV=(sqrt(param.g*param.LWL))*(0.00016);
        du4_dV=(sqrt(param.g*param.LWL))*(-0.008528);
        du5_dV=(sqrt(param.g*param.LWL))*(0.006318);
    elseif (Fn>=0.45 && Fn < 0.50)
        u0=0.061478*Fn-0.0309366;
        u1=-0.005704*Fn+0.002704;
        u2=-0.02428*Fn+0.0114807;
        u3=0.003188*Fn-0.0014078;
        u4=0.00515*Fn-0.0033239;
        u5=-0.007348*Fn+0.0035092;
        du0_dV=(sqrt(param.g*param.LWL))*(0.061478);
        du1_dV=(sqrt(param.g*param.LWL))*(-0.005704);
        du2_dV=(sqrt(param.g*param.LWL))*(-0.02428);
        du3_dV=(sqrt(param.g*param.LWL))*(0.003188);
        du4_dV=(sqrt(param.g*param.LWL))*(0.00515);
        du5_dV=(sqrt(param.g*param.LWL))*(-0.007348);
    elseif (Fn>=0.50 && Fn <=0.55)
        u0=0.035698*Fn-0.0180466;
        u1=-0.004538*Fn+0.002121;
        u2=-0.001024*Fn-0.0001473;
        u3=0.000568*Fn-0.0000978;
        u4=0.005342*Fn-0.0034199;
        u5=0.000948*Fn-0.0006388;
        du0_dV=(sqrt(param.g*param.LWL))*(0.035698);
        du1_dV=(sqrt(param.g*param.LWL))*(-0.004538);
        du2_dV=(sqrt(param.g*param.LWL))*(-0.001024);
        du3_dV=(sqrt(param.g*param.LWL))*(0.000568);
        du4_dV=(sqrt(param.g*param.LWL))*(0.005342);
        du5_dV=(sqrt(param.g*param.LWL))*(0.000948);
    end
        
  
    RrhH_tmp_1=u0+u1*(param.LWL./param.BWL)+u2*(param.BWL./param.TCAN)+...
        u3*((param.BWL./param.TCAN).^2);
    RrhH_tmp_2=u4*(param.XFB./param.LWL)+u5*((param.XFB./param.LWL).^2);
    RrhH_tmp_3=RrhH_tmp_1+RrhH_tmp_2;
    
    
    RrhH20=(param.DIVCAN.*param.g.*param.rho_w) .*RrhH_tmp_3;
    
    rrhH = RrhH20* 6* (phi*pi/180)^1.7;
    drrhH_dphi=RrhH20* 1.7*6* (phi*pi/180)^0.7;
    
    
    dRrhH_dV_tmp_1=du0_dV+du1_dV*(param.LWL./param.BWL)+...
        du2_dV*(param.BWL./param.TCAN);
    dRrhH_dV_tmp_2=du3_dV*((param.BWL./param.TCAN).^2)+...
        du4_dV*(param.XFB./param.LWL)+du5_dV*((param.XFB./param.LWL).^2);
    dRrhH_dV_tmp_3=dRrhH_dV_tmp_1+dRrhH_dV_tmp_2;
    drrhH_dV=(param.DIVCAN.*param.g.*param.rho_w)*6*((phi*pi/180)^1.7)*...
        dRrhH_dV_tmp_3;
    
   
end

function [drrh_dV,drrh_dphi,rrh] = Rrh(V,param)

    Fn = V/sqrt(param.g*param.LWL);
    if Fn>0.60
        Fn=0.60;
    end
    if (Fn>=0 && Fn<0.10)
        a0=-0.014*Fn;
        a1=0.403*Fn;
        a2=0.47*Fn;
        a3=-0.227*Fn;
        a4=-0.119*Fn;
        a5=0.061*Fn;
        a6=-0.086*Fn;
        a7=-0.307*Fn;
        a8=-0.553*Fn;
        da0_dV=(sqrt(param.g*param.LWL))*-0.014;
        da1_dV=(sqrt(param.g*param.LWL))*0.403;
        da2_dV=(sqrt(param.g*param.LWL))*0.47;
        da3_dV=(sqrt(param.g*param.LWL))*(-0.227);
        da4_dV=(sqrt(param.g*param.LWL))*(-0.119);
        da5_dV=(sqrt(param.g*param.LWL))*0.061;
        da6_dV=(sqrt(param.g*param.LWL))*(-0.086);
        da7_dV=(sqrt(param.g*param.LWL))*(-0.307);
        da8_dV=(sqrt(param.g*param.LWL))*(-0.553);        
    elseif (Fn>=0.10 && Fn<0.15)
        a0=0.036*Fn-0.005;
        a1=-4.422*Fn+0.4825;
        a2=2.646*Fn-0.2176;
        a3=0.446*Fn-0.0673;
        a4=0.432*Fn-0.0551;
        a5=0.114*Fn-0.0053;
        a6=0.062*Fn-0.0148;
        a7=4.056*Fn-0.4363;
        a8=-2.35*Fn+0.1797;
        da0_dV=(sqrt(param.g*param.LWL))*0.036;
        da1_dV=(sqrt(param.g*param.LWL))*(-4.422);
        da2_dV=(sqrt(param.g*param.LWL))*(2.646);
        da3_dV=(sqrt(param.g*param.LWL))*(0.446);
        da4_dV=(sqrt(param.g*param.LWL))*0.432;
        da5_dV=(sqrt(param.g*param.LWL))*0.114;
        da6_dV=(sqrt(param.g*param.LWL))*0.062;
        da7_dV=(sqrt(param.g*param.LWL))*4.056;
        da8_dV=(sqrt(param.g*param.LWL))*(-2.35);
    elseif (Fn>=0.15 && Fn<0.20)
        a0=0.02*Fn-0.0026;
        a1=1.474*Fn-0.4019;
        a2=-2.312*Fn+0.5261;
        a3=0.188*Fn-0.0286;
        a4=0.112*Fn-0.0071;
        a5=-0.214*Fn+0.0439;
        a6=0.134*Fn-0.0256;
        a7=-1.4*Fn+0.3821;
        a8=2.16*Fn-0.4968;
        da0_dV=(sqrt(param.g*param.LWL))*0.02;
        da1_dV=(sqrt(param.g*param.LWL))*1.474;
        da2_dV=(sqrt(param.g*param.LWL))*(-2.312);
        da3_dV=(sqrt(param.g*param.LWL))*0.188;
        da4_dV=(sqrt(param.g*param.LWL))*0.112;
        da5_dV=(sqrt(param.g*param.LWL))*(-0.214);
        da6_dV=(sqrt(param.g*param.LWL))*0.134;
        da7_dV=(sqrt(param.g*param.LWL))*(-1.4);
        da8_dV=(sqrt(param.g*param.LWL))*2.16;
    elseif (Fn>=0.20 && Fn<0.25)
        a0=0.026*Fn-0.0038;
        a1=3.068*Fn-0.7207;
        a2=-3.8*Fn+0.8237;
        a3=0.12*Fn-0.015;
        a4=0.242*Fn-0.0331;
        a5=-0.62*Fn+0.1251;
        a6=0.196*Fn-0.038;
        a7=-3.232*Fn+0.7485;
        a8=3.736*Fn-0.812;
        da0_dV=(sqrt(param.g*param.LWL))*0.026;
        da1_dV=(sqrt(param.g*param.LWL))*3.068;
        da2_dV=(sqrt(param.g*param.LWL))*(-3.8);
        da3_dV=(sqrt(param.g*param.LWL))*0.12;
        da4_dV=(sqrt(param.g*param.LWL))*0.242;
        da5_dV=(sqrt(param.g*param.LWL))*(-0.62);
        da6_dV=(sqrt(param.g*param.LWL))*0.196;
        da7_dV=(sqrt(param.g*param.LWL))*(-3.232);
        da8_dV=(sqrt(param.g*param.LWL))*3.736;
    elseif (Fn>=0.25 && Fn<0.30)
        a0=0.058*Fn-0.0118;
        a1=-16.936*Fn+4.2803;
        a2=12.308*Fn-3.2033;
        a3=0.238*Fn-0.0445;
        a4=0.49*Fn-0.0951;
        a5=-0.028*Fn-0.0229;
        a6=0.364*Fn-0.08;
        a7=15.818*Fn-4.014;
        a8=-9.678*Fn+2.5415;
        da0_dV=(sqrt(param.g*param.LWL))*0.058;
        da1_dV=(sqrt(param.g*param.LWL))*(-16.936);
        da2_dV=(sqrt(param.g*param.LWL))*12.308;
        da3_dV=(sqrt(param.g*param.LWL))*0.238;
        da4_dV=(sqrt(param.g*param.LWL))*0.49;
        da5_dV=(sqrt(param.g*param.LWL))*(-0.028);
        da6_dV=(sqrt(param.g*param.LWL))*0.364;
        da7_dV=(sqrt(param.g*param.LWL))*15.818;
        da8_dV=(sqrt(param.g*param.LWL))*(-9.678);
    elseif (Fn>=0.30 && Fn < 0.35)
        a0=-0.048*Fn+0.02;
        a1=13.988*Fn-4.9969;
        a2=-11.408*Fn+3.9115;
        a3=-1.302*Fn+0.4175;
        a4=-0.398*Fn+0.1713;
        a5=-2.336*Fn+0.6695;
        a6=1.09*Fn-0.2978;
        a7=-14.182*Fn+4.986;
        a8=10.412*Fn-3.4855;
        da0_dV=(sqrt(param.g*param.LWL))*(-0.048);
        da1_dV=(sqrt(param.g*param.LWL))*13.988;
        da2_dV=(sqrt(param.g*param.LWL))*(-11.408);
        da3_dV=(sqrt(param.g*param.LWL))*(-1.302);
        da4_dV=(sqrt(param.g*param.LWL))*(-0.398);
        da5_dV=(sqrt(param.g*param.LWL))*(-2.336);
        da6_dV=(sqrt(param.g*param.LWL))*1.09;
        da7_dV=(sqrt(param.g*param.LWL))*(-14.182);
        da8_dV=(sqrt(param.g*param.LWL))*10.412;
    elseif (Fn>=0.35 && Fn<0.40)
        a0=-0.192*Fn+0.0704;
        a1=48.212*Fn-16.9753;
        a2=-28.678*Fn+9.956;
        a3=2.266*Fn-0.8313;
        a4=-2.356*Fn+0.8566;
        a5=-7.736*Fn+2.5595;
        a6=1.756*Fn-0.5309;
        a7=-49.546*Fn+17.3634;
        a8=20.556*Fn-7.0359;
        da0_dV=(sqrt(param.g*param.LWL))*(-0.192);
        da1_dV=(sqrt(param.g*param.LWL))*48.212;
        da2_dV=(sqrt(param.g*param.LWL))*(-28.678);
        da3_dV=(sqrt(param.g*param.LWL))*2.266;
        da4_dV=(sqrt(param.g*param.LWL))*(-2.356);
        da5_dV=(sqrt(param.g*param.LWL))*(-7.736);
        da6_dV=(sqrt(param.g*param.LWL))*1.756;
        da7_dV=(sqrt(param.g*param.LWL))*(-49.546);
        da8_dV=(sqrt(param.g*param.LWL))*20.556;
    elseif (Fn>=0.40 && Fn<0.45)
        a0=-0.214*Fn+0.0792;
        a1=21.844*Fn-6.4281;
        a2=-9.42*Fn+2.2528;
        a3=4.982*Fn-1.9177;
        a4=-1.184*Fn+0.3878;
        a5=-5.388*Fn+1.6203;
        a6=2.474*Fn-0.8181;
        a7=-21.468*Fn+6.1322;
        a8=3.42*Fn-0.1815;
        da0_dV=(sqrt(param.g*param.LWL))*(-0.214);
        da1_dV=(sqrt(param.g*param.LWL))*21.844;
        da2_dV=(sqrt(param.g*param.LWL))*(-9.42);
        da3_dV=(sqrt(param.g*param.LWL))*4.982;
        da4_dV=(sqrt(param.g*param.LWL))*(-1.184);
        da5_dV=(sqrt(param.g*param.LWL))*(-5.388);
        da6_dV=(sqrt(param.g*param.LWL))*2.474;
        da7_dV=(sqrt(param.g*param.LWL))*(-21.468);
        da8_dV=(sqrt(param.g*param.LWL))*3.42;
    elseif (Fn>=0.45 && Fn<0.50)
        a0=-0.06*Fn+0.0099;
        a1=75.118*Fn-30.4014;
        a2=-86.884*Fn+37.1116;
        a3=5.174*Fn-2.0041;
        a4=6.16*Fn-2.917;
        a5=8.154*Fn-4.4736;
        a6=4.142*Fn-1.5687;
        a7=-72.59*Fn+29.1371;
        a8=77.918*Fn-33.7056;
        da0_dV=(sqrt(param.g*param.LWL))*(-0.06);
        da1_dV=(sqrt(param.g*param.LWL))*75.118;
        da2_dV=(sqrt(param.g*param.LWL))*(-86.884);
        da3_dV=(sqrt(param.g*param.LWL))*5.174;
        da4_dV=(sqrt(param.g*param.LWL))*6.16;
        da5_dV=(sqrt(param.g*param.LWL))*8.154;
        da6_dV=(sqrt(param.g*param.LWL))*4.142;
        da7_dV=(sqrt(param.g*param.LWL))*(-72.59);
        da8_dV=(sqrt(param.g*param.LWL))*77.918;
    elseif (Fn>=0.50 && Fn<0.55)
        a0=1.392*Fn-0.7161;
        a1=-111.916*Fn+63.1156;
        a2=5.286*Fn-8.9734;
        a3=5.624*Fn-2.2291;
        a4=20.144*Fn-9.909;
        a5=43.152*Fn-21.9726;
        a6=8.306*Fn-3.6507;
        a7=100.776*Fn-57.5459;
        a8=3.494*Fn+3.5064;
        da0_dV=(sqrt(param.g*param.LWL))*1.392;
        da1_dV=(sqrt(param.g*param.LWL))*(-111.916);
        da2_dV=(sqrt(param.g*param.LWL))*5.286;
        da3_dV=(sqrt(param.g*param.LWL))*5.624;
        da4_dV=(sqrt(param.g*param.LWL))*20.144;
        da5_dV=(sqrt(param.g*param.LWL))*43.152;
        da6_dV=(sqrt(param.g*param.LWL))*8.306;
        da7_dV=(sqrt(param.g*param.LWL))*100.776;
        da8_dV=(sqrt(param.g*param.LWL))*3.494;
    elseif (Fn>=0.55 && Fn<=0.60)
        a0=0.626*Fn-0.2948;
        a1=-137.702*Fn+77.2979;
        a2=98.296*Fn-60.1289;
        a3=2.044*Fn-0.2601;
        a4=8.764*Fn-3.65;
        a5=19.698*Fn-9.0729;
        a6=-1.37*Fn-1.6711;
        a7=136.64*Fn-77.2711;
        a8=-86.384*Fn+52.9393;
        da0_dV=(sqrt(param.g*param.LWL))*0.626;
        da1_dV=(sqrt(param.g*param.LWL))*(-137.702);
        da2_dV=(sqrt(param.g*param.LWL))*98.296;
        da3_dV=(sqrt(param.g*param.LWL))*2.044;
        da4_dV=(sqrt(param.g*param.LWL))*8.764;
        da5_dV=(sqrt(param.g*param.LWL))*19.698;
        da6_dV=(sqrt(param.g*param.LWL))*(-1.37);
        da7_dV=(sqrt(param.g*param.LWL))*136.64;
        da8_dV=(sqrt(param.g*param.LWL))*(-86.384);
    end
    
    
    SC=(1.97+0.171*param.BWL/param.TCAN)*sqrt(param.DIVCAN*param.LWL)*...
        (0.65/param.CMS)^(1/3);
    
    Rrh_tmp_1=param.DIVCAN.*param.g.*param.rho_w;
    Rrh_tmp_2=a1*(param.XFB./param.LWL)+a2*(param.CPL)+...
        a3*(param.DIVCAN.^(2/3)./param.AW);
    Rrh_tmp_3=a4*(param.BWL./param.LWL)+a5*(param.DIVCAN.^(2/3)./SC)+...
        a6*(param.XFB./param.XFF);
    Rrh_tmp_4=a7*((param.XFB./param.LWL).^2)+a8*(param.CPL.^2);
    Rrh_tmp_5=a0+(param.DIVCAN^(1/3)/param.LWL)*(Rrh_tmp_2+Rrh_tmp_3+Rrh_tmp_4);
    rrh=Rrh_tmp_1*Rrh_tmp_5;
  
    drrh_dphi=0;
    
    drrh_dV_tmp_1=param.DIVCAN.*param.g.*param.rho_w;
    drrh_dV_tmp_2=da1_dV*(param.XFB./param.LWL)+da2_dV*(param.CPL)+da3_dV*...
        (param.DIVCAN.^(2/3)./param.AW);
    drrh_dV_tmp_3=da4_dV*(param.BWL./param.LWL)+...
        da5_dV*(param.DIVCAN.^(2/3)./SC)+da6_dV*(param.XFB./param.XFF);
    drrh_dV_tmp_4=da7_dV*((param.XFB./param.LWL).^2)+da8_dV*(param.CPL.^2);
    drrh_dV_tmp_5=da0_dV+(param.DIVCAN^(1/3)/param.LWL)*(drrh_dV_tmp_2+...
        drrh_dV_tmp_3+drrh_dV_tmp_4);
    drrh_dV=drrh_dV_tmp_1*drrh_dV_tmp_5;
   
    
end

function [dri_dV,dri_dphi,ri]= Ri(V,phi,param,f_heel,dFheel_dphi,dFheel_dV)

    if phi>30
        phi=30;
    end
    if (phi>=0 && phi<10)
        A1=0.07437*phi+3.7455;
        dA1_dphi=0.07437;
        A2=-0.12208*phi-3.6246;
        dA2_dphi=-0.12208;
        A3=-2.95e-3*phi+0.0589;
        dA3_dphi=-2.95e-3;
        A4=1.2e-3*phi-0.0296;
        dA4_dphi=1.2e-3;
        B0=0.01925*phi+1.2306;
        dB0_dphi=0.01925;
        B1=-0.05715*phi-0.7256;
        dB1_dphi=-0.05715;
    elseif (phi>=10 && phi< 20)
        A1=-0.053*phi+5.0192;
        dA1_dphi=-0.053;
        A2=0.0865*phi-5.7104;
        dA2_dphi=0.0865;
        A3=-1.1e-4*phi+0.0305;
        dA3_dphi=-1.1e-4;
        A4=-1.97e-3*phi+0.0319;
        dA4_dphi=-1.97e-3;
        B0=0.01219*phi+1.3012;
        dB0_dphi=0.01219;
        B1=-0.02651*phi-1.032;
        dB1_dphi=-0.02651;
    elseif (phi>=20 && phi<=30)
        A1=-0.04701*phi+4.8994;
        dA1_dphi=-0.04701;
        A2=0.10227*phi-6.0258;
        dA2_dphi=0.10227;
        A3=-3.3e-4*phi+0.0349;
        dA3_dphi=-3.3e-4;
        A4=-1.97e-3*phi+0.0319;
        dA4_dphi=-1.97e-3;
        B0=-7.06e-3*phi+1.6862;
        dB0_dphi=-7.06e-3;
        B1=0.02123*phi-1.9868;
        dB1_dphi=0.02123;
    end
    

    Fn = V./sqrt(param.g*param.LWL);
    Te_tmp_1=A1*(param.TCAN/param.T)+A2*((param.TCAN./param.T).^2);
    Te_tmp_2=A3*(param.BWL./param.TCAN)+A4*(param.CHTPK./param.CHRTK);
    Te_tmp_3=B0+B1*Fn;
    Te=param.T*(Te_tmp_1+Te_tmp_2)*Te_tmp_3;    
    ri = f_heel^2/ (pi* Te^2* 0.5* param.rho_w* V^2);
    
    Te_tmp_4=dA1_dphi*(param.TCAN/param.T)+dA2_dphi*((param.TCAN./param.T).^2);
    Te_tmp_5=dA3_dphi*(param.BWL./param.TCAN)+dA4_dphi*(param.CHTPK./param.CHRTK);
    Te_tmp_6=dB0_dphi+Fn*dB1_dphi;
    
    
    
    dTe_dphi=Te_tmp_3*param.T*(Te_tmp_4+Te_tmp_5)+...
        param.T*(Te_tmp_1+Te_tmp_2)*Te_tmp_6;
    dri_dphi=(2/(pi*Te^2*0.5*param.rho_w*V^2))*dFheel_dphi-...
        2*f_heel^2/(pi*0.5*param.rho_w*V^2)*Te^(-3)*dTe_dphi;
    
    
    dFn_dV=1/(sqrt(param.g*param.LWL));
    dTe_dV=param.T*(Te_tmp_1+Te_tmp_2)*B1*dFn_dV;
    dTe2_dV=-2*(1/Te^3)*dTe_dV;
    dV2_dV=-2*(1/V^3);
    dTeV_dV=(1/Te^2)*dV2_dV+(1/V^2)*dTe2_dV;
    dFheel2_dV=2*f_heel*dFheel_dV;
    dri_dV=(f_heel^2/(pi*0.5*param.rho_w))*dTeV_dV+...
        (1/(pi*Te^2*0.5*param.rho_w*V^2))*dFheel2_dV;
  
end

function [B3,C3,m_right] = M_right(phi,param) 

    M1=(param.KM - param.KG)*sin(phi*pi/180)*param.g*param.rho_w*...
        (param.DIVCAN+param.DVK);
    M2 = param.MMVBLCRW* param.g* param.b* cos(phi*pi/180);
    m_right = M1+M2;
    
    dM1_dphi=(param.KM-param.KG)*cos(phi*pi/180)*param.g*param.rho_w*...
        (param.DIVCAN+param.DVK);
    dM2_dphi=-param.MMVBLCRW* param.g* param.b* sin(phi*pi/180);
    C3=dM1_dphi+dM2_dphi;

    B3=0;
    
end

function [B2,B4,C2,C4,m_heel,f_heel,dFheel_dphi,dFheel_dV,f_drive]=...
    aero(V,phi,V_tw,alfa_tw,param)

    V1 = V + V_tw* cos(alfa_tw*pi/180);
    V2 = V_tw* sin(alfa_tw*pi/180)* cos(phi*pi/180);

    V_eff = sqrt( V1.^2 + V2.^2 );
    alfa_eff = (atan2(V2,V1))*180/pi;

    AM = 0.5* param.P* param.E*param.MROACH;             
    AJ = 0.5* sqrt( param.I^2 + param.J^2)* param.LPG;   
    AS = 1.15* param.SL* param.J;                         
    AF = 0.5* param.I* param.J;                           
    AN = AF + AM;                                           

    ZCEM = 0.39* param.P + param.BAD;
    ZCEJ = 0.39* param.I;
    ZCES = 0.59* param.I;
    
     if (alfa_eff >0 && alfa_eff <27)
         Cl_M=0.0638888*alfa_eff;
         m0=0.0638888;
         Cl_J=0.0555555*alfa_eff;
         j0=0.0555555;
         Cl_S=0;
         s0=0;
         Cdp_M=7.4074e-4*alfa_eff;
         m2=7.4074e-4;
         Cdp_J=7.4074e-4*alfa_eff;
         j2=7.4074e-4;
         Cdp_S=0;
         s2=0;
     elseif (alfa_eff >= 27 && alfa_eff<50)
         Cl_M=-9.782e-3*alfa_eff+1.9891304;
         m0=-9.782e-3;
         Cl_J=-0.043478*alfa_eff+2.673913;
         j0=-0.043478;
         Cl_S=0.0652173*alfa_eff-1.760869;
         s0=0.0652173;
         Cdp_M=5.6521e-3*alfa_eff-0.132506;
         m2=5.6521e-3;
         Cdp_J=0.01*alfa_eff-0.25;
         j2=0.01;
         Cdp_S=0.0108695*alfa_eff-0.293478;
         s2=0.0108695;
     elseif (alfa_eff>=50 && alfa_eff <80)
         Cl_M=-0.018333*alfa_eff+2.416666;
         m0=-0.018333;
         Cl_J=-6.666e-3*alfa_eff+0.8333333;
         j0=-6.666e-3;
         Cl_S=-0.016666*alfa_eff+2.3333333;
         s0=-0.016666;
         Cdp_M=0.0216666*alfa_eff-0.933333;
         m2=0.0216666;
         Cdp_J=-3.333e-3*alfa_eff+0.4166666;
         j2=-3.333e-3;
         Cdp_S=0.0216666*alfa_eff-0.833333;
         s2=0.0216666;
     elseif (alfa_eff>=80 && alfa_eff<100)
         Cl_M=-5e-3*alfa_eff+1.35;
         m0=-5e-3;
         Cl_J=-0.015*alfa_eff+1.5;
         j0=-0.015;
         Cl_S=-7.5e-3*alfa_eff+1.6;
         s0=-7.5e-3;
         Cdp_M=0.01*alfa_eff;
         m2=0.01;
         Cdp_J=-7.5e-3*alfa_eff+0.75;
         j2=-7.5e-3;
         Cdp_S=0.015*alfa_eff-0.3;
         s2=0.015;
     elseif (alfa_eff>=100 && alfa_eff<=180)
         Cl_M=-0.010625*alfa_eff+1.9125;
         m0=-0.010625;
         Cl_J=0;
         j0=0;
         Cl_S=-0.010625*alfa_eff+1.9125;
         s0=-0.010625;
         Cdp_M=-1.25e-3*alfa_eff+1.125;
         m2=-1.25e-3;
         Cdp_J=0;
         j2=0;
         Cdp_S=-6.75e-3*alfa_eff+1.875; 
         s2=-6.75e-3;
     end         
                 
    switch param.SAILSET
        case 1    %Main only
            kM=1;
            kS=0;
            kJ=0;            
        case 3    %Main and Jib
            kM=1;
            kS=0;
            kJ=1;            
        case 5    %Main and Spi
            kM=1;
            kS=1;
            kJ=0; 
        case 7    %Main, Jib and Spi
            kM=1;
            kS=1;
            kJ=1; 
    end

    Cl = param.F*(kM*Cl_M*AM+kJ*Cl_J*AJ+kS*Cl_S*AS)/AN;
    Cdp=(kM*Cdp_M*AM + kJ*Cdp_J* AJ + kS*Cdp_S* AS )/ AN;
    ZCE = (kM*ZCEM* AM + kJ*ZCEJ* AJ + kS*ZCES* AS)/ (kM*AM + kJ*AJ + kS*AS);

    if alfa_eff < 45
        AR = (1.1*( param.EHM + param.AVGFREB)).^2./AN;
    else
        AR = (1.1*( param.EHM )).^2./AN;
    end

    CdI = Cl^2* ( 1/(pi*AR) + 0.005); %induced resistance
    Cd0 = 1.13* ( (param.B*param.AVGFREB) + (param.EHM*param.EMDC) )/ AN;

    Cd = Cdp + Cd0 + CdI;

    L = 0.5 * param.rho_a * V_eff.^2* AN* Cl;
    D = 0.5 * param.rho_a * V_eff.^2* AN* Cd;

    f_drive = L * sin(alfa_eff*pi/180) - D * cos(alfa_eff*pi/180);   
    f_heel = L * cos(alfa_eff*pi/180) + D * sin(alfa_eff*pi/180);
    m_heel = f_heel*(ZCE + param.T - param.ZCBK);  

    
    
    
    dV2_dphi=(-V_tw*sin(alfa_tw*pi/180)*sin(phi*pi/180));
    dalfaeff_dphi=(V1/V_eff^2)*dV2_dphi;
    dCdpJ_dphi=j2*dalfaeff_dphi;
    dCdpS_dphi=s2*dalfaeff_dphi;
    dCdpM_dphi=m2*dalfaeff_dphi;
    dClj_dphi=j0*dalfaeff_dphi;
    dCls_dphi=s0*dalfaeff_dphi;
    dClM_dphi=m0*dalfaeff_dphi;
    dCl_dphi=kM*AM*dClM_dphi+kS*AS*dCls_dphi+kJ*AJ*dClj_dphi;
    dCdp_dphi=(1/AN)*(kM*AM*dCdpM_dphi+kS*AS*dCdpS_dphi+kJ*AJ*dCdpJ_dphi);
    dCd0_dphi=0;
    dCdI_dphi=2*Cl*((1/(pi*AR))+0.005)*dCl_dphi;
    dCd_dphi=dCdp_dphi+dCd0_dphi+dCdI_dphi;
    dVeff_dphi=V2*(V1^2+V2^2)*(-V_tw*sin(alfa_tw*pi/180)*sin(phi*pi/180));
    dD_dphi=0.5*param.rho_a*AN*(V_eff^2*dCd_dphi+2*V_eff*dVeff_dphi);
    dcosalfaeff_dphi=-sin(alfa_eff*pi/180)*(V1/V_eff^2)*...
        (-V_tw*sin(alfa_tw*pi/180)*sin(phi*pi/180));
    dsinalfaeff_dphi=cos(alfa_eff*pi/180)*(V1/V_eff^2)*...
        (-V_tw*sin(alfa_tw*pi/180)*sin(phi*pi/180));
  
    
    dV2Cl_dphi=2*V_eff*Cl*dVeff_dphi+V_eff^2*dCl_dphi;
    dL_dphi=0.5*param.rho_a*AN*dV2Cl_dphi;
    dDcosalpha_dphi=cos(alfa_eff*pi/180)*dD_dphi+D*dcosalfaeff_dphi;
    dLsinalpha_dphi=sin(alfa_eff*pi/180)*dL_dphi+L*dsinalfaeff_dphi;
    C2=dLsinalpha_dphi-dDcosalpha_dphi;
    
    dDsinalfaeff_dphi=sin(alfa_eff*pi/180)*dD_dphi+D*dsinalfaeff_dphi;
    dLcosalfaeff_dphi=cos(alfa_eff*pi/180)*dL_dphi+L*dcosalfaeff_dphi;
    dFheel_dphi=dLcosalfaeff_dphi+dDsinalfaeff_dphi;
    C4=(ZCE+param.T-param.ZCBK)*dFheel_dphi;
    
    dalfaeff_dV=-(V2/(V1^2+V2^2));
    dCdpJ_dV=j2*dalfaeff_dV;
    dCdpS_dV=s2*dalfaeff_dV;
    dCdpM_dV=m2*dalfaeff_dV;
    dCdp_dV=(1/AN)*(kM*AM*dCdpM_dV+kS*AS*dCdpS_dV+kJ*AJ*dCdpJ_dV);
    dCd0_dV=0;
    dClS_dV=s0*dalfaeff_dV;
    dClJ_dV=j0*dalfaeff_dV;
    dClM_dV=m0*dalfaeff_dV;
    dCl_dV=(param.F/AN)*(kM*AM*dClM_dV+kS*AS*dClS_dV+kJ*AJ*dClJ_dV);
    dCdI_dV=2*Cl*((1/(pi*AR))+0.005)*dCl_dV;
    dCd_dV=dCdp_dV+dCd0_dV+dCdI_dV;
    dV22_dV=0;
    dV12_dV=2*V1;
    dVeff_dV=0.5*((V1^2+V2^2)^(0.5))*(dV12_dV+dV22_dV);
    dVeff2Cd_dV=2*V_eff*Cd*dVeff_dV+V_eff^2*dCd_dV;
    dD_dV=0.5*param.rho_a*AN*dVeff2Cd_dV;
    dcosalfaeff_dV=-sin(alfa_eff*pi/180)*dalfaeff_dV;
    dDcosalfaeff_dV=D*dcosalfaeff_dV+cos(alfa_eff*pi/180)*dD_dV;
    dsinalfaeff_dV=cos(alfa_eff*pi/180)*dalfaeff_dV;
    dVeff2Cl_dV=2*V_eff*Cl*dVeff_dV+V_eff^2*dCl_dV;
    dL_dV=0.5*param.rho_a*AN*dVeff2Cl_dV;
    dLsinalfaeff_dV=L*dsinalfaeff_dV+sin(alfa_eff*pi/180)*dL_dV;
    B2=dLsinalfaeff_dV-dDcosalfaeff_dV;
    
    dDsinalfaeff_dV=D*dsinalfaeff_dV+sin(alfa_eff*pi/180)*dD_dV;
    dLcosalfaeff_dV=L*dcosalfaeff_dV+cos(alfa_eff*pi/180)*dL_dV;
    dFheel_dV=dLcosalfaeff_dV+dDsinalfaeff_dV;
    B4=(ZCE+param.T-param.ZCBK)*dFheel_dV;
    
  
end

function [param] =load_data()


    fid = fopen('param.txt'); 
    fseek(fid,0,'eof');
    endstatus = ftell(fid);
    frewind(fid);
    position = ftell(fid);
    while position < endstatus
        position = ftell(fid);
        entry = fscanf(fid,'%s  ',[1 1]); 

        switch entry
            case 'DIVCAN',      % Displ vol canoe body [m^3]
                                value = fscanf(fid,'%f\n',[1 1]);
                                param.DIVCAN     = value; 
            case 'LWL',         % Length at waterline [m]
                                value = fscanf(fid,'%f\n',[1 1]);
                                param.LWL        = value; 
            case 'BWL',         % Beam at waterline [m]
                                value = fscanf(fid,'%f\n',[1 1]);
                                param.BWL        = value; 
            case 'B',           % Beam (max) [m]
                                value = fscanf(fid,'%f\n',[1 1]);
                                param.B          = value; 
            case 'AVGFREB',     % Avg freeboard [m]
                                value = fscanf(fid,'%f\n',[1 1]);
                                param.AVGFREB    = value; 
            case 'XFB',         % Long ctr buoyancy [m]
                                value = fscanf(fid,'%f\n',[1 1]);
                                param.XFB        = value;  
            case 'XFF',         % Long ctf floatation [m]
                                value = fscanf(fid,'%f\n',[1 1]);
                                param.XFF        = value; 
            case 'CPL',         % Long prismatic coef [-]
                                value = fscanf(fid,'%f\n',[1 1]);
                                param.CPL        = value; 
            case 'HULLFF',      % Hull form factor [-] 
                                value = fscanf(fid,'%f\n',[1 1]);
                                param.HULLFF     = value; 
            case 'AW',          % Area of water plane [m^2]
                                value = fscanf(fid,'%f\n',[1 1]);
                                param.AW         = value; 
            case 'CMS',         % Midship section coef [-]
                                value = fscanf(fid,'%f\n',[1 1]);
                                param.CMS        = value; 
            case 'T',           % Draft of hull [m]
                                value = fscanf(fid,'%f\n',[1 1]);
                                param.T          = value; 
            case 'TCAN',        % Draft of canoe body [m]
                                value = fscanf(fid,'%f\n',[1 1]);
                                param.TCAN       = value; 
            case 'ALT',         % Total lateral area [m^2]
                                value = fscanf(fid,'%f\n',[1 1]);
                                param.ALT        = value; 
            case 'KG',          % Center of gravity above keel [m]
                                value = fscanf(fid,'%f\n',[1 1]);
                                param.KG         = value; 
            case 'KM',          % Transv metacenter above keel [m]
                                value = fscanf(fid,'%f\n',[1 1]);
                                param.KM         = value; 
            case 'DVK',         % Displ vol of keel [m^3]
                                value = fscanf(fid,'%f\n',[1 1]);
                                param.DVK        = value; 
            case 'APK',         value = fscanf(fid,'%f\n',[1 1]);
                                param.APK        = value;
            case 'ASK',         % Area of skeg [m^2]
                                value = fscanf(fid,'%f\n',[1 1]);
                                param.ASK        = value; 
            case 'SK',          % Wet surf area of keel [m^2]
                                value = fscanf(fid,'%f\n',[1 1]);
                                param.SK         = value; 
            case 'ZCBK',        value = fscanf(fid,'%f\n',[1 1]);
                                param.ZCBK       = value;
            case 'CHMEK',       % Mean chord length of keel [m]
                                value = fscanf(fid,'%f\n',[1 1]);
                                param.CHMEK      = value; 
            case 'CHRTK',       % Root chord length of keel [m]
                                value = fscanf(fid,'%f\n',[1 1]);
                                param.CHRTK      = value; 
            case 'CHTPK',       % Tip chord length of keel [m]
                                value = fscanf(fid,'%f\n',[1 1]);
                                param.CHTPK      = value; 
            case 'KEELFF',      % Keel form factor [-]
                                value = fscanf(fid,'%f\n',[1 1]);
                                param.KEELFF     = value; 
            case 'DELTTK',      % Thickness ratio of keel [-]
                                value = fscanf(fid,'%f\n',[1 1]);
                                param.DELTTK     = value; 
            case 'TAK',         % Draft @ aft perp of keel [m]
                                value = fscanf(fid,'%f\n',[1 1]);
                                param.TAK        = value; 
            case 'DVR',         % Displ vol of rudder [m^3]
                                value = fscanf(fid,'%f\n',[1 1]);
                                param.DVR        = value; 
            case 'APR',         value = fscanf(fid,'%f\n',[1 1]);
                                param.APR        = value;
            case 'SR',          % Wet surf area of rudder [m^2]
                                value = fscanf(fid,'%f\n',[1 1]);
                                param.SR         = value; 
            case 'CHMER',       % Mean chord length of rudder [m]
                                value = fscanf(fid,'%f\n',[1 1]);
                                param.CHMER      = value; 
            case 'CHRTR',       % Root chord length of rudder [m]
                                value = fscanf(fid,'%f\n',[1 1]);
                                param.CHRTR      = value; 
            case 'CHTPR',       % Tip chord length of rudder [m]
                                value = fscanf(fid,'%f\n',[1 1]);
                                param.CHTPR      = value; 
            case 'DELTTR',      % Thickness ratio of rudder [-]
                                value = fscanf(fid,'%f\n',[1 1]);
                                param.DELTTR     = value; 
            case 'RUDDFF',      % Rudder form factor [-]
                                value = fscanf(fid,'%f\n',[1 1]);
                                param.RUDDFF     = value; 
            case 'SAILSET',     % 3=main+jib; 5=main+spin; 7=main+jib+spin
                                value = fscanf(fid,'%f\n',[1 1]);
                                param.SAILSET    = value; 
            case 'P',           % Mainsail height [m]
                                value = fscanf(fid,'%f\n',[1 1]);
                                param.P          = value; 
            case 'E',           % Mainsail base [m]
                                value = fscanf(fid,'%f\n',[1 1]);
                                param.E          = value; 
            case 'MROACH',	    % Correction for mainsail roach [-]
                                value = fscanf(fid,'%f\n',[1 1]);
                                param.MROACH     = value; 
            case 'MFLB',        % Full main battens in main [0/1]
                                value = fscanf(fid,'%f\n',[1 1]);
                                param.MFLB       = value; 
            case 'BAD',         % Height of main boom above sheer [m]
                                value = fscanf(fid,'%f\n',[1 1]);
                                param.BAD        = value; 
            case 'I',           % Fore-triangle height [m]
                                value = fscanf(fid,'%f\n',[1 1]);
                                param.I          = value; 
            case 'J',           % Fore-triangle base [m]
                                value = fscanf(fid,'%f\n',[1 1]);
                                param.J          = value; 
            case 'LPG',         % Length perpendicular of jib [m]
                                value = fscanf(fid,'%f\n',[1 1]);
                                param.LPG        = value; 
            case 'SL',          % Spinnaker leech length [m]
                                value = fscanf(fid,'%f\n',[1 1]);
                                param.SL         = value; 
            case 'EHM',         % Mast height above sheer [m]
                                value = fscanf(fid,'%f\n',[1 1]);
                                param.EHM        = value; 
            case 'EMDC',        % Average mast diameter [m]
                                value = fscanf(fid,'%f\n',[1 1]);
                                param.EMDC       = value; 
            case 'MMVBLCRW',    % Mass of moveable crew [kg]
                                value = fscanf(fid,'%f\n',[1 1]);
                                param.MMVBLCRW   = value; 
            case 'rho_w',       value = fscanf(fid,'%f\n',[1 1]);
                                param.rho_w  = value;
            case 'ni_w',        value = fscanf(fid,'%f\n',[1 1]);
                                param.ni_w   = value;
            case 'rho_a',       value = fscanf(fid,'%f\n',[1 1]);
                                param.rho_a  = value;
            case 'g',           value = fscanf(fid,'%f\n',[1 1]);
                                param.g      = value;
            case 'b',           value = fscanf(fid,'%f\n',[1 1]);
                                param.b=value;
            case 'F',           value = fscanf(fid,'%f\n',[1 1]);
                                param.F=value;
            case 'V_tw',
                                value = fscanf(fid,'%f\n',[1 inf]);
                                if  length(value) == 3 && value(2) < value(1)  
                                    param.V_tw = value(1):value(2):value(3);
                                else
                                    param.V_tw = value;
                                end
            case 'alfa_tw'
                                value = fscanf(fid,'%f\n',[1 inf]);
                                if      length(value) == 3 && value(2) < value(1) 
                                    param.alfa_tw = value(1):value(2):value(3);
                                else
                                    param.alfa_tw = value;
                                end
        end
    end
    fclose(fid);
end

function testing_data(param)

    if param.LWL/param.BWL  <= 2.73 || param.LWL/param.BWL >=  5.00                   
        Warning('LWL/BWL  is out of valuable interval')        
        pause
    end
    if param.BWL/param.TCAN <= 2.46 || param.LWL/param.BWL >= 19.38                   
        Warning('BWL/TCAN is out of valuable interval')          
        pause 
    end
    if param.LWL/param.DIVCAN^(1/3) <= 4.34 || param.LWL/param.DIVCAN^(1/3) >= 8.50
        Warning('LWL/DIVCAN^(1/3) is out of valuable interval')
        pause
    end
    %if geom.XFB < 0 
    %if geom.XFF
    if param.CPL <= 0.52 || param.CPL >= 0.6
        Warning('CPL is out of valuable interval')
        pause
    end
    if param.CMS <= 0.65 || param.CMS >= 0.78
        Warning('CMS is out of valuable interval')
        pause
    end
    if param.AW/param.DIVCAN^(2/3) <= 3.78 || param.AW/param.DIVCAN^(2/3) >= 12.67
        Warning('Loading Factor is out of valuable data')
        pause
    end

end

function [] = postprocessing(x,param)

  
    k = 0;
    for i = 1:length(param.V_tw)
        for j = 1:length(param.alfa_tw)

            k = k+1;
  
            V_tw=x(1,(i-1)*length(param.V_tw)+j);
            alfa_tw=x(2,(i-1)*length(param.V_tw)+j);
            V=x(3,(i-1)*length(param.V_tw)+j);
            phi=x(4,(i-1)*length(param.V_tw)+j);
            
            VMG = V * cos(alfa_tw*pi/180);

            if isnan(V)
                [f_side,f_drive,m_heel,ri,rrh,rrhH,rrk,rrkH,rvh,rvhH,rvk,rvr]=...
                    deal(NaN);
            else
                [~,~,~,~,m_heel,f_heel,dFheel_dphi,dFheel_dV,f_drive]=...
                    aero(V,phi,V_tw,alfa_tw,param);
                [~,~,ri]= Ri(V,phi,param,f_heel,dFheel_dphi,dFheel_dV);
                [~,~,rrh]= Rrh(V,param);
                [~,~,rrhH]= RrhH(V,phi,param);
                [~,~,rrk]= Rrk(V,param);
                [~,~,rrkH]= RrkH(V,phi,param);
                [~,~,rvh]= Rvh(V,param);
                [~,~,rvhH]= RvhH(V,phi,param);
                [~,~,rvk]=Rvk(V,param);
                [~,~,rvr]=Rvr(V,param);
               
                f_side = f_heel*cos(phi*pi/180);
               
                

            end
            r_tot = ri + rrh + rrhH + rrk + rrkH + rvh + rvhH + rvk + rvr;

            Rall(k,:) = real([rrh+rrhH rvh+rvhH rrk+rrkH rvk ri rvr]);
            Rperc(k,:) = Rall(k,:)./r_tot;
            data(k,:) = real([V_tw alfa_tw V VMG phi param.b param.F f_drive...
                f_side m_heel r_tot rrh+rrhH rvh+rvhH rrk+rrkH rvk ri rvr]);
        end
    end

 
    for ii = i:-1:1
        [~,Iup(ii)] = max(data(j*ii-(j-1):j*ii,4));
        [~,Idn(ii)] = min(data(j*ii-(j-1):j*ii,4));
        Iup(ii) = Iup(ii) + j*(ii-1);
        Idn(ii) = Idn(ii) + j*(ii-1);
        tupwind(ii)   =  100./data(Iup(ii)  ,4);
        tdnwind(ii) = -100./data(Idn(ii),4);
        tcourse(ii)   = tupwind(ii) + tdnwind(ii);

        V1=data(Iup(ii),3)+data(Iup(ii),1).*cos(data(Iup(ii),2)*pi/180);
        V2=data(Iup(ii),1).*sin(data(Iup(ii),2)*pi/180).*cos(data(Iup(ii),5)*...
            pi/180);  
        alfa_eff(Iup(ii))=atan2(V2,V1)*180/pi;

        V1=data(Idn(ii),3)+data(Idn(ii),1).*cos(data(Idn(ii),2)*pi/180);
        V2=data(Idn(ii),1).*sin(data(Idn(ii),2)*pi/180).*cos(data(Idn(ii),5)*...
            pi/180);  
        alfa_eff(Idn(ii))=atan2(V2,V1)*180/pi;
    end


    %saving data
    fid = fopen('results.txt','w');
    fprintf(fid,' :::: VPP results :::: \n');
    fprintf(fid,'      V_tw   alfa_tw         V       VMG       phi         b',...
        '         F    Fdrive     Fside     Mheel      Rtot  Rrh+RrhH',...
        '  Rvh+RvhH  Rrk+RrkH       Rvk        Ri       Rvr\n');
    fprintf(fid,'     [m/s]     [deg]     [m/s]     [m/s]     [deg]       [m]',...
        '[-]       [N]       [N]      [Nm]       [N]       [N]       [N]',...
        '[N]       [N]       [N]       [N]\n');
    fprintf(fid,' %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f \n',data');
    for ii = 1:i
        fprintf(fid,'\n\n :::: Best VMG :::: \n');
        fprintf(fid,'      V_tw   alfa_tw         V       VMG       phi',...
            'b         F    Fdrive     Fside     Mheel      Rtot  Rrh+RrhH',...
            'Rvh+RvhH  Rrk+RrkH       Rvk        Ri       Rvr\n');
        fprintf(fid,'     [m/s]     [deg]     [m/s]     [m/s]     [deg]',...
            '[m]       [-]       [N]       [N]      [Nm]       [N]       [N]',...
            '[N]       [N]       [N]       [N]       [N]\n');
        fprintf(fid,' %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f \n',data(Iup(ii),  :)');
        fprintf(fid,' %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f \n',data(Idn(ii),:)');
        fprintf(fid,'\n\n :::: Course Time (100 m) :::: \n');
        fprintf(fid,'     upwind    %.2f [s]\n',tupwind(ii));
        fprintf(fid,'     downwind  %.2f [s]\n',tdnwind(ii));
        fprintf(fid,'     total     %.2f [s]\n',tcourse(ii));
    end
    fclose(fid);

    save data
    csvwrite('results.csv',data)

    edit results.txt;

    % plotting data
    close all
    
    ms2kt = 3600/1852;
    n2lbf = 2.2/32.2/.3048;
    figure('Name','V vs. alfa_tw')   
    for ii = i:-1:1
        plot3(repmat(param.V_tw(ii),1,j)*ms2kt,data(j*ii-(j-1):j*ii,2),...
            data(j*ii-(j-1):j*ii,3:4)*ms2kt)
        if ii == i
            hold on
        elseif ii == 1
            hold off
        end     
    end
    grid on
    set(gca,'YDir','Reverse')
    legend('V','VMG')
    xlabel('V_{tw}[kts]')
    ylabel('\alpha_{tw}[deg]')
    zlabel('kts')
    
    figure('Name','R vs. alfa_tw')  
    for ii = i:-1:1
        plot3(repmat(param.V_tw(ii),1,j)*ms2kt,data(j*ii-(j-1):j*ii,2),...
            cumsum(Rall(j*ii-(j-1):j*ii,:)*n2lbf,2))
        %100*Rperc(j*ii-(j-1):j*ii,:))
        if ii == i
            hold on
        elseif ii == 1
            hold off
        end     
    end
    grid on
    set(gca,'YDir','Reverse')
    legend('Hull - residuary','Hull - viscuos','Keel - residuary',...
        'Keel - viscous','Induced','Rudder - viscous')
    xlabel('V_{tw}[kts]')
    ylabel('\alpha_{tw}[deg]')
    %zlabel('%')
    zlabel('lbf')
    
    figure('Name','F vs. alfa_tw')
    for ii = i:-1:1
        plot3(repmat(param.V_tw(ii),1,j)*ms2kt,data(j*ii-(j-1):j*ii,2),...
            data(j*ii-(j-1):j*ii,8:9)*n2lbf)
        if ii == i
            hold on
        elseif ii == 1
            hold off
        end     
    end
    grid on
    set(gca,'YDir','Reverse')
    legend('Fdrive','Fside')
    xlabel('V_{tw}[kt]')
    ylabel('\alpha_{tw}[deg]')
    zlabel('lbf')
    
    figure('Name','Polar')
    polar(reshape(data(:,2),j,[])*pi/180,reshape(data(:,3),j,[])*ms2kt)
    hold on
    polar(reshape(data([Iup Idn],2),i,[])'*pi/180,...
       reshape(data([Iup Idn],3),i,[])'*ms2kt,'*')
    hold off
    view(90,-90)
    legend(char(num2str(param.V_tw'*ms2kt,'V_{tw} = %2.0f kts')),...
        'Location','West')
    for ii = 1:i
        [tx,ty] = pol2cart(data(Iup(ii),2)*pi/180,data(Iup(ii),3)*ms2kt);
        text(tx,ty,num2str([alfa_eff(Iup(ii)),data(Iup(ii),3)*ms2kt],...
            '(%3.0f^o,%2.1f kt)'),...
            'VerticalAlignment','Bottom','HorizontalAlign','Left','Fontsize',6)
        [tx,ty] = pol2cart(data(Idn(ii),2)*pi/180,data(Idn(ii),3)*ms2kt);
        text(tx,ty,num2str([alfa_eff(Idn(ii)),data(Idn(ii),3)*ms2kt],...
            '(%3.0f^o,%2.1f kt)'),...
            'VerticalAlignment','Top','HorizontalAlign','Left','Fontsize',6)
    end
    figure(5)
    for ii = i:-1:1
        subplot(3,1,4-ii)
        aii = data(j*ii-(j-1):j*ii,2);
        vii = data(j*ii-(j-1):j*ii,3);
        Rii = Rall(j*ii-(j-1):j*ii,:)*n2lbf;
        nn = ~isnan(vii);
        h = area(aii(nn),Rii(nn,:));
        set(h,'linestyle','none')
        if ii == 1
            xlabel('\alpha_{tw}[deg]')
        end
        if ii == 1
            legend('Hull - residuary','Hull - viscuos','Keel - residuary',...
                'Keel - viscous','Induced','Rudder - viscous')
        end
        title(['V_{tw} = ', num2str(param.V_tw(ii)*ms2kt), ' kts'])
        ylabel('lbf')
    end
    figure(6)
    sailcoefplot
end

function sailcoefplot


    LE_Cl = [
        0  0      0   0
        27  1.725  1.5 0
        50  1.5  0.5 1.5
        80  0.95 0.3 1.0
        100 0.85 0.0 0.85
        180 0    0   0];
    LE_Cdp = [
        0   0    0    0
        27  0.02 0.02 0
        50  0.15 0.25 0.25
        80  0.8  0.15 0.9
        100 1.0  0.0  1.2
        180 0.9  0.0  0.66];

    alfa_eff = linspace(0,180);

    Cl_M = pchip(LE_Cl(:,1),LE_Cl(:,2),alfa_eff);
    Cl_J = pchip(LE_Cl(:,1),LE_Cl(:,3),alfa_eff);
    Cl_S = pchip(LE_Cl(:,1),LE_Cl(:,4),alfa_eff);

    Cdp_M = pchip(LE_Cdp(:,1),LE_Cdp(:,2),alfa_eff);
    Cdp_J = pchip(LE_Cdp(:,1),LE_Cdp(:,3),alfa_eff);
    Cdp_S = pchip(LE_Cdp(:,1),LE_Cdp(:,4),alfa_eff);

    subplot(2,1,1)
    plot(alfa_eff,[Cl_M;Cl_J;Cl_S])
    legend('Main','Jib','Spin')
    ylabel('C_L')
    subplot(2,1,2)
    plot(alfa_eff,[Cdp_M;Cdp_J;Cdp_S])
    ylabel('C_D')
    xlabel('\alpha_{eff}')

end




