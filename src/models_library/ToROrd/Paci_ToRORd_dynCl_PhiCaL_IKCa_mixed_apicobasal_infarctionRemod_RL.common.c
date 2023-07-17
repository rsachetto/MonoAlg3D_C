if((is_paci > 0) && (is_paci < 4)) {
    
} else if (is_paci == 0) {
    //     Cardiac model ToR-ORd
    //     Copyright (C) 2019 Jakub Tomek. Contact: jakub.tomek.mff@gmail.com
    //
    //     This program is free software: you can redistribute it and/or modify
    //     it under the terms of the GNU General Public License as published by
    //     the Free Software Foundation, either version 3 of the License, or
    //     (at your option) any later version.
    //
    //     This program is distributed in the hope that it will be useful,
    //     but WITHOUT ANY WARRANTY; without even the implied warranty of
    //     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    //     GNU General Public License for more details.
    //
    //     You should have received a copy of the GNU General Public License
    //     along with this program.  If not, see <https://www.gnu.org/licenses/>.

    ////////////////////////////////////////////////////////////////////////////////

    //     Adapted to MonoAlg3D by Leto Riebel, Feb 2022
    //
    //     This model combines:
    //     - Updated PhiCaL for the dynamic clorid model with transmural heterogeneities
    //       as per https://github.com/jtmff/torord/blob/master/matlab/otherModelCodes/model_Torord_dynCl_updatedPhiCaL.m
    //     - Additional IKCa current as per Zhou and Wang et al., 2022, https://www.biorxiv.org/content/10.1101/2022.02.15.480392v1
    //     - Apicobasal heterogeneities in IKs
    //     - Infarct remodelling in healthy, border, infarct, and remote zone through the ischemic, and acute
    //       (Wang et al., 2021, Europace) and acute-chronic, and chronic stages (Zhou and Wang et al., 2022)
    //     - Option for Rush-Larsen - Euler method as implemented by Lucas Berg on GitHub
    //
    //     layer: 0=fast-endocardial, 1=endocardial, 2=midmyocardial, 3=epicardial
    //     infarct zone: 0=healthy, 1=infarct, 2=border
    //     infarct stage: 0=healthy, 1=acute, 2=acute-chronic, 3=chronic

    //Physical constants
    const real R = 8314.0;
    const real T = 310.0;
    const real F = 96485.0;

    //Cell geometry
    const real L = 0.01;
    const real rad = 0.0011;
    const real vcell = 1000*3.14*rad*rad*L;
    const real Ageo = 2*3.14*rad*rad+2*3.14*rad*L;
    const real Acap = 2*Ageo;
    const real vmyo = 0.68*vcell;
    const real vnsr = 0.0552*vcell;
    const real vjsr = 0.0048*vcell;
    const real vss = 0.02*vcell;

    //Concentrations
    const real nao = 140.0;
    const real cao = 1.8;
    real ko = 5.0;
    if(infarct_zone == 1) {
        ko = Ko_mult_IZ;
    }
    const real clo = 150.0;

    ////////////////////////////////////////////////////////////////////////////////

    //CaMK constants
    const real KmCaMK = 0.15;
    real aCaMK = 0.05;
    if(infarct_zone == 0) {
        aCaMK = aCaMK * aCaMK_mult_RZ;
    } else if(infarct_zone == 1) {
        aCaMK = aCaMK * aCaMK_mult_IZ;
    } else if(infarct_zone == 2) {
        aCaMK = aCaMK * aCaMK_mult_BZ;
    }
    const real bCaMK = 0.00068;
    const real CaMKo = 0.05;
    const real KmCaM = 0.0015;

    //Update CaMK
    real CaMKb = CaMKo*(1.0-CaMKt)/(1.0+KmCaM/cass);
    real CaMKa = CaMKb+CaMKt;
    rDY[36] = aCaMK*CaMKb*(CaMKb+CaMKt)-bCaMK*CaMKt;

    ////////////////////////////////////////////////////////////////////////////////

    //Reversal potentials
    real ENa = (R*T/F)*log(nao/nai);
    real EK = (R*T/F)*log(ko/ki);
    real PKNa = 0.01833;
    real EKs = (R*T/F)*log((ko+PKNa*nao)/(ki+PKNa*nai));

    //Convenient shorthand calculations
    real vffrt = v*F*F/(R*T);
    real vfrt = v*F/(R*T);
    //real frt = F/(R*T);

    ////////////////////////////////////////////////////////////////////////////////

    real fINap = (1.0/(1.0+KmCaMK/CaMKa));
    real fINaLp = (1.0/(1.0+KmCaMK/CaMKa));
    real fItop = (1.0/(1.0+KmCaMK/CaMKa));
    real fICaLp = (1.0/(1.0+KmCaMK/CaMKa));

    //INa
    //The Grandi implementation updated with INa phosphorylation.
    //m gate
    real mss = 1 / (pow((1 + exp( -(56.86 + v) / 9.03 )),2));
    real taum = 0.1292 * exp(-pow(((v+45.79)/15.54),2)) + 0.06487 * exp(-pow(((v-4.823)/51.12),2));
    rDY[9] = (mss - m) / taum;
    a_[9] = (-1.0 / taum);
    b_[9] = (mss / taum);
    //rDY[9] = mss + (m - mss)*exp(-dt/taum); //Rush-Larsen
    //h gate
    real ah = (v >= -40) * (0) + (v < -40) * (0.057 * exp( -(v + 80) / 6.8 ));
    real bh = (v >= -40) * (0.77 / (0.13*(1 + exp( -(v + 10.66) / 11.1 )))) + (v < -40) * ((2.7 * exp( 0.079 * v) + 3.1e5 * exp(0.3485 * v)));
    real tauh = 1 / (ah + bh);
    real hss = 1 / (pow((1 + exp( (v + 71.55)/7.43 )),2));
    rDY[11] = (hss - h) / tauh;
    a_[11] = (-1.0 / tauh);
    b_[11] = (hss / tauh);
    //rDY[11] = hss + (h - hss)*exp(-dt/tauh); //Rush-Larsen
    //j gate
    real aj = (v >= -40) * (0) + (v < -40) * (((-2.5428e4*exp(0.2444*v) - 6.948e-6 * exp(-0.04391*v)) * (v + 37.78)) / (1 + exp( 0.311 * (v + 79.23) )));
    real bj = (v >= -40) * ((0.6 * exp( 0.057 * v)) / (1 + exp( -0.1 * (v + 32) ))) + (v < -40) * ((0.02424 * exp( -0.01052 * v )) / (1 + exp( -0.1378 * (v + 40.14) )));
    real tauj = 1 / (aj + bj);
    real jss = 1 / (pow((1 + exp( (v + 71.55)/7.43 )),2));
    rDY[12] = (jss - j) / tauj;
    a_[12] = (-1.0 / tauj);
    b_[12] = (jss / tauj);
    //rDY[12] = jss + (j - jss)*exp(-dt/tauj); //Rush-Larsen
    //h phosphorylated
    real hssp = 1 / (pow((1 + exp( (v + 71.55 + 6)/7.43 )),2));
    rDY[10] = (hssp - hp) / tauh;
    a_[10] = (-1.0 / tauh);
    b_[10] = (hssp / tauh);
    //rDY[10] = hssp + (hp - hssp)*exp(-dt/tauh); //Rush-Larsen
    //j phosphorylated
    real taujp = 1.46 * tauj;
    rDY[13] = (jss - jp) / taujp;
    a_[13] = (-1.0 / taujp);
    b_[13] = (jss / taujp);
    //rDY[13] = jss + (jp - jss)*exp(-dt/taujp); //Rush-Larsen
    //Current
    real GNa = 11.7802 * INa_mult;
    if(infarct_zone == 0) {
        GNa = GNa * INa_mult_RZ;
    } else if(infarct_zone == 1) {
        GNa = GNa * INa_mult_IZ;
    } else if(infarct_zone == 2) {
        GNa = GNa * INa_mult_BZ;
    }
    real INa = GNa*(v-ENa)*pow(m,3)*((1.0-fINap)*h*j+fINap*hp*jp);

    //INaL
    //mL gate
    real mLss = 1.0/(1.0+exp((-(v+42.85))/5.264));
    real tm = 0.1292 * exp(-pow(((v+45.79)/15.54),2)) + 0.06487 * exp(-pow(((v-4.823)/51.12),2));
    real tmL = tm;
    rDY[14] = (mLss-mL)/tmL;
    a_[14] = (-1.0 / tmL);
    b_[14] = (mLss / tmL);
    //rDY[14] = mLss + (mL - mLss)*exp(-dt/tmL); //Rush-Larsen
    //hL gate
    real hLss = 1.0/(1.0+exp((v+87.61)/7.488));
    real thL = 200.0;
    rDY[15] = (hLss-hL)/thL;
    a_[15] = (-1.0 / thL);
    b_[15] = (hLss / thL);
    //rDY[15] = hLss + (hL - hLss)*exp(-dt/thL); //Rush-Larsen
    //hL phosphorylated
    real hLssp = 1.0/(1.0+exp((v+93.81)/7.488));
    real thLp = 3.0*thL;
    rDY[16] = (hLssp-hLp)/thLp;
    a_[16] = (-1.0 / thLp);
    b_[16] = (hLssp / thLp);
    //rDY[16] = hLssp + (hLp - hLssp)*exp(-dt/thLp); //Rush-Larsen
    //Current
    real GNaL = 0.0279 * INaL_mult;
    if(layer == 3) {
        GNaL = GNaL*0.6;
    }
    if(infarct_zone == 0) {
        GNaL = GNaL * INaL_mult_RZ;
    } else if(infarct_zone == 1) {
        GNaL = GNaL * INaL_mult_IZ;
    } else if(infarct_zone == 2) {
        GNaL = GNaL * INaL_mult_BZ;
    }
    real INaL = GNaL*(v-ENa)*mL*((1.0-fINaLp)*hL+fINaLp*hLp);

    //ITo
    //a gate
    real ass = 1.0/(1.0+exp((-(v-14.34))/14.82));
    real ta = 1.0515/(1.0/(1.2089*(1.0+exp(-(v-18.4099)/29.3814)))+3.5/(1.0+exp((v+100.0)/29.3814)));
    rDY[17] = (ass-a)/ta;
    a_[17] = (-1.0 / ta);
    b_[17] = (ass / ta);
    //rDY[17] = ass + (a - ass)*exp(-dt/ta); //Rush-Larsen
    real iss = 1.0/(1.0+exp((v+43.94)/5.711));
    real delta_epi = 1.0;
    if(layer == 3) {
        delta_epi = 1.0-(0.95/(1.0+exp((v+70.0)/5.0)));
    }
    //iF and iS gate
    real tiF = 4.562+1/(0.3933*exp((-(v+100.0))/100.0)+0.08004*exp((v+50.0)/16.59));
    real tiS = 23.62+1/(0.001416*exp((-(v+96.52))/59.05)+1.780e-8*exp((v+114.1)/8.079));
    tiF = tiF*delta_epi;
    tiS = tiS*delta_epi;
    real AiF = 1.0/(1.0+exp((v-213.6)/151.2));
    real AiS = 1.0-AiF;
    rDY[18] = (iss-iF)/tiF;
    a_[18] = (-1.0 / tiF);
    b_[18] = (iss / tiF);
    //rDY[18] = iss + (iF - iss)*exp(-dt/tiF); //Rush-Larsen
    rDY[19] = (iss-iS)/tiS;
    a_[19] = (-1.0 / tiS);
    b_[19] = (iss / tiS);
    //rDY[19] = iss + (iS - iss)*exp(-dt/tiS); //Rush-Larsen
    real i = AiF*iF+AiS*iS;
    //a phosphorylated
    real assp = 1.0/(1.0+exp((-(v-24.34))/14.82));
    rDY[20] = (assp-ap)/ta;
    a_[20] = (-1.0 / ta);
    b_[20] = (assp / ta);
    //rDY[20] = assp + (ap - assp)*exp(-dt/ta); //Rush-Larsen
    //iF and iS phosphorylated
    real dti_develop = 1.354+1.0e-4/(exp((v-167.4)/15.89)+exp(-(v-12.23)/0.2154));
    real dti_recover = 1.0-0.5/(1.0+exp((v+70.0)/20.0));
    real tiFp = dti_develop*dti_recover*tiF;
    real tiSp = dti_develop*dti_recover*tiS;
    rDY[21] = (iss-iFp)/tiFp;
    a_[21] = (-1.0 / tiFp);
    b_[21] = (iss / tiFp);
    //rDY[21] = iss + (iFp - iss)*exp(-dt/tiFp); //Rush-Larsen
    rDY[22] = (iss-iSp)/tiSp;
    a_[22] = (-1.0 / tiSp);
    b_[22] = (iss / tiSp);
    //rDY[22] = iss + (iSp - iss)*exp(-dt/tiSp); //Rush-Larsen
    real ip = AiF*iFp+AiS*iSp;
    //Current
    real Gto = 0.16 * Ito_mult;
    if(layer == 2 || layer == 3) {
        Gto = Gto*2.0;
    }
    if(infarct_zone == 0) {
        Gto = Gto * Ito_mult_RZ;
    } else if(infarct_zone == 1) {
        Gto = Gto * Ito_mult_IZ;
    } else if(infarct_zone == 2) {
        Gto = Gto * Ito_mult_BZ;
    }
    real Ito = Gto*(v-EK)*((1.0-fItop)*a*i+fItop*ap*ip);

    //ICaL
    //A variant updated by jakub, using a changed activation curve
    //It computes both ICaL in subspace and myoplasm (_i)
    //d gate
    real dss = 1.0763*exp(-1.0070*exp(-0.0829*(v)));  //magyar
    if(v > 31.4978){ //activation cannot be greater than 1
        dss = 1;
    }
    real td = 0.6+1.0/(exp(-0.05*(v+6.0))+exp(0.09*(v+14.0)));
    rDY[23] = (dss-d)/td;
    a_[23] = (-1.0 / td);
    b_[23] = (dss / td);
    //rDY[23] = dss + (d - dss)*exp(-dt/td); //Rush-Larsen
    //fs and ff gate
    real fss = 1.0/(1.0+exp((v+19.58)/3.696));
    real tff = 7.0+1.0/(0.0045*exp(-(v+20.0)/10.0)+0.0045*exp((v+20.0)/10.0));
    real tfs = 1000.0+1.0/(0.000035*exp(-(v+5.0)/4.0)+0.000035*exp((v+5.0)/6.0));
    real Aff = 0.6;
    real Afs = 1.0-Aff;
    rDY[24] = (fss-ff)/tff;
    a_[24] = (-1.0 / tff);
    b_[24] = (fss / tff);
    //rDY[24] = fss + (ff - fss)*exp(-dt/tff); //Rush-Larsen
    rDY[25] = (fss-fs)/tfs;
    a_[25] = (-1.0 / tfs);
    b_[25] = (fss / tfs);
    //rDY[25] = fss + (fs - fss)*exp(-dt/tfs); //Rush-Larsen
    real f = Aff*ff+Afs*fs;
    //fcas and fcaf gate
    real fcass = fss;
    real tfcaf = 7.0+1.0/(0.04*exp(-(v-4.0)/7.0)+0.04*exp((v-4.0)/7.0));
    real tfcas = 100.0+1.0/(0.00012*exp(-v/3.0)+0.00012*exp(v/7.0));
    real Afcaf = 0.3+0.6/(1.0+exp((v-10.0)/10.0));
    real Afcas = 1.0-Afcaf;
    rDY[26] = (fcass-fcaf)/tfcaf;
    a_[26] = (-1.0 / tfcaf);
    b_[26] = (fcass / tfcaf);
    //rDY[26] = fcass + (fcaf - fcass)*exp(-dt/tfcaf); //Rush-Larsen
    rDY[27] = (fcass-fcas)/tfcas;
    a_[27] = (-1.0 / tfcas);
    b_[27] = (fcass / tfcas);
    //rDY[27] = fcass + (fcas - fcass)*exp(-dt/tfcas); //Rush-Larsen
    real fca = Afcaf*fcaf+Afcas*fcas;
    //jca gate
    real tjca = 72.5;
    real jcass = 1.0/(1.0+exp((v+18.08)/(2.7916)));
    rDY[28] = (jcass-jca)/tjca;
    a_[28] = (-1.0 / tjca);
    b_[28] = (jcass / tjca);
    //rDY[28] = jcass + (jca - jcass)*exp(-dt/tjca); //Rush-Larsen
    //ff phosphorylated
    real tffp = 2.5*tff;
    rDY[31] = (fss-ffp)/tffp;
    a_[31] = (-1.0 / tffp);
    b_[31] = (fss / tffp);
    //rDY[31] = fss + (ffp - fss)*exp(-dt/tffp); //Rush-Larsen
    real fp = Aff*ffp+Afs*fs;
    //fcaf phosphorylated
    real tfcafp = 2.5*tfcaf;
    rDY[32] = (fcass-fcafp)/tfcafp;
    a_[32] = (-1.0 / tfcafp);
    b_[32] = (fcass / tfcafp);
    //rDY[32] = fcass + (fcafp - fcass)*exp(-dt/tfcafp); //Rush-Larsen
    real fcap = Afcaf*fcafp+Afcas*fcas;
    //SS nca
    real Kmn = 0.002;
    real k2n = 500.0;
    real km2n = jca*1;
    real anca = 1.0/(k2n/km2n+pow((1.0+Kmn/cass),4));
    rDY[29] = anca*k2n-nca*km2n;
    //myoplasmic nca
    real anca_i = 1.0/(k2n/km2n+pow((1.0+Kmn/cai),4));
    rDY[30] = anca_i*k2n-nca_i*km2n;
    //SS driving force
    real Io = 0.5*(nao + ko + clo + 4*cao)/1000; //ionic strength outside. /1000 is for things being in micromolar
    real Ii = 0.5*(nass + kss + clss + 4*cass)/1000; // ionic strength outside. /1000 is for things being in micromolar
    //The ionic strength is too high for basic DebHuc. We'll use Davies
    real dielConstant = 74.0; //water at 37C.
    real temp = 310.0; //body temp in kelvins.
    real constA = (1.82e6)*pow((dielConstant*temp),(-1.5));
    real gamma_cai = pow(10,(-constA * 4 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii)));
    real gamma_cao = pow(10,(-constA * 4 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io)));
    real gamma_nai = pow(10,(-constA * 1 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii)));
    real gamma_nao = pow(10,(-constA * 1 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io)));
    real gamma_ki = pow(10,(-constA * 1 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii)));
    real gamma_kao = pow(10,(-constA * 1 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io)));
    real PhiCaL_ss = 4.0*vffrt*(gamma_cai*cass*exp(2.0*vfrt)-gamma_cao*cao)/(exp(2.0*vfrt)-1.0);
    real PhiCaNa_ss = 1.0*vffrt*(gamma_nai*nass*exp(1.0*vfrt)-gamma_nao*nao)/(exp(1.0*vfrt)-1.0);
    real PhiCaK_ss = 1.0*vffrt*(gamma_ki*kss*exp(1.0*vfrt)-gamma_kao*ko)/(exp(1.0*vfrt)-1.0);
    //Myo driving force
    Io = 0.5*(nao + ko + clo + 4*cao)/1000 ; //ionic strength outside. /1000 is for things being in micromolar
    Ii = 0.5*(nai + ki + cli + 4*cai)/1000 ; //ionic strength outside. /1000 is for things being in micromolar
    gamma_cai = pow(10,(-constA * 4 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii)));
    gamma_cao = pow(10,(-constA * 4 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io)));
    gamma_nai = pow(10,(-constA * 1 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii)));
    gamma_nao = pow(10,(-constA * 1 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io)));
    gamma_ki = pow(10,(-constA * 1 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii)));
    gamma_kao = pow(10,(-constA * 1 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io)));
    real gammaCaoMyo = gamma_cao;
    real gammaCaiMyo = gamma_cai;
    real PhiCaL_i = 4.0*vffrt*(gamma_cai*cai*exp(2.0*vfrt)-gamma_cao*cao)/(exp(2.0*vfrt)-1.0);
    real PhiCaNa_i = 1.0*vffrt*(gamma_nai*nai*exp(1.0*vfrt)-gamma_nao*nao)/(exp(1.0*vfrt)-1.0);
    real PhiCaK_i = 1.0*vffrt*(gamma_ki*ki*exp(1.0*vfrt)-gamma_kao*ko)/(exp(1.0*vfrt)-1.0);
    //The rest
    real PCa=8.3757e-05 * 1.8969 * ICaL_mult;
    if(layer == 3) {
        PCa = PCa*1.2;
    } else if(layer == 2) {
        PCa = PCa*2;
    }
    if(infarct_zone == 0) {
        PCa = PCa * ICaL_mult_RZ;
    } else if(infarct_zone == 1) {
        PCa = PCa * ICaL_mult_IZ;
    } else if(infarct_zone == 2) {
        PCa = PCa * ICaL_mult_BZ;
    }
    real PCap = 1.1*PCa;
    real PCaNa = 0.00125*PCa*1.1737/1.8969;
    real PCaK = 3.574e-4*PCa*1.1737/1.8969;
    real PCaNap = 0.00125*PCap;
    real PCaKp = 3.574e-4*PCap;
    //SS
    real ICaL_ss = (1.0-fICaLp)*PCa*PhiCaL_ss*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCap*PhiCaL_ss*d*(fp*(1.0-nca)+jca*fcap*nca);
    real ICaNa_ss = (1.0-fICaLp)*PCaNa*PhiCaNa_ss*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCaNap*PhiCaNa_ss*d*(fp*(1.0-nca)+jca*fcap*nca);
    real ICaK_ss = (1.0-fICaLp)*PCaK*PhiCaK_ss*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCaKp*PhiCaK_ss*d*(fp*(1.0-nca)+jca*fcap*nca);
    //Myo
    real ICaL_i = (1.0-fICaLp)*PCa*PhiCaL_i*d*(f*(1.0-nca_i)+jca*fca*nca_i)+fICaLp*PCap*PhiCaL_i*d*(fp*(1.0-nca_i)+jca*fcap*nca_i);
    real ICaNa_i = (1.0-fICaLp)*PCaNa*PhiCaNa_i*d*(f*(1.0-nca_i)+jca*fca*nca_i)+fICaLp*PCaNap*PhiCaNa_i*d*(fp*(1.0-nca_i)+jca*fcap*nca_i);
    real ICaK_i = (1.0-fICaLp)*PCaK*PhiCaK_i*d*(f*(1.0-nca_i)+jca*fca*nca_i)+fICaLp*PCaKp*PhiCaK_i*d*(fp*(1.0-nca_i)+jca*fcap*nca_i);
    //And we weight ICaL (in ss) and ICaL_i
    real ICaL_fractionSS    = 0.8;
    ICaL_i = ICaL_i * (1-ICaL_fractionSS);
    ICaNa_i = ICaNa_i * (1-ICaL_fractionSS);
    ICaK_i = ICaK_i * (1-ICaL_fractionSS);
    ICaL_ss = ICaL_ss * ICaL_fractionSS;
    ICaNa_ss = ICaNa_ss * ICaL_fractionSS;
    ICaK_ss = ICaK_ss * ICaL_fractionSS;
    //Current
    real ICaL = ICaL_ss + ICaL_i;
    real ICaNa = ICaNa_ss + ICaNa_i;
    real ICaK = ICaK_ss + ICaK_i;
    real ICaL_tot = ICaL + ICaNa + ICaK;

    //IKr
    //Variant based on Lu-Vandenberg
    //transition rates
    //from c0 to c1 in l-v model,
    real alpha = 0.1161 * exp(0.2990 * vfrt);
    //from c1 to c0 in l-v/
    real beta =  0.2442 * exp(-1.604 * vfrt);
    //from c1 to c2 in l-v/
    real alpha1 = 1.25 * 0.1235;
    //from c2 to c1 in l-v/
    real beta1 =  0.1911;
    //from c2 to o/           c1 to o
    real alpha2 = 0.0578 * exp(0.9710 * vfrt);
    //from o to c2/
    real beta2 = 0.349e-3 * exp(-1.062 * vfrt);
    //from o to i
    real alphai = 0.2533 * exp(0.5953 * vfrt);
    //from i to o
    real betai = 1.25 * 0.0522 * exp(-0.8209 * vfrt);
    //from c2 to i (from c1 in orig)
    real alphac2ToI = 0.52e-4 * exp(1.525 * vfrt);
    //from i to c2
    real betaItoC2 = (beta2 * betai * alphac2ToI)/(alpha2 * alphai);
    //transitions themselves
    //for reason of backward compatibility of naming of an older version of a
    //MM IKr, c3 in code is c0 in article diagram, c2 is c1, c1 is c2.
    real c0_temp = ikr_c1 * beta - ikr_c0 * alpha; //delta for c0;
    real c1_temp = ikr_c0 * alpha + ikr_c2 * beta1 - ikr_c1 * (beta + alpha1); //c1;
    real c2_temp = ikr_c1 * alpha1 + ikr_o * beta2 + ikr_i * betaItoC2 - ikr_c2 * (beta1 + alpha2 + alphac2ToI); //subtraction is into c2, to o, to i.;
    if(c0_temp < -1)
        c0_temp = -1;
    else if(c0_temp > 1)
        c0_temp = 1;
    if(c1_temp < -1)
        c1_temp = -1;
    else if(c1_temp > 1)
        c1_temp = 1;
    if(c2_temp < -1)
        c2_temp = -1;
    else if(c2_temp > 1)
        c2_temp = 1;
    rDY[37] = c0_temp;
    rDY[38] = c1_temp;
    rDY[39] = c2_temp;
    rDY[40] = ikr_c2 * alpha2 + ikr_i * betai - ikr_o * (beta2 + alphai);
    rDY[41] = ikr_c2 * alphac2ToI + ikr_o * alphai - ikr_i * (betaItoC2 + betai);
    //Current
    real GKr = 0.0321 * sqrt(ko/5) * IKr_mult; //1st element compensates for change to ko (sqrt(5/5.4)* 0.0362)
    if(layer == 3) {
        GKr = GKr*1.3;
    } else if(layer == 2) {
        GKr=GKr*0.8;
    }
    if(infarct_zone == 0) {
        GKr = GKr * IKr_mult_RZ;
    } else if(infarct_zone == 1) {
        GKr = GKr * IKr_mult_IZ;
    } else if(infarct_zone == 2) {
        GKr = GKr * IKr_mult_BZ;
    }
    real IKr = GKr * ikr_o * (v-EK);

    //IKs
    //xs1 gate
    real xs1ss = 1.0/(1.0+exp((-(v+11.60))/8.932));
    real txs1 = 817.3+1.0/(2.326e-4*exp((v+48.28)/17.80)+0.001292*exp((-(v+210.0))/230.0));
    rDY[33] = (xs1ss-xs1)/txs1;
    a_[33] = (-1.0 / txs1);
    b_[33] = (xs1ss / txs1);
    //rDY[33] = xs1ss + (xs1 - xs1ss)*exp(-dt/txs1); //Rush-Larsen
    //xs2 gate
    real xs2ss = xs1ss;
    real txs2 = 1.0/(0.01*exp((v-50.0)/20.0)+0.0193*exp((-(v+66.54))/31.0));
    rDY[34] = (xs2ss-xs2)/txs2;
    a_[34] = (-1.0 / txs2);
    b_[34] = (xs2ss / txs2);
    //rDY[34] = xs2ss + (xs2 - xs2ss)*exp(-dt/txs2); //Rush-Larsen
    real KsCa = 1.0+0.6/(1.0+pow((3.8e-5/cai),1.4));
    //Current
    real GKs = 0.0011 * (pow(0.2, (2*apicobasal - 1))) * IKs_mult; // * apicobasal * IKs_mult //Apicobasal gradient normalised to produce scaling factor between 0.2 and 5.0
    if(layer == 3) {
        GKs = GKs*1.4;
    }
    if(infarct_zone == 0) {
        GKs = GKs * IKs_mult_RZ;
    } else if(infarct_zone == 1) {
        GKs = GKs * IKs_mult_IZ;
    } else if(infarct_zone == 2) {
        GKs = GKs * IKs_mult_BZ;
    }
    real IKs = GKs*KsCa*xs1*xs2*(v-EKs);

    //IK1
    real aK1 = 4.094/(1+exp(0.1217*(v-EK-49.934)));
    real bK1 = (15.72*exp(0.0674*(v-EK-3.257))+exp(0.0618*(v-EK-594.31)))/(1+exp(-0.1629*(v-EK+14.207)));
    real K1ss = aK1/(aK1+bK1);
    real GK1= 0.6992 * IK1_mult;
    if(layer == 3) {
        GK1 = GK1*1.2;
    } else if(layer == 2) {
        GK1=GK1*1.3;
    }
    if(infarct_zone == 0) {
        GK1 = GK1 * IK1_mult_RZ;
    } else if(infarct_zone == 1) {
        GK1 = GK1 * IK1_mult_IZ;
    } else if(infarct_zone == 2) {
        GK1 = GK1 * IK1_mult_BZ;
    }
    real IK1 = GK1*sqrt(ko/5)*K1ss*(v-EK);

    //INaCa
    //Myo
    real zca = 2.0;
    real kna1 = 15.0;
    real kna2 = 5.0;
    real kna3 = 88.12;
    real kasymm = 12.5;
    real wna = 6.0e4;
    real wca = 6.0e4;
    real wnaca = 5.0e3;
    real kcaon = 1.5e6;
    real kcaoff = 5.0e3;
    real qna = 0.5224;
    real qca = 0.1670;
    real hca = exp((qca*v*F)/(R*T));
    real hna = exp((qna*v*F)/(R*T));
    real h1 = 1+nai/kna3*(1+hna);
    real h2 = (nai*hna)/(kna3*h1);
    real h3 = 1.0/h1;
    real h4 = 1.0+nai/kna1*(1+nai/kna2);
    real h5 = nai*nai/(h4*kna1*kna2);
    real h6 = 1.0/h4;
    real h7 = 1.0+nao/kna3*(1.0+1.0/hna);
    real h8 = nao/(kna3*hna*h7);
    real h9 = 1.0/h7;
    real h10 = kasymm+1.0+nao/kna1*(1.0+nao/kna2);
    real h11 = nao*nao/(h10*kna1*kna2);
    real h12 = 1.0/h10;
    real k1 = h12*cao*kcaon;
    real k2 = kcaoff;
    real k3p = h9*wca;
    real k3pp = h8*wnaca;
    real k3 = k3p+k3pp;
    real k4p = h3*wca/hca;
    real k4pp = h2*wnaca;
    real k4 = k4p+k4pp;
    real k5 = kcaoff;
    real k6 = h6*cai*kcaon;
    real k7 = h5*h2*wna;
    real k8 = h8*h11*wna;
    real x1 = k2*k4*(k7+k6)+k5*k7*(k2+k3);
    real x2 = k1*k7*(k4+k5)+k4*k6*(k1+k8);
    real x3 = k1*k3*(k7+k6)+k8*k6*(k2+k3);
    real x4 = k2*k8*(k4+k5)+k3*k5*(k1+k8);
    real E1 = x1/(x1+x2+x3+x4);
    real E2 = x2/(x1+x2+x3+x4);
    real E3 = x3/(x1+x2+x3+x4);
    real E4 = x4/(x1+x2+x3+x4);
    real KmCaAct = 150.0e-6;
    real allo = 1.0/(1.0+pow((KmCaAct/cai),2));
    real zna = 1.0;
    real JncxNa = 3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
    real JncxCa = E2*k2-E1*k1;
    //Current
    real Gncx= 0.0034;
    if(layer == 3) {
        Gncx=Gncx*1.1;
    } else if(layer == 2) {
        Gncx=Gncx*1.4;
    }
    real INaCa_fractionSS = 0.35;
    real INaCa_i = (1-INaCa_fractionSS)*Gncx*allo*(zna*JncxNa+zca*JncxCa);
    //INaCa_ss
    h1 = 1+nass/kna3*(1+hna);
    h2 = (nass*hna)/(kna3*h1);
    h3 = 1.0/h1;
    h4 = 1.0+nass/kna1*(1+nass/kna2);
    h5 = nass*nass/(h4*kna1*kna2);
    h6 = 1.0/h4;
    h7 = 1.0+nao/kna3*(1.0+1.0/hna);
    h8 = nao/(kna3*hna*h7);
    h9 = 1.0/h7;
    h10 = kasymm+1.0+nao/kna1*(1+nao/kna2);
    h11 = nao*nao/(h10*kna1*kna2);
    h12 = 1.0/h10;
    k1 = h12*cao*kcaon;
    k2 = kcaoff;
    k3p = h9*wca;
    k3pp = h8*wnaca;
    k3 = k3p+k3pp;
    k4p = h3*wca/hca;
    k4pp = h2*wnaca;
    k4 = k4p+k4pp;
    k5 = kcaoff;
    k6 = h6*cass*kcaon;
    k7 = h5*h2*wna;
    k8 = h8*h11*wna;
    x1 = k2*k4*(k7+k6)+k5*k7*(k2+k3);
    x2 = k1*k7*(k4+k5)+k4*k6*(k1+k8);
    x3 = k1*k3*(k7+k6)+k8*k6*(k2+k3);
    x4 = k2*k8*(k4+k5)+k3*k5*(k1+k8);
    E1 = x1/(x1+x2+x3+x4);
    E2 = x2/(x1+x2+x3+x4);
    E3 = x3/(x1+x2+x3+x4);
    E4 = x4/(x1+x2+x3+x4);
    KmCaAct = 150.0e-6;
    allo = 1.0/(1.0+pow((KmCaAct/cass),2.0));
    JncxNa = 3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
    JncxCa = E2*k2-E1*k1;
    //Current
    real INaCa_ss = INaCa_fractionSS*Gncx*allo*(zna*JncxNa+zca*JncxCa);

    //INaK
    real k1p = 949.5;
    real k1m = 182.4;
    real k2p = 687.2;
    real k2m = 39.4;
    k3p = 1899.0;
    real k3m = 79300.0;
    k4p = 639.0;
    real k4m = 40.0;
    real Knai0 = 9.073;
    real Knao0 = 27.78;
    real delta = -0.1550;
    real Knai = Knai0*exp((delta*v*F)/(3.0*R*T));
    real Knao = Knao0*exp(((1.0-delta)*v*F)/(3.0*R*T));
    real Kki = 0.5;
    real Kko = 0.3582;
    real MgADP = 0.05;
    real MgATP = 9.8;
    real Kmgatp = 1.698e-7;
    real H = 1.0e-7;
    real eP = 4.2;
    real Khp = 1.698e-7;
    real Knap = 224.0;
    real Kxkur = 292.0;
    real P = eP/(1.0+H/Khp+nai/Knap+ki/Kxkur);
    real a1 = (k1p*pow((nai/Knai),3))/(pow((1.0+nai/Knai),3)+pow((1.0+ki/Kki),2)-1.0);
    real b1 = k1m*MgADP;
    real a2 = k2p;
    real b2 = (k2m*pow((nao/Knao),3))/(pow((1.0+nao/Knao),3)+pow((1.0+ko/Kko),2)-1.0);
    real a3 = (k3p*pow((ko/Kko),2))/(pow((1.0+nao/Knao),3)+pow((1.0+ko/Kko),2)-1.0);
    real b3 = (k3m*P*H)/(1.0+MgATP/Kmgatp);
    real a4 = (k4p*MgATP/Kmgatp)/(1.0+MgATP/Kmgatp);
    real b4 = (k4m*pow((ki/Kki),2))/(pow((1.0+nai/Knai),3)+pow((1.0+ki/Kki),2)-1.0);
    x1 = a4*a1*a2+b2*b4*b3+a2*b4*b3+b3*a1*a2;
    x2 = b2*b1*b4+a1*a2*a3+a3*b1*b4+a2*a3*b4;
    x3 = a2*a3*a4+b3*b2*b1+b2*b1*a4+a3*a4*b1;
    x4 = b4*b3*b2+a3*a4*a1+b2*a4*a1+b3*b2*a1;
    E1 = x1/(x1+x2+x3+x4);
    E2 = x2/(x1+x2+x3+x4);
    E3 = x3/(x1+x2+x3+x4);
    E4 = x4/(x1+x2+x3+x4);
    real zk = 1.0;
    real JnakNa = 3.0*(E1*a3-E2*b3);
    real JnakK = 2.0*(E4*b1-E3*a1);
    real Pnak = 15.4509;
    if(layer == 3) {
        Pnak = Pnak*0.9;
    } else if(layer == 2) {
        Pnak = Pnak*0.7;
    }
    real INaK = Pnak*(zna*JnakNa+zk*JnakK);

    //Minor/background currents
    //IKb
    real xkb = 1.0/(1.0+exp(-(v-10.8968)/(23.9871)));
    real GKb = 0.0189 * 0.9; //scaled by 0.9 to account for IKCa as per Zhou and Wang et al., 2022
    if(layer == 3) {
        GKb = GKb*0.6;
    }
    real IKb = GKb*xkb*(v-EK);
    //INab
    real PNab = 1.9239e-09;
    real INab = PNab*vffrt*(nai*exp(vfrt)-nao)/(exp(vfrt)-1.0);
    //ICab
    real PCab = (5.9194e-08) * 1.8969;
    if(infarct_zone == 0) {
        PCab = PCab * ICab_mult_RZ;
    } else if(infarct_zone == 1) {
        PCab = PCab * ICab_mult_IZ;
    } else if(infarct_zone == 2) {
        PCab = PCab * ICab_mult_BZ;
    }
    real ICab = PCab*4.0*vffrt*(gammaCaiMyo*cai*exp(2.0*vfrt)-gammaCaoMyo*cao)/(exp(2.0*vfrt)-1.0);
    //IpCa
    real GpCa = 5e-04;
    real IpCa = GpCa*cai/(0.0005+cai);
    //Chloride
    //I_ClCa: Ca-activated Cl Current, I_Clbk: background Cl Current
    real ecl = (R*T/F)*log(cli/clo);
    real eclss = (R*T/F)*log(clss/clo);
    real Fjunc = 1; //fraction in SS and in myoplasm - as per literature, I(Ca)Cl is in junctional subspace
    real Fsl = 1-Fjunc; //fraction in SS and in myoplasm
    real GClCa = 0.2843;
    if(infarct_zone == 0) {
        GClCa = GClCa * IClCa_mult_RZ;
    } else if(infarct_zone == 2) {
        GClCa = GClCa * IClCa_mult_BZ;
    }
    real GClB = 1.98e-3;
    real KdClCa = 0.1;
    real I_ClCa_junc = Fjunc*GClCa/(1+KdClCa/cass)*(v-eclss);
    real I_ClCa_sl = Fsl*GClCa/(1+KdClCa/cai)*(v-ecl);
    real I_ClCa = I_ClCa_junc+I_ClCa_sl;
    real I_Clbk = GClB*(v-ecl);
    //IKCa
    //As in Zhou and Wang et al., 2022
    real GKCa = 0.003;
    if(infarct_zone == 0) {
        GKCa = GKCa * IKCa_mult_RZ;
    } else if(infarct_zone == 2) {
        GKCa = GKCa * IKCa_mult_BZ;
    }
    real IKCan = 3.5;
    real KdIKCa = 6.05e-04; //(3.45e-04 for HF)
    real IKCa_fractionSS = 0.8;
    real IKCa_ss = GKCa * IKCa_fractionSS * pow(cass, IKCan)/(pow(cass, IKCan) + pow(KdIKCa, IKCan))*(v-EK);
    real IKCa_i = GKCa * (1 - IKCa_fractionSS) * pow(cai, IKCan)/(pow(cai, IKCan) + pow(KdIKCa, IKCan))*(v-EK);
    real IKCa = IKCa_ss + IKCa_i;
    // IKATP
    //real A_atp = 2.0;
    //real K_atp = 0.25;
    //real K_o_n = 5.4;
    //real fkatp = 0.0;
    //real gkatp = 4.4;
    //real akik = pow((ko/K_o_n),0.24);
    //real bkik = 1/(1 + pow((A_atp/K_atp),2));
    //real IKatp = fkatp * gkatp * akik * bkik * (v-EK);

    //Calcium handling
    //Jrel
    real fJrelp = (1.0/(1.0+KmCaMK/CaMKa));
    real jsrMidpoint = 1.7;
    real bt = 4.75;
    real a_rel = 0.5*bt;
    real Jrel_inf = a_rel*(-ICaL)/(1.0+pow((jsrMidpoint/cajsr),8));
    if(layer == 2) {
        Jrel_inf = Jrel_inf*1.7;
    }
    real tau_rel = bt/(1.0+0.0123/cajsr);
    if(tau_rel < 0.001) {
        tau_rel = 0.001;
    }
    rDY[35] = (Jrel_inf-Jrel_np)/tau_rel;
    a_[35] = (-1.0 / tau_rel);
    b_[35] = (Jrel_inf / tau_rel);
    //rDY[35] = Jrel_inf + (Jrel_np - Jrel_inf)*exp(-dt/tau_rel); //Rush-Larsen
    //Jrel phosphorylated
    real btp = 1.25*bt;
    real a_relp = 0.5*btp;
    real Jrel_infp = a_relp*(-ICaL)/(1.0+pow((jsrMidpoint/cajsr),8));
    if(layer == 2) {
        Jrel_infp = Jrel_infp*1.7;
    }
    real tau_relp = btp/(1.0+0.0123/cajsr);
    if(infarct_zone == 0) {
        tau_relp = tau_relp * tau_relp_mult_RZ;
    } else if(infarct_zone == 1) {
        tau_relp = tau_relp * tau_relp_mult_IZ;
    } else if(infarct_zone == 2) {
        tau_relp = tau_relp * tau_relp_mult_BZ;
    }
    if (tau_relp < 0.001) {
        tau_relp = 0.001;
    }
    rDY[42] = (Jrel_infp-Jrel_p)/tau_relp;
    a_[42] = (-1.0 / tau_relp);
    b_[42] = (Jrel_infp / tau_relp);
    //rDY[42] = Jrel_infp + (Jrel_p - Jrel_infp)*exp(-dt/tau_relp); //Rush-Larsen
    //Flux
    real Jrel = 1.5378 * ((1.0-fJrelp)*Jrel_np+fJrelp*Jrel_p);

    //Jup and Jleak
    real Jup_scaling = 1.0;
    if(infarct_zone == 0) {
        Jup_scaling = Jup_scaling * Iup_mult_RZ;
    } else if(infarct_zone == 2) {
        Jup_scaling = Jup_scaling * Iup_mult_BZ;
    }
    real fJupp = (1.0/(1.0+KmCaMK/CaMKa));
    real Jupnp = Jup_scaling * 0.005425*cai/(cai+0.00092);
    real Jupp = Jup_scaling * 2.75*0.005425*cai/(cai+0.00092-0.00017);
    if(layer == 3) {
        Jupnp = Jupnp*1.3;
        Jupp = Jupp*1.3;
    }
    //Fluxes
    real Jleak = Jup_scaling * 0.0048825*cansr/15.0;
    real Jup = ((1.0-fJupp)*Jupnp+fJupp*Jupp-Jleak);

    //Tranlocation flux
    real Jtr = (cansr-cajsr)/60;

    ////////////////////////////////////////////////////////////////////////////////

    //Update the membrane voltage
    // Istim
     real Istim = stim_current;

    // Itot
     real Itot = INa+INaL+Ito+ICaL+ICaNa+ICaK+IKr+IKs+IK1+INaCa_i+INaCa_ss+INaK+INab+IKb+IpCa+ICab+I_ClCa+I_Clbk+IKCa;
    rDY[0] = -(Itot + Istim);

    ////////////////////////////////////////////////////////////////////////////////

    //Calculate diffusion fluxes
    real JdiffNa = (nass-nai)/2.0;
    real JdiffK = (kss-ki)/2.0;
    real JdiffCl = (clss-cli)/2.0;
    real Jdiff = (cass-cai)/0.2;

    ////////////////////////////////////////////////////////////////////////////////

    //Calcium buffer constants
    real cmdnmax = 0.05;
    if(layer == 3) {
        cmdnmax=cmdnmax*1.3;
    }
    const real kmcmdn = 0.00238;
    const real trpnmax = 0.07;
    const real kmtrpn = 0.0005;
    const real BSRmax = 0.047;
    const real KmBSR = 0.00087;
    const real BSLmax = 1.124;
    const real KmBSL = 0.0087;
    const real csqnmax = 10.0;
    const real kmcsqn = 0.8;

    //Update intracellular concentrations, using buffers for cai, cass, cajsr
    rDY[1] = -(ICaNa_i+INa+INaL+3.0*INaCa_i+3.0*INaK+INab)*Acap/(F*vmyo)+JdiffNa*vss/vmyo;
    rDY[2] = -(ICaNa_ss+3.0*INaCa_ss)*Acap/(F*vss)-JdiffNa;
    rDY[3] = -(ICaK_i+Ito+IKr+IKs+IK1+IKb+Istim-2.0*INaK+IKCa_i)*Acap/(F*vmyo)+JdiffK*vss/vmyo;
    rDY[4] = -(ICaK_ss+IKCa_ss)*Acap/(F*vss)-JdiffK;
    real Bcai = 1.0/(1.0+cmdnmax*kmcmdn/pow((kmcmdn+cai),2)+trpnmax*kmtrpn/pow((kmtrpn+cai),2));
    rDY[5] = Bcai*(-(ICaL_i + IpCa+ICab-2.0*INaCa_i)*Acap/(2.0*F*vmyo)-Jup*vnsr/vmyo+Jdiff*vss/vmyo);
    real Bcass = 1.0/(1.0+BSRmax*KmBSR/pow((KmBSR+cass),2)+BSLmax*KmBSL/pow((KmBSL+cass),2));
    rDY[6] = Bcass*(-(ICaL_ss-2.0*INaCa_ss)*Acap/(2.0*F*vss)+Jrel*vjsr/vss-Jdiff);
    rDY[43] = - (I_Clbk + I_ClCa_sl)*Acap/(-1*F*vmyo)+JdiffCl*vss/vmyo;
    rDY[44] = - I_ClCa_junc*Acap/(-1*F*vss)-JdiffCl;
    rDY[7] = Jup-Jtr*vjsr/vnsr;
    real Bcajsr = 1.0/(1.0+csqnmax*kmcsqn/pow((kmcsqn+cajsr),2));
    rDY[8] = Bcajsr*(Jtr-Jrel);
}
