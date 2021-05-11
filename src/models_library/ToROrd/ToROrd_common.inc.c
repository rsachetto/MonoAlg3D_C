#define IFNUMBER_1(name)if((celltype==1.000000000000000e+00f)) { (name) = (cmdnmax_b*1.300000000000000e+00f);    }  else{ (name) = cmdnmax_b;    }
#define IFNUMBER_2(name)if((v_old_>=(-4.000000000000000e+01f))) { (name) = 0.000000000000000e+00f;    }  else{ (name) = (5.700000000000000e-02f*exp(((-(v_old_+8.000000000000000e+01f))/6.800000000000000e+00f)));    }
#define IFNUMBER_3(name)if((v_old_>=(-4.000000000000000e+01f))) { (name) = (7.700000000000000e-01f/(1.300000000000000e-01f*(1.000000000000000e+00f+exp(((-(v_old_+1.066000000000000e+01f))/1.110000000000000e+01f)))));    }  else{ (name) = ((2.700000000000000e+00f*exp((7.900000000000000e-02f*v_old_)))+(3.100000e+05*exp((3.485000000000000e-01f*v_old_))));    }
#define IFNUMBER_4(name)if((v_old_>=(-4.000000000000000e+01f))) { (name) = 0.000000000000000e+00f;    }  else{ (name) = (((((-2.542800e+04)*exp((2.444000000000000e-01f*v_old_)))-(6.948000e-06*exp(((-4.391000000000000e-02f)*v_old_))))*(v_old_+3.778000000000000e+01f))/(1.000000000000000e+00f+exp((3.110000000000000e-01f*(v_old_+7.923000000000000e+01f)))));    }
#define IFNUMBER_5(name)if((v_old_>=(-4.000000000000000e+01f))) { (name) = ((6.000000000000000e-01f*exp((5.700000000000000e-02f*v_old_)))/(1.000000000000000e+00f+exp(((-1.000000000000000e-01f)*(v_old_+3.200000000000000e+01f)))));    }  else{ (name) = ((2.424000000000000e-02f*exp(((-1.052000000000000e-02f)*v_old_)))/(1.000000000000000e+00f+exp(((-1.378000000000000e-01f)*(v_old_+4.014000000000000e+01f)))));    }
#define IFNUMBER_6(name)if((celltype==1.000000000000000e+00f)) { (name) = (GNaL_b*6.000000000000000e-01f);    }  else{ (name) = GNaL_b;    }
#define IFNUMBER_7(name)if((celltype==1.000000000000000e+00f)) { (name) = (1.000000000000000e+00f-(9.500000000000000e-01f/(1.000000000000000e+00f+exp(((v_old_+EKshift+7.000000000000000e+01f)/5.000000000000000e+00f)))));    }  else{ (name) = 1.000000000000000e+00f;    }
#define IFNUMBER_8(name)if((celltype==1.000000000000000e+00f)) { (name) = (Gto_b*2.000000000000000e+00f);    }  else if((celltype==2.000000000000000e+00f)){ (name) = (Gto_b*2.000000000000000e+00f);    } else{ (name) = Gto_b;    }
#define IFNUMBER_9(name)if((v_old_>=3.149780000000000e+01f)) { (name) = 1.000000000000000e+00f;    }  else{ (name) = (1.076300000000000e+00f*exp(((-1.007000000000000e+00f)*exp(((-8.290000000000000e-02f)*v_old_)))));    }
#define IFNUMBER_10(name)if((celltype==1.000000000000000e+00f)) { (name) = (PCa_b*1.200000000000000e+00f);    }  else if((celltype==2.000000000000000e+00f)){ (name) = (PCa_b*2.000000000000000e+00f);    } else{ (name) = PCa_b;    }
#define IFNUMBER_11(name)if((celltype==1.000000000000000e+00f)) { (name) = (GKr_b*1.300000000000000e+00f);    }  else if((celltype==2.000000000000000e+00f)){ (name) = (GKr_b*8.000000000000000e-01f);    } else{ (name) = GKr_b;    }
#define IFNUMBER_12(name)if((celltype==1.000000000000000e+00f)) { (name) = (GKs_b*1.400000000000000e+00f);    }  else{ (name) = GKs_b;    }
#define IFNUMBER_13(name)if((celltype==1.000000000000000e+00f)) { (name) = (GK1_b*1.200000000000000e+00f);    }  else if((celltype==2.000000000000000e+00f)){ (name) = (GK1_b*1.300000000000000e+00f);    } else{ (name) = GK1_b;    }
#define IFNUMBER_14(name)if((celltype==1.000000000000000e+00f)) { (name) = (Gncx_b*1.100000000000000e+00f);    }  else if((celltype==2.000000000000000e+00f)){ (name) = (Gncx_b*1.400000000000000e+00f);    } else{ (name) = Gncx_b;    }
#define IFNUMBER_15(name)if((celltype==1.000000000000000e+00f)) { (name) = (Pnak_b*9.000000000000000e-01f);    }  else if((celltype==2.000000000000000e+00f)){ (name) = (Pnak_b*7.000000000000000e-01f);    } else{ (name) = Pnak_b;    }
#define IFNUMBER_16(name)if((celltype==1.000000000000000e+00f)) { (name) = (GKb_b*6.000000000000000e-01f);    }  else{ (name) = GKb_b;    }
#define IFNUMBER_17(name)if((celltype==2.000000000000000e+00f)) { (name) = (calc_Jrel_inf_b*1.700000000000000e+00f);    }  else{ (name) = calc_Jrel_inf_b;    }
#define IFNUMBER_18(name)if((calc_tau_rel_b<1.000000000000000e-03f)) { (name) = 1.000000000000000e-03f;    }  else{ (name) = calc_tau_rel_b;    }
#define IFNUMBER_19(name)if((celltype==2.000000000000000e+00f)) { (name) = (calc_Jrel_infp_b*1.700000000000000e+00f);    }  else{ (name) = calc_Jrel_infp_b;    }
#define IFNUMBER_20(name)if((calc_tau_relp_b<1.000000000000000e-03f)) { (name) = 1.000000000000000e-03f;    }  else{ (name) = calc_tau_relp_b;    }
#define IFNUMBER_21(name)if((celltype==1.000000000000000e+00f)) { (name) = 1.300000000000000e+00f;    }  else{ (name) = 1.000000000000000e+00f;    }

// Parameters
const real rad = 1.100000000000000e-03f;
const real L = 1.000000000000000e-02f;
const real F = 9.648500000000000e+04f;
const real R = 8.314000000000000e+03f;
const real T = 3.100000000000000e+02f;
const real CaMKo = 5.000000000000000e-02f;
const real KmCaM = 1.500000000000000e-03f;
const real aCaMK = 5.000000000000000e-02f;
const real bCaMK = 6.800000000000000e-04f;
const real cmdnmax_b = 5.000000000000000e-02f;
const real celltype = 0.000000000000000e+00f;
const real kmcmdn = 2.380000000000000e-03f;
const real trpnmax = 7.000000000000001e-02f;
const real kmtrpn = 5.000000000000000e-04f;
const real BSRmax = 4.700000000000000e-02f;
const real KmBSR = 8.700000000000000e-04f;
const real BSLmax = 1.124000000000000e+00f;
const real KmBSL = 8.699999999999999e-03f;
const real csqnmax = 1.000000000000000e+01f;
const real kmcsqn = 8.000000000000000e-01f;
const real zna = 1.000000000000000e+00f;
const real nao = 1.400000000000000e+02f;
const real zk = 1.000000000000000e+00f;
const real ko = 5.000000000000000e+00f;
const real PKNa = 1.833000000000000e-02f;
const real zcl = -1.000000000000000e+00f;
const real clo = 1.500000000000000e+02f;
const real cli = 2.400000000000000e+01f;
const real K_o_n = 5.000000000000000e+00f;
const real A_atp = 2.000000000000000e+00f;
const real K_atp = 2.500000000000000e-01f;
const real fkatp = 0.000000000000000e+00f;
const real gkatp = 4.319500000000000e+00f;
const real KmCaMK = 1.500000000000000e-01f;
const real GNa = 1.178020000000000e+01f;
const real thL = 2.000000000000000e+02f;
const real GNaL_b = 2.790000000000000e-02f;
const real EKshift = 0.000000000000000e+00f;
const real Gto_b = 1.600000000000000e-01f;
const real offset = 0.000000000000000e+00f;
const real vShift = 0.000000000000000e+00f;
const real Aff = 6.000000000000000e-01f;
const real tjca = 7.500000000000000e+01f;
const real k2n = 5.000000000000000e+02f;
const real Kmn = 2.000000000000000e-03f;
const real cao = 1.800000000000000e+00f;
const real dielConstant = 7.400000000000000e+01f;
const real PCa_b = 8.375700000000000e-05f;
const real ICaL_fractionSS = 8.000000000000000e-01f;
const real beta_1 = 1.911000000000000e-01f;
const real alpha_1 = 1.543750000000000e-01f;
const real GKr_b = 3.210000000000000e-02f;
const real GKs_b = 1.100000000000000e-03f;
const real GK1_b = 6.992000000000000e-01f;
const real qca = 1.670000000000000e-01f;
const real qna = 5.224000000000000e-01f;
const real kna3 = 8.812000000000000e+01f;
const real kna1 = 1.500000000000000e+01f;
const real kna2 = 5.000000000000000e+00f;
const real kasymm = 1.250000000000000e+01f;
const real kcaon = 1.500000000000000e+06f;
const real kcaoff = 5.000000000000000e+03f;
const real wca = 6.000000000000000e+04f;
const real wnaca = 5.000000000000000e+03f;
const real wna = 6.000000000000000e+04f;
const real KmCaAct = 1.500000000000000e-04f;
const real Gncx_b = 3.400000000000000e-03f;
const real INaCa_fractionSS = 3.500000000000000e-01f;
const real zca = 2.000000000000000e+00f;
const real Knai0 = 9.073000000000000e+00f;
const real delta = -1.550000000000000e-01f;
const real Knao0 = 2.778000000000000e+01f;
const real eP = 4.200000000000000e+00f;
const real H = 1.000000000000000e-07f;
const real Khp = 1.698000000000000e-07f;
const real Knap = 2.240000000000000e+02f;
const real Kxkur = 2.920000000000000e+02f;
const real k1p = 9.495000000000000e+02f;
const real Kki = 5.000000000000000e-01f;
const real k1m = 1.824000000000000e+02f;
const real MgADP = 5.000000000000000e-02f;
const real k2p = 6.872000000000000e+02f;
const real k2m = 3.940000000000000e+01f;
const real Kko = 3.582000000000000e-01f;
const real k3p = 1.899000000000000e+03f;
const real k3m = 7.930000000000000e+04f;
const real MgATP = 9.800000000000001e+00f;
const real Kmgatp = 1.698000000000000e-07f;
const real k4p = 6.390000000000000e+02f;
const real k4m = 4.000000000000000e+01f;
const real Pnak_b = 1.545090000000000e+01f;
const real GKb_b = 1.890000000000000e-02f;
const real PNab = 1.923900000000000e-09f;
const real PCab = 5.919400000000000e-08f;
const real GpCa = 5.000000000000000e-04f;
const real KmCap = 5.000000000000000e-04f;
const real Fjunc = 1.000000000000000e+00f;
const real GClCa = 2.843000000000000e-01f;
const real KdClCa = 1.000000000000000e-01f;
const real GClb = 1.980000000000000e-03f;
const real tauNa = 2.000000000000000e+00f;
const real tauK = 2.000000000000000e+00f;
const real tauCa = 2.000000000000000e-01f;
const real bt = 4.750000000000000e+00f;
const real cajsr_half = 1.700000000000000e+00f;
const real Jrel_b = 1.537800000000000e+00f;
const real Jup_b = 1.000000000000000e+00f;

real calc_vcell = (1.000000000000000e+03f * 3.140000000000000e+00f * rad * rad * L); // 0
real calc_Ageo = ((2.000000000000000e+00f * 3.140000000000000e+00f * rad * rad) +
                  (2.000000000000000e+00f * 3.140000000000000e+00f * rad * L)); // 1

real calc_vffrt = ((v_old_ * F * F) / (R * T)); // 8
real calc_vfrt = ((v_old_ * F) / (R * T));      // 9
real calc_CaMKb =
    ((CaMKo * (1.000000000000000e+00f - CaMKt_old_)) / (1.000000000000000e+00f + (KmCaM / cass_old_))); // 11
real calc_cmdnmax = 0.0f;
IFNUMBER_1(calc_cmdnmax); // 14
real calc_Bcass =
    (1.000000000000000e+00f /
     (1.000000000000000e+00f + ((BSRmax * KmBSR) / pow((KmBSR + cass_old_), 2.000000000000000e+00f)) +
      ((BSLmax * KmBSL) / pow((KmBSL + cass_old_), 2.000000000000000e+00f)))); // 21
real calc_Bcajsr =
    (1.000000000000000e+00f /
     (1.000000000000000e+00f + ((csqnmax * kmcsqn) / pow((kmcsqn + cajsr_old_), 2.000000000000000e+00f)))); // 24
real calc_ENa = (((R * T) / (zna * F)) * log((nao / nai_old_)));                                            // 26
real calc_EK = (((R * T) / (zk * F)) * log((ko / ki_old_)));                                                // 27
real calc_EKs = (((R * T) / (zk * F)) * log(((ko + (PKNa * nao)) / (ki_old_ + (PKNa * nai_old_)))));        // 28
real calc_ECl = (((R * T) / (zcl * F)) * log((clo / cli)));                                                 // 29
real calc_akik = pow((ko / K_o_n), 2.400000000000000e-01f);                                                 // 30
real calc_bkik =
    (1.000000000000000e+00f / (1.000000000000000e+00f + pow((A_atp / K_atp), 2.000000000000000e+00f))); // 31
real calc_mss =
    (1.000000000000000e+00f /
     pow((1.000000000000000e+00f + exp(((-(v_old_ + 5.686000000000000e+01f)) / 9.029999999999999e+00f))),
          2.000000000000000e+00f)); // 33
real calc_tm =
    ((1.292000000000000e-01f *
      exp((-pow(((v_old_ + 4.579000000000000e+01f) / 1.554000000000000e+01f), 2.000000000000000e+00f)))) +
     (6.487000000000000e-02f *
      exp((-pow(((v_old_ - 4.823000000000000e+00f) / 5.112000000000000e+01f), 2.000000000000000e+00f))))); // 34
real calc_hss = (1.000000000000000e+00f /
                 pow((1.000000000000000e+00f + exp(((v_old_ + 7.155000000000000e+01f) / 7.430000000000000e+00f))),
                      2.000000000000000e+00f)); // 36
real calc_ah = 0.0f;
IFNUMBER_2(calc_ah); // 37
real calc_bh = 0.0f;
IFNUMBER_3(calc_bh); // 38
real calc_aj = 0.0f;
IFNUMBER_4(calc_aj); // 41
real calc_bj = 0.0f;
IFNUMBER_5(calc_bj); // 42
real calc_hssp =
    (1.000000000000000e+00f /
     pow((1.000000000000000e+00f + exp(((v_old_ + 7.755000000000000e+01f) / 7.430000000000000e+00f))),
          2.000000000000000e+00f)); // 46
real calc_mLss =
    (1.000000000000000e+00f /
     (1.000000000000000e+00f + exp(((-(v_old_ + 4.285000000000000e+01f)) / 5.264000000000000e+00f)))); // 52
real calc_tmL =
    ((1.292000000000000e-01f *
      exp((-pow(((v_old_ + 4.579000000000000e+01f) / 1.554000000000000e+01f), 2.000000000000000e+00f)))) +
     (6.487000000000000e-02f *
      exp((-pow(((v_old_ - 4.823000000000000e+00f) / 5.112000000000000e+01f), 2.000000000000000e+00f))))); // 53
real calc_hLss =
    (1.000000000000000e+00f /
     (1.000000000000000e+00f + exp(((v_old_ + 8.761000000000000e+01f) / 7.488000000000000e+00f)))); // 55
real calc_hLssp =
    (1.000000000000000e+00f /
     (1.000000000000000e+00f + exp(((v_old_ + 9.381000000000000e+01f) / 7.488000000000000e+00f)))); // 57
real calc_thLp = (3.000000000000000e+00f * thL);                                                     // 58
real calc_GNaL = 0.0f;
IFNUMBER_6(calc_GNaL); // 60
real calc_ass = (1.000000000000000e+00f /
                 (1.000000000000000e+00f +
                  exp(((-((v_old_ + EKshift) - 1.434000000000000e+01f)) / 1.482000000000000e+01f)))); // 63
real calc_ta =
    (1.051500000000000e+00f /
     ((1.000000000000000e+00f /
       (1.208900000000000e+00f * (1.000000000000000e+00f + exp(((-((v_old_ + EKshift) - 1.840990000000000e+01f)) /
                                                                 2.938140000000000e+01f))))) +
      (3.500000000000000e+00f / (1.000000000000000e+00f + exp(((v_old_ + EKshift + 1.000000000000000e+02f) /
                                                                2.938140000000000e+01f)))))); // 64
real calc_iss =
    (1.000000000000000e+00f /
     (1.000000000000000e+00f + exp(((v_old_ + EKshift + 4.394000000000000e+01f) / 5.711000000000000e+00f)))); // 66
real calc_delta_epi = 0.0f;
IFNUMBER_7(calc_delta_epi); // 67
real calc_tiF_b =
    (4.562000000000000e+00f +
     (1.000000000000000e+00f /
      ((3.933000000000000e-01f * exp(((-(v_old_ + EKshift + 1.000000000000000e+02f)) / 1.000000000000000e+02f))) +
       (8.004000000000000e-02f *
        exp(((v_old_ + EKshift + 5.000000000000000e+01f) / 1.659000000000000e+01f)))))); // 68
real calc_tiS_b =
    (2.362000000000000e+01f +
     (1.000000000000000e+00f /
      ((1.416000000000000e-03f * exp(((-(v_old_ + EKshift + 9.652000000000000e+01f)) / 5.905000000000000e+01f))) +
       (1.780000e-08 * exp(((v_old_ + EKshift + 1.141000000000000e+02f) / 8.079000000000001e+00f)))))); // 69
real calc_AiF = (1.000000000000000e+00f /
                 (1.000000000000000e+00f +
                  exp((((v_old_ + EKshift) - 2.136000000000000e+02f) / 1.512000000000000e+02f)))); // 72
real calc_assp = (1.000000000000000e+00f /
                  (1.000000000000000e+00f +
                   exp(((-((v_old_ + EKshift) - 2.434000000000000e+01f)) / 1.482000000000000e+01f)))); // 77
real calc_dti_develop =
    (1.354000000000000e+00f +
     (1.000000e-04 / (exp((((v_old_ + EKshift) - 1.674000000000000e+02f) / 1.589000000000000e+01f)) +
                      exp(((-((v_old_ + EKshift) - 1.223000000000000e+01f)) / 2.154000000000000e-01f))))); // 79
real calc_dti_recover =
    (1.000000000000000e+00f -
     (5.000000000000000e-01f / (1.000000000000000e+00f + exp(((v_old_ + EKshift + 7.000000000000000e+01f) /
                                                               2.000000000000000e+01f))))); // 80
real calc_Gto = 0.0f;
IFNUMBER_8(calc_Gto); // 86
real calc_dss = 0.0f;
IFNUMBER_9(calc_dss); // 89
real calc_td = (offset + 6.000000000000000e-01f +
                (1.000000000000000e+00f /
                 (exp(((-5.000000000000000e-02f) * (v_old_ + vShift + 6.000000000000000e+00f))) +
                  exp((9.000000000000000e-02f * (v_old_ + vShift + 1.400000000000000e+01f)))))); // 90
real calc_fss =
    (1.000000000000000e+00f /
     (1.000000000000000e+00f + exp(((v_old_ + 1.958000000000000e+01f) / 3.696000000000000e+00f)))); // 92
real calc_tff =
    (7.000000000000000e+00f +
     (1.000000000000000e+00f /
      ((4.500000000000000e-03f * exp(((-(v_old_ + 2.000000000000000e+01f)) / 1.000000000000000e+01f))) +
       (4.500000000000000e-03f * exp(((v_old_ + 2.000000000000000e+01f) / 1.000000000000000e+01f)))))); // 93
real calc_tfs =
    (1.000000000000000e+03f +
     (1.000000000000000e+00f /
      ((3.500000000000000e-05f * exp(((-(v_old_ + 5.000000000000000e+00f)) / 4.000000000000000e+00f))) +
       (3.500000000000000e-05f * exp(((v_old_ + 5.000000000000000e+00f) / 6.000000000000000e+00f)))))); // 94
real calc_Afs = (1.000000000000000e+00f - Aff);                                                          // 95
real calc_tfcaf =
    (7.000000000000000e+00f +
     (1.000000000000000e+00f /
      ((4.000000000000000e-02f * exp(((-(v_old_ - 4.000000000000000e+00f)) / 7.000000000000000e+00f))) +
       (4.000000000000000e-02f * exp(((v_old_ - 4.000000000000000e+00f) / 7.000000000000000e+00f)))))); // 100
real calc_tfcas =
    (1.000000000000000e+02f +
     (1.000000000000000e+00f / ((1.200000000000000e-04f * exp(((-v_old_) / 3.000000000000000e+00f))) +
                                (1.200000000000000e-04f * exp((v_old_ / 7.000000000000000e+00f)))))); // 101
real calc_Afcaf =
    (3.000000000000000e-01f +
     (6.000000000000000e-01f /
      (1.000000000000000e+00f + exp(((v_old_ - 1.000000000000000e+01f) / 1.000000000000000e+01f))))); // 102
real calc_jcass =
    (1.000000000000000e+00f /
     (1.000000000000000e+00f + exp(((v_old_ + 1.808000000000000e+01f) / 2.791600000000000e+00f)))); // 107
real calc_km2n = (jca_old_ * 1.000000000000000e+00f);                                                // 115
real calc_Io =
    ((5.000000000000000e-01f * (nao + ko + clo + (4.000000000000000e+00f * cao))) / 1.000000000000000e+03f); // 118
real calc_Iss = ((5.000000000000000e-01f * (nass_old_ + kss_old_ + cli + (4.000000000000000e+00f * cass_old_))) /
                 1.000000000000000e+03f);                                                // 119
real calc_constA = (1.820000e+06 * pow((dielConstant * T), (-1.500000000000000e+00f))); // 120
real calc_PCa = 0.0f;
IFNUMBER_10(calc_PCa); // 130
real calc_Ii = ((5.000000000000000e-01f * (nai_old_ + ki_old_ + cli + (4.000000000000000e+00f * cai_old_))) /
                1.000000000000000e+03f); // 142
real calc_GKr = 0.0f;
IFNUMBER_11(calc_GKr); // 168
real calc_xs1ss =
    (1.000000000000000e+00f /
     (1.000000000000000e+00f + exp(((-(v_old_ + 1.160000000000000e+01f)) / 8.932000000000000e+00f)))); // 170
real calc_txs1 =
    (8.173000000000000e+02f +
     (1.000000000000000e+00f /
      ((2.326000e-04 * exp(((v_old_ + 4.828000000000000e+01f) / 1.780000000000000e+01f))) +
       (1.292000000000000e-03f * exp(((-(v_old_ + 2.100000000000000e+02f)) / 2.300000000000000e+02f)))))); // 171
real calc_txs2 =
    (1.000000000000000e+00f /
     ((1.000000000000000e-02f * exp(((v_old_ - 5.000000000000000e+01f) / 2.000000000000000e+01f))) +
      (1.930000000000000e-02f * exp(((-(v_old_ + 6.654000000000001e+01f)) / 3.100000000000000e+01f))))); // 174
real calc_KsCa = (1.000000000000000e+00f +
                  (6.000000000000000e-01f /
                   (1.000000000000000e+00f + pow((3.800000e-05 / cai_old_), 1.400000000000000e+00f)))); // 176
real calc_GKs = 0.0f;
IFNUMBER_12(calc_GKs); // 177
real calc_GK1 = 0.0f;
IFNUMBER_13(calc_GK1); // 182
real calc_h4_i =
    (1.000000000000000e+00f + ((nai_old_ / kna1) * (1.000000000000000e+00f + (nai_old_ / kna2)))); // 189
real calc_h10_i =
    (kasymm + 1.000000000000000e+00f + ((nao / kna1) * (1.000000000000000e+00f + (nao / kna2)))); // 195
real calc_k2_i = kcaoff;                                                                          // 199
real calc_k5_i = kcaoff;                                                                          // 206
real calc_allo_i =
    (1.000000000000000e+00f / (1.000000000000000e+00f + pow((KmCaAct / cai_old_), 2.000000000000000e+00f))); // 218
real calc_Gncx = 0.0f;
IFNUMBER_14(calc_Gncx); // 221
real calc_h4_ss =
    (1.000000000000000e+00f + ((nass_old_ / kna1) * (1.000000000000000e+00f + (nass_old_ / kna2)))); // 226
real calc_h10_ss =
    (kasymm + 1.000000000000000e+00f + ((nao / kna1) * (1.000000000000000e+00f + (nao / kna2)))); // 232
real calc_k2_ss = kcaoff;                                                                         // 236
real calc_k5_ss = kcaoff;                                                                         // 243
real calc_allo_ss = (1.000000000000000e+00f /
                     (1.000000000000000e+00f + pow((KmCaAct / cass_old_), 2.000000000000000e+00f))); // 255
real calc_P = (eP / (1.000000000000000e+00f + (H / Khp) + (nai_old_ / Knap) + (ki_old_ / Kxkur)));    // 261
real calc_b1 = (k1m * MgADP);                                                                         // 263
real calc_a2 = k2p;                                                                                   // 264
real calc_a4 = (((k4p * MgATP) / Kmgatp) / (1.000000000000000e+00f + (MgATP / Kmgatp)));              // 268
real calc_Pnak = 0.0f;
IFNUMBER_15(calc_Pnak); // 280
real calc_xkb =
    (1.000000000000000e+00f /
     (1.000000000000000e+00f + exp(((-(v_old_ - 1.089680000000000e+01f)) / 2.398710000000000e+01f)))); // 282
real calc_GKb = 0.0f;
IFNUMBER_16(calc_GKb);                                                                         // 283
real calc_IpCa = ((GpCa * cai_old_) / (KmCap + cai_old_));                                     // 287
real calc_JdiffNa = ((nass_old_ - nai_old_) / tauNa);                                          // 292
real calc_JdiffK = ((kss_old_ - ki_old_) / tauK);                                              // 293
real calc_Jdiff = ((cass_old_ - cai_old_) / tauCa);                                            // 294
real calc_a_rel = ((5.000000000000000e-01f * bt) / 1.000000000000000e+00f);                    // 295
real calc_tau_rel_b = (bt / (1.000000000000000e+00f + (1.230000000000000e-02f / cajsr_old_))); // 298
real calc_btp = (1.250000000000000e+00f * bt);                                                 // 301
real calc_upScale = 0.0f;
IFNUMBER_21(calc_upScale);                                                          // 310
real calc_Jleak = ((4.882500000000000e-03f * cansr_old_) / 1.500000000000000e+01f); // 314
real calc_Jtr = ((cansr_old_ - cajsr_old_) / 6.000000000000000e+01f);               // 316
real calc_Acap = (2.000000000000000e+00f * calc_Ageo);                              // 2
real calc_vmyo = (6.800000000000000e-01f * calc_vcell);                             // 3
real calc_vnsr = (5.520000000000000e-02f * calc_vcell);                             // 4
real calc_vjsr = (4.800000000000000e-03f * calc_vcell);                             // 5
real calc_vss = (2.000000000000000e-02f * calc_vcell);                              // 6
real calc_CaMKa = (calc_CaMKb + CaMKt_old_);                                        // 12
real calc_Bcai =
    (1.000000000000000e+00f /
     (1.000000000000000e+00f + ((calc_cmdnmax * kmcmdn) / pow((kmcmdn + cai_old_), 2.000000000000000e+00f)) +
      ((trpnmax * kmtrpn) / pow((kmtrpn + cai_old_), 2.000000000000000e+00f)))); // 19
real calc_I_katp = (fkatp * gkatp * calc_akik * calc_bkik * (v_old_ - calc_EK));  // 32
real calc_th = (1.000000000000000e+00f / (calc_ah + calc_bh));                    // 39
real calc_jss = calc_hss;                                                         // 43
real calc_tj = (1.000000000000000e+00f / (calc_aj + calc_bj));                    // 44
real calc_tiF = (calc_tiF_b * calc_delta_epi);                                    // 70
real calc_tiS = (calc_tiS_b * calc_delta_epi);                                    // 71
real calc_AiS = (1.000000000000000e+00f - calc_AiF);                              // 73
real calc_f = ((Aff * ff_old_) + (calc_Afs * fs_old_));                           // 98
real calc_fcass = calc_fss;                                                       // 99
real calc_Afcas = (1.000000000000000e+00f - calc_Afcaf);                          // 103
real calc_tffp = (2.500000000000000e+00f * calc_tff);                             // 109
real calc_fp = ((Aff * ffp_old_) + (calc_Afs * fs_old_));                         // 111
real calc_tfcafp = (2.500000000000000e+00f * calc_tfcaf);                         // 112
real calc_anca_ss =
    (1.000000000000000e+00f /
     ((k2n / calc_km2n) + pow((1.000000000000000e+00f + (Kmn / cass_old_)), 4.000000000000000e+00f))); // 116
real calc_gamma_cass = exp(((-calc_constA) * 4.000000000000000e+00f *
                             ((pow(calc_Iss, 1.0 / 2.0) / (1.000000000000000e+00f + pow(calc_Iss, 1.0 / 2.0))) -
                              (3.000000000000000e-01f * calc_Iss)))); // 121
real calc_gamma_cao = exp(((-calc_constA) * 4.000000000000000e+00f *
                            ((pow(calc_Io, 1.0 / 2.0) / (1.000000000000000e+00f + pow(calc_Io, 1.0 / 2.0))) -
                             (3.000000000000000e-01f * calc_Io)))); // 122
real calc_gamma_nass = exp(((-calc_constA) * 1.000000000000000e+00f *
                             ((pow(calc_Iss, 1.0 / 2.0) / (1.000000000000000e+00f + pow(calc_Iss, 1.0 / 2.0))) -
                              (3.000000000000000e-01f * calc_Iss)))); // 123
real calc_gamma_nao = exp(((-calc_constA) * 1.000000000000000e+00f *
                            ((pow(calc_Io, 1.0 / 2.0) / (1.000000000000000e+00f + pow(calc_Io, 1.0 / 2.0))) -
                             (3.000000000000000e-01f * calc_Io)))); // 124
real calc_gamma_kss = exp(((-calc_constA) * 1.000000000000000e+00f *
                            ((pow(calc_Iss, 1.0 / 2.0) / (1.000000000000000e+00f + pow(calc_Iss, 1.0 / 2.0))) -
                             (3.000000000000000e-01f * calc_Iss)))); // 125
real calc_gamma_ko = exp(((-calc_constA) * 1.000000000000000e+00f *
                           ((pow(calc_Io, 1.0 / 2.0) / (1.000000000000000e+00f + pow(calc_Io, 1.0 / 2.0))) -
                            (3.000000000000000e-01f * calc_Io)))); // 126
real calc_PCap = (1.100000000000000e+00f * calc_PCa);              // 131
real calc_PCaNa = (1.250000000000000e-03f * calc_PCa);             // 132
real calc_PCaK = (3.574000e-04 * calc_PCa);                        // 133
real calc_anca_i = (1.000000000000000e+00f / ((k2n / calc_km2n) + pow((1.000000000000000e+00f + (Kmn / cai_old_)),
                                                                       4.000000000000000e+00f))); // 140
real calc_gamma_cai = exp(((-calc_constA) * 4.000000000000000e+00f *
                            ((pow(calc_Ii, 1.0 / 2.0) / (1.000000000000000e+00f + pow(calc_Ii, 1.0 / 2.0))) -
                             (3.000000000000000e-01f * calc_Ii)))); // 143
real calc_gamma_nai = exp(((-calc_constA) * 1.000000000000000e+00f *
                            ((pow(calc_Ii, 1.0 / 2.0) / (1.000000000000000e+00f + pow(calc_Ii, 1.0 / 2.0))) -
                             (3.000000000000000e-01f * calc_Ii)))); // 144
real calc_gamma_ki = exp(((-calc_constA) * 1.000000000000000e+00f *
                           ((pow(calc_Ii, 1.0 / 2.0) / (1.000000000000000e+00f + pow(calc_Ii, 1.0 / 2.0))) -
                            (3.000000000000000e-01f * calc_Ii))));                                         // 145
real calc_alpha = (1.161000000000000e-01f * exp((2.990000000000000e-01f * calc_vfrt)));                   // 155
real calc_beta = (2.442000000000000e-01f * exp(((-1.604000000000000e+00f) * calc_vfrt)));                 // 156
real calc_alpha_2 = (5.780000000000000e-02f * exp((9.710000000000000e-01f * calc_vfrt)));                 // 157
real calc_beta_2 = (0.349000e-03 * exp(((-1.062000000000000e+00f) * calc_vfrt)));                         // 158
real calc_alpha_i = (2.533000000000000e-01f * exp((5.953000000000001e-01f * calc_vfrt)));                 // 159
real calc_beta_i = (6.525000000000000e-02f * exp(((-8.209000000000000e-01f) * calc_vfrt)));               // 160
real calc_alpha_C2ToI = (0.520000e-04 * exp((1.525000000000000e+00f * calc_vfrt)));                       // 161
real calc_IKr = (calc_GKr * pow((ko / 5.000000000000000e+00f), 1.0 / 2.0) * O_old_ * (v_old_ - calc_EK)); // 169
real calc_xs2ss = calc_xs1ss;                                                                              // 173
real calc_IKs = (calc_GKs * calc_KsCa * xs1_old_ * xs2_old_ * (v_old_ - calc_EKs));                        // 178
real calc_aK1 = (4.094000000000000e+00f /
                 (1.000000000000000e+00f +
                  exp((1.217000000000000e-01f * ((v_old_ - calc_EK) - 4.993400000000000e+01f))))); // 179
real calc_bK1 =
    (((1.572000000000000e+01f * exp((6.740000000000000e-02f * ((v_old_ - calc_EK) - 3.257000000000000e+00f)))) +
      exp((6.180000000000000e-02f * ((v_old_ - calc_EK) - 5.943099999999999e+02f)))) /
     (1.000000000000000e+00f +
      exp(((-1.629000000000000e-01f) * ((v_old_ - calc_EK) + 1.420700000000000e+01f)))));                  // 180
real calc_hca = exp((qca * calc_vfrt));                                                                    // 184
real calc_hna = exp((qna * calc_vfrt));                                                                    // 185
real calc_h5_i = ((nai_old_ * nai_old_) / (calc_h4_i * kna1 * kna2));                                       // 190
real calc_h6_i = (1.000000000000000e+00f / calc_h4_i);                                                      // 191
real calc_h11_i = ((nao * nao) / (calc_h10_i * kna1 * kna2));                                               // 196
real calc_h12_i = (1.000000000000000e+00f / calc_h10_i);                                                    // 197
real calc_h5_ss = ((nass_old_ * nass_old_) / (calc_h4_ss * kna1 * kna2));                                   // 227
real calc_h6_ss = (1.000000000000000e+00f / calc_h4_ss);                                                    // 228
real calc_h11_ss = ((nao * nao) / (calc_h10_ss * kna1 * kna2));                                             // 233
real calc_h12_ss = (1.000000000000000e+00f / calc_h10_ss);                                                  // 234
real calc_Knai = (Knai0 * exp(((delta * calc_vfrt) / 3.000000000000000e+00f)));                            // 259
real calc_Knao = (Knao0 * exp((((1.000000000000000e+00f - delta) * calc_vfrt) / 3.000000000000000e+00f))); // 260
real calc_b3 = ((k3m * calc_P * H) / (1.000000000000000e+00f + (MgATP / Kmgatp)));                          // 267
real calc_IKb = (calc_GKb * calc_xkb * (v_old_ - calc_EK));                                                 // 284
real calc_INab = ((PNab * calc_vffrt * ((nai_old_ * exp(calc_vfrt)) - nao)) /
                  (exp(calc_vfrt) - 1.000000000000000e+00f)); // 285
real calc_IClCa_junc =
    (((Fjunc * GClCa) / (1.000000000000000e+00f + (KdClCa / cass_old_))) * (v_old_ - calc_ECl)); // 288
real calc_IClCa_sl =
    ((((1.000000000000000e+00f - Fjunc) * GClCa) / (1.000000000000000e+00f + (KdClCa / cai_old_))) *
     (v_old_ - calc_ECl));                     // 289
real calc_IClb = (GClb * (v_old_ - calc_ECl)); // 291
real calc_tau_rel = 0.0f;
IFNUMBER_18(calc_tau_rel);                                                                                    // 299
real calc_a_relp = ((5.000000000000000e-01f * calc_btp) / 1.000000000000000e+00f);                            // 302
real calc_tau_relp_b = (calc_btp / (1.000000000000000e+00f + (1.230000000000000e-02f / cajsr_old_)));         // 305
real calc_Jupnp = ((calc_upScale * 5.425000000000000e-03f * cai_old_) / (cai_old_ + 9.200000000000000e-04f)); // 311
real calc_Jupp = ((calc_upScale * 2.750000000000000e+00f * 5.425000000000000e-03f * cai_old_) /
                  ((cai_old_ + 9.200000000000000e-04f) - 1.700000000000000e-04f));              // 312
real calc_tjp = (1.460000000000000e+00f * calc_tj);                                             // 48
real calc_fINap = (1.000000000000000e+00f / (1.000000000000000e+00f + (KmCaMK / calc_CaMKa)));  // 50
real calc_fINaLp = (1.000000000000000e+00f / (1.000000000000000e+00f + (KmCaMK / calc_CaMKa))); // 61
real calc_i = ((calc_AiF * iF_old_) + (calc_AiS * iS_old_));                                    // 76
real calc_tiFp = (calc_dti_develop * calc_dti_recover * calc_tiF);                              // 81
real calc_tiSp = (calc_dti_develop * calc_dti_recover * calc_tiS);                              // 82
real calc_ip = ((calc_AiF * iFp_old_) + (calc_AiS * iSp_old_));                                 // 85
real calc_fItop = (1.000000000000000e+00f / (1.000000000000000e+00f + (KmCaMK / calc_CaMKa)));  // 87
real calc_fca = ((calc_Afcaf * fcaf_old_) + (calc_Afcas * fcas_old_));                          // 106
real calc_fcap = ((calc_Afcaf * fcafp_old_) + (calc_Afcas * fcas_old_));                        // 114
real calc_PhiCaL_ss =
    ((4.000000000000000e+00f * calc_vffrt *
      ((calc_gamma_cass * cass_old_ * exp((2.000000000000000e+00f * calc_vfrt))) - (calc_gamma_cao * cao))) /
     (exp((2.000000000000000e+00f * calc_vfrt)) - 1.000000000000000e+00f)); // 127
real calc_PhiCaNa_ss =
    ((1.000000000000000e+00f * calc_vffrt *
      ((calc_gamma_nass * nass_old_ * exp((1.000000000000000e+00f * calc_vfrt))) - (calc_gamma_nao * nao))) /
     (exp((1.000000000000000e+00f * calc_vfrt)) - 1.000000000000000e+00f)); // 128
real calc_PhiCaK_ss =
    ((1.000000000000000e+00f * calc_vffrt *
      ((calc_gamma_kss * kss_old_ * exp((1.000000000000000e+00f * calc_vfrt))) - (calc_gamma_ko * ko))) /
     (exp((1.000000000000000e+00f * calc_vfrt)) - 1.000000000000000e+00f));                    // 129
real calc_PCaNap = (1.250000000000000e-03f * calc_PCap);                                        // 134
real calc_PCaKp = (3.574000e-04 * calc_PCap);                                                   // 135
real calc_fICaLp = (1.000000000000000e+00f / (1.000000000000000e+00f + (KmCaMK / calc_CaMKa))); // 136
real calc_PhiCaL_i =
    ((4.000000000000000e+00f * calc_vffrt *
      ((calc_gamma_cai * cai_old_ * exp((2.000000000000000e+00f * calc_vfrt))) - (calc_gamma_cao * cao))) /
     (exp((2.000000000000000e+00f * calc_vfrt)) - 1.000000000000000e+00f)); // 146
real calc_PhiCaNa_i =
    ((1.000000000000000e+00f * calc_vffrt *
      ((calc_gamma_nai * nai_old_ * exp((1.000000000000000e+00f * calc_vfrt))) - (calc_gamma_nao * nao))) /
     (exp((1.000000000000000e+00f * calc_vfrt)) - 1.000000000000000e+00f)); // 147
real calc_PhiCaK_i =
    ((1.000000000000000e+00f * calc_vffrt *
      ((calc_gamma_ki * ki_old_ * exp((1.000000000000000e+00f * calc_vfrt))) - (calc_gamma_ko * ko))) /
     (exp((1.000000000000000e+00f * calc_vfrt)) - 1.000000000000000e+00f));                             // 148
real calc_beta_ItoC2 = ((calc_beta_2 * calc_beta_i * calc_alpha_C2ToI) / (calc_alpha_2 * calc_alpha_i)); // 162
real calc_K1ss = (calc_aK1 / (calc_aK1 + calc_bK1));                                                     // 181
real calc_h1_i = (1.000000000000000e+00f + ((nai_old_ / kna3) * (1.000000000000000e+00f + calc_hna)));   // 186
real calc_h7_i = (1.000000000000000e+00f +
                  ((nao / kna3) * (1.000000000000000e+00f + (1.000000000000000e+00f / calc_hna))));      // 192
real calc_k1_i = (calc_h12_i * cao * kcaon);                                                             // 198
real calc_k6_i = (calc_h6_i * cai_old_ * kcaon);                                                         // 207
real calc_h1_ss = (1.000000000000000e+00f + ((nass_old_ / kna3) * (1.000000000000000e+00f + calc_hna))); // 223
real calc_h7_ss = (1.000000000000000e+00f +
                   ((nao / kna3) * (1.000000000000000e+00f + (1.000000000000000e+00f / calc_hna)))); // 229
real calc_k1_ss = (calc_h12_ss * cao * kcaon);                                                       // 235
real calc_k6_ss = (calc_h6_ss * cass_old_ * kcaon);                                                  // 244
real calc_a1 = ((k1p * pow((nai_old_ / calc_Knai), 3.000000000000000e+00f)) /
                ((pow((1.000000000000000e+00f + (nai_old_ / calc_Knai)), 3.000000000000000e+00f) +
                  pow((1.000000000000000e+00f + (ki_old_ / Kki)), 2.000000000000000e+00f)) -
                 1.000000000000000e+00f)); // 262
real calc_b2 = ((k2m * pow((nao / calc_Knao), 3.000000000000000e+00f)) /
                ((pow((1.000000000000000e+00f + (nao / calc_Knao)), 3.000000000000000e+00f) +
                  pow((1.000000000000000e+00f + (ko / Kko)), 2.000000000000000e+00f)) -
                 1.000000000000000e+00f)); // 265
real calc_a3 = ((k3p * pow((ko / Kko), 2.000000000000000e+00f)) /
                ((pow((1.000000000000000e+00f + (nao / calc_Knao)), 3.000000000000000e+00f) +
                  pow((1.000000000000000e+00f + (ko / Kko)), 2.000000000000000e+00f)) -
                 1.000000000000000e+00f)); // 266
real calc_b4 = ((k4m * pow((ki_old_ / Kki), 2.000000000000000e+00f)) /
                ((pow((1.000000000000000e+00f + (nai_old_ / calc_Knai)), 3.000000000000000e+00f) +
                  pow((1.000000000000000e+00f + (ki_old_ / Kki)), 2.000000000000000e+00f)) -
                 1.000000000000000e+00f)); // 269
real calc_ICab =
    ((PCab * 4.000000000000000e+00f * calc_vffrt *
      ((calc_gamma_cai * cai_old_ * exp((2.000000000000000e+00f * calc_vfrt))) - (calc_gamma_cao * cao))) /
     (exp((2.000000000000000e+00f * calc_vfrt)) - 1.000000000000000e+00f)); // 286
real calc_IClCa = (calc_IClCa_junc + calc_IClCa_sl);                         // 290
real calc_tau_relp = 0.0f;
IFNUMBER_20(calc_tau_relp);                                                                     // 306
real calc_fJrelp = (1.000000000000000e+00f / (1.000000000000000e+00f + (KmCaMK / calc_CaMKa))); // 308
real calc_fJupp = (1.000000000000000e+00f / (1.000000000000000e+00f + (KmCaMK / calc_CaMKa)));  // 313
real calc_INa =
    (GNa * (v_old_ - calc_ENa) * pow(m_old_, 3.000000000000000e+00f) *
     (((1.000000000000000e+00f - calc_fINap) * h_old_ * j_old_) + (calc_fINap * hp_old_ * jp_old_))); // 51
real calc_INaL = (calc_GNaL * (v_old_ - calc_ENa) * mL_old_ *
                  (((1.000000000000000e+00f - calc_fINaLp) * hL_old_) + (calc_fINaLp * hLp_old_))); // 62
real calc_Ito =
    (calc_Gto * (v_old_ - calc_EK) *
     (((1.000000000000000e+00f - calc_fItop) * a_old_ * calc_i) + (calc_fItop * ap_old_ * calc_ip))); // 88
real calc_ICaL_ss =
    (ICaL_fractionSS *
     (((1.000000000000000e+00f - calc_fICaLp) * calc_PCa * calc_PhiCaL_ss * d_old_ *
       ((calc_f * (1.000000000000000e+00f - nca_ss_old_)) + (jca_old_ * calc_fca * nca_ss_old_))) +
      (calc_fICaLp * calc_PCap * calc_PhiCaL_ss * d_old_ *
       ((calc_fp * (1.000000000000000e+00f - nca_ss_old_)) + (jca_old_ * calc_fcap * nca_ss_old_))))); // 137
real calc_ICaNa_ss =
    (ICaL_fractionSS *
     (((1.000000000000000e+00f - calc_fICaLp) * calc_PCaNa * calc_PhiCaNa_ss * d_old_ *
       ((calc_f * (1.000000000000000e+00f - nca_ss_old_)) + (jca_old_ * calc_fca * nca_ss_old_))) +
      (calc_fICaLp * calc_PCaNap * calc_PhiCaNa_ss * d_old_ *
       ((calc_fp * (1.000000000000000e+00f - nca_ss_old_)) + (jca_old_ * calc_fcap * nca_ss_old_))))); // 138
real calc_ICaK_ss =
    (ICaL_fractionSS *
     (((1.000000000000000e+00f - calc_fICaLp) * calc_PCaK * calc_PhiCaK_ss * d_old_ *
       ((calc_f * (1.000000000000000e+00f - nca_ss_old_)) + (jca_old_ * calc_fca * nca_ss_old_))) +
      (calc_fICaLp * calc_PCaKp * calc_PhiCaK_ss * d_old_ *
       ((calc_fp * (1.000000000000000e+00f - nca_ss_old_)) + (jca_old_ * calc_fcap * nca_ss_old_))))); // 139
real calc_ICaL_i =
    ((1.000000000000000e+00f - ICaL_fractionSS) *
     (((1.000000000000000e+00f - calc_fICaLp) * calc_PCa * calc_PhiCaL_i * d_old_ *
       ((calc_f * (1.000000000000000e+00f - nca_i_old_)) + (jca_old_ * calc_fca * nca_i_old_))) +
      (calc_fICaLp * calc_PCap * calc_PhiCaL_i * d_old_ *
       ((calc_fp * (1.000000000000000e+00f - nca_i_old_)) + (jca_old_ * calc_fcap * nca_i_old_))))); // 149
real calc_ICaNa_i =
    ((1.000000000000000e+00f - ICaL_fractionSS) *
     (((1.000000000000000e+00f - calc_fICaLp) * calc_PCaNa * calc_PhiCaNa_i * d_old_ *
       ((calc_f * (1.000000000000000e+00f - nca_i_old_)) + (jca_old_ * calc_fca * nca_i_old_))) +
      (calc_fICaLp * calc_PCaNap * calc_PhiCaNa_i * d_old_ *
       ((calc_fp * (1.000000000000000e+00f - nca_i_old_)) + (jca_old_ * calc_fcap * nca_i_old_))))); // 150
real calc_ICaK_i =
    ((1.000000000000000e+00f - ICaL_fractionSS) *
     (((1.000000000000000e+00f - calc_fICaLp) * calc_PCaK * calc_PhiCaK_i * d_old_ *
       ((calc_f * (1.000000000000000e+00f - nca_i_old_)) + (jca_old_ * calc_fca * nca_i_old_))) +
      (calc_fICaLp * calc_PCaKp * calc_PhiCaK_i * d_old_ *
       ((calc_fp * (1.000000000000000e+00f - nca_i_old_)) + (jca_old_ * calc_fcap * nca_i_old_)))));          // 151
real calc_IK1 = (calc_GK1 * pow((ko / 5.000000000000000e+00f), 1.0 / 2.0) * calc_K1ss * (v_old_ - calc_EK)); // 183
real calc_h2_i = ((nai_old_ * calc_hna) / (kna3 * calc_h1_i));                                                // 187
real calc_h3_i = (1.000000000000000e+00f / calc_h1_i);                                                        // 188
real calc_h8_i = (nao / (kna3 * calc_hna * calc_h7_i));                                                       // 193
real calc_h9_i = (1.000000000000000e+00f / calc_h7_i);                                                        // 194
real calc_h2_ss = ((nass_old_ * calc_hna) / (kna3 * calc_h1_ss));                                             // 224
real calc_h3_ss = (1.000000000000000e+00f / calc_h1_ss);                                                      // 225
real calc_h8_ss = (nao / (kna3 * calc_hna * calc_h7_ss));                                                     // 230
real calc_h9_ss = (1.000000000000000e+00f / calc_h7_ss);                                                      // 231
real calc_x1 = ((calc_a4 * calc_a1 * calc_a2) + (calc_b2 * calc_b4 * calc_b3) + (calc_a2 * calc_b4 * calc_b3) +
                (calc_b3 * calc_a1 * calc_a2)); // 270
real calc_x2 = ((calc_b2 * calc_b1 * calc_b4) + (calc_a1 * calc_a2 * calc_a3) + (calc_a3 * calc_b1 * calc_b4) +
                (calc_a2 * calc_a3 * calc_b4)); // 271
real calc_x3 = ((calc_a2 * calc_a3 * calc_a4) + (calc_b3 * calc_b2 * calc_b1) + (calc_b2 * calc_b1 * calc_a4) +
                (calc_a3 * calc_a4 * calc_b1)); // 272
real calc_x4 = ((calc_b4 * calc_b3 * calc_b2) + (calc_a3 * calc_a4 * calc_a1) + (calc_b2 * calc_a4 * calc_a1) +
                (calc_b3 * calc_b2 * calc_a1)); // 273
real calc_Jrel =
    (Jrel_b * (((1.000000000000000e+00f - calc_fJrelp) * Jrel_np_old_) + (calc_fJrelp * Jrel_p_old_))); // 309
real calc_Jup = (Jup_b * ((((1.000000000000000e+00f - calc_fJupp) * calc_Jupnp) + (calc_fJupp * calc_Jupp)) -
                          calc_Jleak));                             // 315
real calc_ICaL = (calc_ICaL_ss + calc_ICaL_i);                      // 152
real calc_ICaNa = (calc_ICaNa_ss + calc_ICaNa_i);                   // 153
real calc_ICaK = (calc_ICaK_ss + calc_ICaK_i);                      // 154
real calc_k3p_i = (calc_h9_i * wca);                                // 200
real calc_k3pp_i = (calc_h8_i * wnaca);                             // 201
real calc_k4p_i = ((calc_h3_i * wca) / calc_hca);                   // 203
real calc_k4pp_i = (calc_h2_i * wnaca);                             // 204
real calc_k7_i = (calc_h5_i * calc_h2_i * wna);                     // 208
real calc_k8_i = (calc_h8_i * calc_h11_i * wna);                    // 209
real calc_k3p_ss = (calc_h9_ss * wca);                              // 237
real calc_k3pp_ss = (calc_h8_ss * wnaca);                           // 238
real calc_k4p_ss = ((calc_h3_ss * wca) / calc_hca);                 // 240
real calc_k4pp_ss = (calc_h2_ss * wnaca);                           // 241
real calc_k7_ss = (calc_h5_ss * calc_h2_ss * wna);                  // 245
real calc_k8_ss = (calc_h8_ss * calc_h11_ss * wna);                 // 246
real calc_E1 = (calc_x1 / (calc_x1 + calc_x2 + calc_x3 + calc_x4)); // 274
real calc_E2 = (calc_x2 / (calc_x1 + calc_x2 + calc_x3 + calc_x4)); // 275
real calc_E3 = (calc_x3 / (calc_x1 + calc_x2 + calc_x3 + calc_x4)); // 276
real calc_E4 = (calc_x4 / (calc_x1 + calc_x2 + calc_x3 + calc_x4)); // 277
real calc_Jrel_inf_b = ((((-calc_a_rel) * calc_ICaL_ss) / 1.000000000000000e+00f) /
                        (1.000000000000000e+00f + pow((cajsr_half / cajsr_old_), 8.000000000000000e+00f))); // 296
real calc_Jrel_infp_b = ((((-calc_a_relp) * calc_ICaL_ss) / 1.000000000000000e+00f) /
                         (1.000000000000000e+00f + pow((cajsr_half / cajsr_old_), 8.000000000000000e+00f))); // 303
real calc_k3_i = (calc_k3p_i + calc_k3pp_i);                                                                  // 202
real calc_k4_i = (calc_k4p_i + calc_k4pp_i);                                                                  // 205
real calc_k3_ss = (calc_k3p_ss + calc_k3pp_ss);                                                               // 239
real calc_k4_ss = (calc_k4p_ss + calc_k4pp_ss);                                                               // 242
real calc_JnakNa = (3.000000000000000e+00f * ((calc_E1 * calc_a3) - (calc_E2 * calc_b3)));                    // 278
real calc_JnakK = (2.000000000000000e+00f * ((calc_E4 * calc_b1) - (calc_E3 * calc_a1)));                     // 279
real calc_Jrel_inf = 0.0f;
IFNUMBER_17(calc_Jrel_inf); // 297
real calc_Jrel_infp = 0.0f;
IFNUMBER_19(calc_Jrel_infp); // 304
real calc_x1_i =
    ((calc_k2_i * calc_k4_i * (calc_k7_i + calc_k6_i)) + (calc_k5_i * calc_k7_i * (calc_k2_i + calc_k3_i))); // 210
real calc_x2_i =
    ((calc_k1_i * calc_k7_i * (calc_k4_i + calc_k5_i)) + (calc_k4_i * calc_k6_i * (calc_k1_i + calc_k8_i))); // 211
real calc_x3_i =
    ((calc_k1_i * calc_k3_i * (calc_k7_i + calc_k6_i)) + (calc_k8_i * calc_k6_i * (calc_k2_i + calc_k3_i))); // 212
real calc_x4_i =
    ((calc_k2_i * calc_k8_i * (calc_k4_i + calc_k5_i)) + (calc_k3_i * calc_k5_i * (calc_k1_i + calc_k8_i))); // 213
real calc_x1_ss = ((calc_k2_ss * calc_k4_ss * (calc_k7_ss + calc_k6_ss)) +
                   (calc_k5_ss * calc_k7_ss * (calc_k2_ss + calc_k3_ss))); // 247
real calc_x2_ss = ((calc_k1_ss * calc_k7_ss * (calc_k4_ss + calc_k5_ss)) +
                   (calc_k4_ss * calc_k6_ss * (calc_k1_ss + calc_k8_ss))); // 248
real calc_x3_ss = ((calc_k1_ss * calc_k3_ss * (calc_k7_ss + calc_k6_ss)) +
                   (calc_k8_ss * calc_k6_ss * (calc_k2_ss + calc_k3_ss))); // 249
real calc_x4_ss = ((calc_k2_ss * calc_k8_ss * (calc_k4_ss + calc_k5_ss)) +
                   (calc_k3_ss * calc_k5_ss * (calc_k1_ss + calc_k8_ss)));            // 250
real calc_INaK = (calc_Pnak * ((zna * calc_JnakNa) + (zk * calc_JnakK)));             // 281
real calc_E1_i = (calc_x1_i / (calc_x1_i + calc_x2_i + calc_x3_i + calc_x4_i));       // 214
real calc_E2_i = (calc_x2_i / (calc_x1_i + calc_x2_i + calc_x3_i + calc_x4_i));       // 215
real calc_E3_i = (calc_x3_i / (calc_x1_i + calc_x2_i + calc_x3_i + calc_x4_i));       // 216
real calc_E4_i = (calc_x4_i / (calc_x1_i + calc_x2_i + calc_x3_i + calc_x4_i));       // 217
real calc_E1_ss = (calc_x1_ss / (calc_x1_ss + calc_x2_ss + calc_x3_ss + calc_x4_ss)); // 251
real calc_E2_ss = (calc_x2_ss / (calc_x1_ss + calc_x2_ss + calc_x3_ss + calc_x4_ss)); // 252
real calc_E3_ss = (calc_x3_ss / (calc_x1_ss + calc_x2_ss + calc_x3_ss + calc_x4_ss)); // 253
real calc_E4_ss = (calc_x4_ss / (calc_x1_ss + calc_x2_ss + calc_x3_ss + calc_x4_ss)); // 254
real calc_JncxNa_i =
    (((3.000000000000000e+00f * ((calc_E4_i * calc_k7_i) - (calc_E1_i * calc_k8_i))) + (calc_E3_i * calc_k4pp_i)) -
     (calc_E2_i * calc_k3pp_i));                                          // 219
real calc_JncxCa_i = ((calc_E2_i * calc_k2_i) - (calc_E1_i * calc_k1_i)); // 220
real calc_JncxNa_ss = (((3.000000000000000e+00f * ((calc_E4_ss * calc_k7_ss) - (calc_E1_ss * calc_k8_ss))) +
                        (calc_E3_ss * calc_k4pp_ss)) -
                       (calc_E2_ss * calc_k3pp_ss));                           // 256
real calc_JncxCa_ss = ((calc_E2_ss * calc_k2_ss) - (calc_E1_ss * calc_k1_ss)); // 257
real calc_INaCa_i = ((1.000000000000000e+00f - INaCa_fractionSS) * calc_Gncx * calc_allo_i *
                     ((zna * calc_JncxNa_i) + (zca * calc_JncxCa_i))); // 222
real calc_INaCa_ss =
    (INaCa_fractionSS * calc_Gncx * calc_allo_ss * ((zna * calc_JncxNa_ss) + (zca * calc_JncxCa_ss))); // 258

rDY_[0] = (-(calc_INa + calc_INaL + calc_Ito + calc_ICaL + calc_ICaNa + calc_ICaK + calc_IKr + calc_IKs + calc_IK1 +
calc_INaCa_i + calc_INaCa_ss + calc_INaK + calc_INab + calc_IKb + calc_IpCa + calc_ICab + calc_IClCa +
calc_IClb + calc_I_katp + stim_current));
rDY_[1] = ((aCaMK * calc_CaMKb * (calc_CaMKb + CaMKt_old_)) - (bCaMK * CaMKt_old_));
rDY_[2] = ((((-(calc_INa + calc_INaL + (3.000000000000000e+00f * calc_INaCa_i) + calc_ICaNa_i +
(3.000000000000000e+00f * calc_INaK) + calc_INab)) *
calc_Acap) /
(F * calc_vmyo)) +
((calc_JdiffNa * calc_vss) / calc_vmyo));
rDY_[3] =
((((-(calc_ICaNa_ss + (3.000000000000000e+00f * calc_INaCa_ss))) * calc_Acap) / (F * calc_vss)) - calc_JdiffNa);

//rDY_[4] = ((((-(((calc_Ito + calc_IKr + calc_IKs + calc_IK1 + calc_IKb + calc_I_katp + stim_current) - (2.000000000000000e+00f * calc_INaK)) + calc_ICaK_i))*calc_Acap)/(F * calc_vmyo)) +((calc_JdiffK * calc_vss) / calc_vmyo));
rDY_[4] =-(calc_ICaK_i+calc_Ito+calc_IKr+calc_IKs+calc_IK1+calc_IKb+stim_current-2.0*calc_INaK)*calc_Acap/(F*calc_vmyo)+calc_JdiffK*calc_vss/calc_vmyo;

rDY_[5] = ((((-calc_ICaK_ss) * calc_Acap) / (F * calc_vss)) - calc_JdiffK);
rDY_[6] = (calc_Bcai *
(((((-((calc_ICaL_i + calc_IpCa + calc_ICab) - (2.000000000000000e+00f * calc_INaCa_i))) * calc_Acap) /
(2.000000000000000e+00f * F * calc_vmyo)) -
((calc_Jup * calc_vnsr) / calc_vmyo)) +
((calc_Jdiff * calc_vss) / calc_vmyo)));
rDY_[7] = (calc_Bcass * (((((-(calc_ICaL_ss - (2.000000000000000e+00f * calc_INaCa_ss))) * calc_Acap) /
(2.000000000000000e+00f * F * calc_vss)) +
((calc_Jrel * calc_vjsr) / calc_vss)) -
calc_Jdiff));
rDY_[8] = (calc_Jup - ((calc_Jtr * calc_vjsr) / calc_vnsr));
rDY_[9] = (calc_Bcajsr * (calc_Jtr - calc_Jrel));
rDY_[10] = ((calc_mss - m_old_) / calc_tm);
rDY_[11] = ((calc_hss - h_old_) / calc_th);
rDY_[12] = ((calc_jss - j_old_) / calc_tj);
rDY_[13] = ((calc_hssp - hp_old_) / calc_th);
rDY_[14] = ((calc_jss - jp_old_) / calc_tjp);
rDY_[15] = ((calc_mLss - mL_old_) / calc_tmL);
rDY_[16] = ((calc_hLss - hL_old_) / thL);
rDY_[17] = ((calc_hLssp - hLp_old_) / calc_thLp);
rDY_[18] = ((calc_ass - a_old_) / calc_ta);
rDY_[19] = ((calc_iss - iF_old_) / calc_tiF);
rDY_[20] = ((calc_iss - iS_old_) / calc_tiS);
rDY_[21] = ((calc_assp - ap_old_) / calc_ta);
rDY_[22] = ((calc_iss - iFp_old_) / calc_tiFp);
rDY_[23] = ((calc_iss - iSp_old_) / calc_tiSp);
rDY_[24] = ((calc_dss - d_old_) / calc_td);
rDY_[25] = ((calc_fss - ff_old_) / calc_tff);
rDY_[26] = ((calc_fss - fs_old_) / calc_tfs);
rDY_[27] = ((calc_fcass - fcaf_old_) / calc_tfcaf);
rDY_[28] = ((calc_fcass - fcas_old_) / calc_tfcas);
rDY_[29] = ((calc_jcass - jca_old_) / tjca);
rDY_[30] = ((calc_fss - ffp_old_) / calc_tffp);
rDY_[31] = ((calc_fcass - fcafp_old_) / calc_tfcafp);
rDY_[32] = ((calc_anca_ss * k2n) - (nca_ss_old_ * calc_km2n));
rDY_[33] = ((calc_anca_i * k2n) - (nca_i_old_ * calc_km2n));
rDY_[34] = ((calc_beta * C2_old_) - (calc_alpha * C3_old_));
rDY_[35] = (((calc_alpha * C3_old_) + (beta_1 * C1_old_)) - ((calc_beta + alpha_1) * C2_old_));
rDY_[36] = (((alpha_1 * C2_old_) + (calc_beta_2 * O_old_) + (calc_beta_ItoC2 * I_old_)) -
((beta_1 + calc_alpha_2 + calc_alpha_C2ToI) * C1_old_));
rDY_[37] = (((calc_alpha_2 * C1_old_) + (calc_beta_i * I_old_)) - ((calc_beta_2 + calc_alpha_i) * O_old_));
rDY_[38] = (((calc_alpha_C2ToI * C1_old_) + (calc_alpha_i * O_old_)) - ((calc_beta_ItoC2 + calc_beta_i) * I_old_));
rDY_[39] = ((calc_xs1ss - xs1_old_) / calc_txs1);
rDY_[40] = ((calc_xs2ss - xs2_old_) / calc_txs2);
rDY_[41] = ((calc_Jrel_inf - Jrel_np_old_) / calc_tau_rel);
rDY_[42] = ((calc_Jrel_infp - Jrel_p_old_) / calc_tau_relp);