/*
   There are a total of 226 entries in the algebraic variable array.
   There are a total of 45 entries in each of the rate and state variable arrays.
   There are a total of 162 entries in the constant variable array.
 */
/*
 * VOI is time in component environment (millisecond).
 * STATES[0] is v in component membrane (millivolt).
 * STATES[1] is CaMKt in component CaMK (millimolar).
 * STATES[2] is cass in component intracellular_ions (millimolar).
 * STATES[3] is nai in component intracellular_ions (millimolar).
 * STATES[4] is nass in component intracellular_ions (millimolar).
 * STATES[5] is ki in component intracellular_ions (millimolar).
 * STATES[6] is kss in component intracellular_ions (millimolar).
 * STATES[7] is cansr in component intracellular_ions (millimolar).
 * STATES[8] is cajsr in component intracellular_ions (millimolar).
 * STATES[9] is cai in component intracellular_ions (millimolar).
 * STATES[10] is cli in component intracellular_ions (millimolar).
 * STATES[11] is clss in component intracellular_ions (millimolar).
 * STATES[12] is m in component INa (dimensionless).
 * STATES[13] is h in component INa (dimensionless).
 * STATES[14] is j in component INa (dimensionless).
 * STATES[15] is hp in component INa (dimensionless).
 * STATES[16] is jp in component INa (dimensionless).
 * STATES[17] is mL in component INaL (dimensionless).
 * STATES[18] is hL in component INaL (dimensionless).
 * STATES[19] is hLp in component INaL (dimensionless).
 * STATES[20] is a in component Ito (dimensionless).
 * STATES[21] is iF in component Ito (dimensionless).
 * STATES[22] is iS in component Ito (dimensionless).
 * STATES[23] is ap in component Ito (dimensionless).
 * STATES[24] is iFp in component Ito (dimensionless).
 * STATES[25] is iSp in component Ito (dimensionless).
 * STATES[26] is d in component ICaL (dimensionless).
 * STATES[27] is ff in component ICaL (dimensionless).
 * STATES[28] is fs in component ICaL (dimensionless).
 * STATES[29] is fcaf in component ICaL (dimensionless).
 * STATES[30] is fcas in component ICaL (dimensionless).
 * STATES[31] is jca in component ICaL (dimensionless).
 * STATES[32] is ffp in component ICaL (dimensionless).
 * STATES[33] is fcafp in component ICaL (dimensionless).
 * STATES[34] is nca_ss in component ICaL (dimensionless).
 * STATES[35] is nca_i in component ICaL (dimensionless).
 * STATES[36] is C1 in component IKr (dimensionless).
 * STATES[37] is C2 in component IKr (dimensionless).
 * STATES[38] is C3 in component IKr (dimensionless).
 * STATES[39] is I in component IKr (dimensionless).
 * STATES[40] is O in component IKr (dimensionless).
 * STATES[41] is xs1 in component IKs (dimensionless).
 * STATES[42] is xs2 in component IKs (dimensionless).
 * STATES[43] is Jrel_np in component ryr (millimolar_per_millisecond).
 * STATES[44] is Jrel_p in component ryr (millimolar_per_millisecond).
 * CONSTANTS[0] is celltype in component environment (dimensionless).
 * CONSTANTS[1] is nao in component extracellular (millimolar).
 * CONSTANTS[2] is cao in component extracellular (millimolar).
 * CONSTANTS[3] is ko in component extracellular (millimolar).
 * CONSTANTS[4] is clo in component extracellular (millimolar).
 * CONSTANTS[5] is R in component physical_constants (joule_per_kilomole_kelvin).
 * CONSTANTS[6] is T in component physical_constants (kelvin).
 * CONSTANTS[7] is F in component physical_constants (coulomb_per_mole).
 * CONSTANTS[8] is zna in component physical_constants (dimensionless).
 * CONSTANTS[9] is zca in component physical_constants (dimensionless).
 * CONSTANTS[10] is zk in component physical_constants (dimensionless).
 * CONSTANTS[11] is zcl in component physical_constants (dimensionless).
 * CONSTANTS[12] is L in component cell_geometry (centimeter).
 * CONSTANTS[13] is rad in component cell_geometry (centimeter).
 * CONSTANTS[112] is vcell in component cell_geometry (microliter).
 * CONSTANTS[129] is Ageo in component cell_geometry (centimeter_squared).
 * CONSTANTS[135] is Acap in component cell_geometry (centimeter_squared).
 * CONSTANTS[141] is vmyo in component cell_geometry (microliter).
 * CONSTANTS[142] is vnsr in component cell_geometry (microliter).
 * CONSTANTS[143] is vjsr in component cell_geometry (microliter).
 * CONSTANTS[144] is vss in component cell_geometry (microliter).
 * ALGEBRAIC[25] is vffrt in component membrane (coulomb_per_mole).
 * ALGEBRAIC[28] is vfrt in component membrane (dimensionless).
 * ALGEBRAIC[70] is INa in component INa (microA_per_microF).
 * ALGEBRAIC[72] is INaL in component INaL (microA_per_microF).
 * ALGEBRAIC[78] is Ito in component Ito (microA_per_microF).
 * ALGEBRAIC[114] is ICaL in component ICaL (microA_per_microF).
 * ALGEBRAIC[115] is ICaNa in component ICaL (microA_per_microF).
 * ALGEBRAIC[116] is ICaK in component ICaL (microA_per_microF).
 * ALGEBRAIC[117] is IKr in component IKr (microA_per_microF).
 * ALGEBRAIC[119] is IKs in component IKs (microA_per_microF).
 * ALGEBRAIC[123] is IK1 in component IK1 (microA_per_microF).
 * ALGEBRAIC[155] is INaCa_i in component INaCa (microA_per_microF).
 * ALGEBRAIC[185] is INaCa_ss in component INaCa (microA_per_microF).
 * ALGEBRAIC[204] is INaK in component INaK (microA_per_microF).
 * ALGEBRAIC[207] is INab in component INab (microA_per_microF).
 * ALGEBRAIC[206] is IKb in component IKb (microA_per_microF).
 * ALGEBRAIC[211] is IpCa in component IpCa (microA_per_microF).
 * ALGEBRAIC[209] is ICab in component ICab (microA_per_microF).
 * ALGEBRAIC[216] is IClCa in component ICl (microA_per_microF).
 * ALGEBRAIC[218] is IClb in component ICl (microA_per_microF).
 * ALGEBRAIC[68] is I_katp in component I_katp (microA_per_microF).
 * ALGEBRAIC[11] is Istim in component membrane (microA_per_microF).
 * CONSTANTS[14] is i_Stim_Start in component membrane (millisecond).
 * CONSTANTS[15] is i_Stim_End in component membrane (millisecond).
 * CONSTANTS[16] is i_Stim_Amplitude in component membrane (microA_per_microF).
 * CONSTANTS[17] is i_Stim_Period in component membrane (millisecond).
 * CONSTANTS[18] is i_Stim_PulseDuration in component membrane (millisecond).
 * CONSTANTS[19] is KmCaMK in component CaMK (millimolar).
 * CONSTANTS[20] is aCaMK in component CaMK (per_millimolar_per_millisecond).
 * CONSTANTS[21] is bCaMK in component CaMK (per_millisecond).
 * CONSTANTS[22] is CaMKo in component CaMK (dimensionless).
 * CONSTANTS[23] is KmCaM in component CaMK (millimolar).
 * ALGEBRAIC[43] is CaMKb in component CaMK (millimolar).
 * ALGEBRAIC[49] is CaMKa in component CaMK (millimolar).
 * CONSTANTS[24] is cmdnmax_b in component intracellular_ions (millimolar).
 * CONSTANTS[113] is cmdnmax in component intracellular_ions (millimolar).
 * CONSTANTS[25] is kmcmdn in component intracellular_ions (millimolar).
 * CONSTANTS[26] is trpnmax in component intracellular_ions (millimolar).
 * CONSTANTS[27] is kmtrpn in component intracellular_ions (millimolar).
 * CONSTANTS[28] is BSRmax in component intracellular_ions (millimolar).
 * CONSTANTS[29] is KmBSR in component intracellular_ions (millimolar).
 * CONSTANTS[30] is BSLmax in component intracellular_ions (millimolar).
 * CONSTANTS[31] is KmBSL in component intracellular_ions (millimolar).
 * CONSTANTS[32] is csqnmax in component intracellular_ions (millimolar).
 * CONSTANTS[33] is kmcsqn in component intracellular_ions (millimolar).
 * ALGEBRAIC[93] is ICaL_ss in component ICaL (microA_per_microF).
 * ALGEBRAIC[94] is ICaNa_ss in component ICaL (microA_per_microF).
 * ALGEBRAIC[97] is ICaK_ss in component ICaL (microA_per_microF).
 * ALGEBRAIC[111] is ICaL_i in component ICaL (microA_per_microF).
 * ALGEBRAIC[112] is ICaNa_i in component ICaL (microA_per_microF).
 * ALGEBRAIC[113] is ICaK_i in component ICaL (microA_per_microF).
 * ALGEBRAIC[210] is JdiffNa in component diff (millimolar_per_millisecond).
 * ALGEBRAIC[213] is Jdiff in component diff (millimolar_per_millisecond).
 * ALGEBRAIC[224] is Jup in component SERCA (millimolar_per_millisecond).
 * ALGEBRAIC[208] is JdiffK in component diff (millimolar_per_millisecond).
 * ALGEBRAIC[221] is JdiffCl in component diff (millimolar_per_millisecond).
 * ALGEBRAIC[217] is Jrel in component ryr (millimolar_per_millisecond).
 * ALGEBRAIC[225] is Jtr in component trans_flux (millimolar_per_millisecond).
 * ALGEBRAIC[53] is Bcai in component intracellular_ions (dimensionless).
 * ALGEBRAIC[59] is Bcajsr in component intracellular_ions (dimensionless).
 * ALGEBRAIC[56] is Bcass in component intracellular_ions (dimensionless).
 * ALGEBRAIC[214] is IClCa_sl in component ICl (microA_per_microF).
 * ALGEBRAIC[212] is IClCa_junc in component ICl (microA_per_microF).
 * CONSTANTS[34] is PKNa in component reversal_potentials (dimensionless).
 * ALGEBRAIC[63] is ENa in component reversal_potentials (millivolt).
 * ALGEBRAIC[64] is EK in component reversal_potentials (millivolt).
 * ALGEBRAIC[65] is EKs in component reversal_potentials (millivolt).
 * ALGEBRAIC[66] is ECl in component reversal_potentials (millivolt).
 * ALGEBRAIC[67] is EClss in component reversal_potentials (millivolt).
 * CONSTANTS[35] is gkatp in component I_katp (milliS_per_microF).
 * CONSTANTS[36] is fkatp in component I_katp (dimensionless).
 * CONSTANTS[37] is K_o_n in component I_katp (millimolar).
 * CONSTANTS[38] is A_atp in component I_katp (millimolar).
 * CONSTANTS[39] is K_atp in component I_katp (millimolar).
 * CONSTANTS[114] is akik in component I_katp (dimensionless).
 * CONSTANTS[115] is bkik in component I_katp (dimensionless).
 * ALGEBRAIC[0] is mss in component INa (dimensionless).
 * ALGEBRAIC[13] is tm in component INa (millisecond).
 * ALGEBRAIC[1] is hss in component INa (dimensionless).
 * ALGEBRAIC[14] is ah in component INa (dimensionless).
 * ALGEBRAIC[29] is bh in component INa (dimensionless).
 * ALGEBRAIC[37] is th in component INa (millisecond).
 * ALGEBRAIC[38] is jss in component INa (dimensionless).
 * ALGEBRAIC[15] is aj in component INa (dimensionless).
 * ALGEBRAIC[30] is bj in component INa (dimensionless).
 * ALGEBRAIC[44] is tj in component INa (millisecond).
 * ALGEBRAIC[45] is hssp in component INa (dimensionless).
 * ALGEBRAIC[50] is tjp in component INa (millisecond).
 * ALGEBRAIC[69] is fINap in component INa (dimensionless).
 * CONSTANTS[40] is GNa in component INa (milliS_per_microF).
 * ALGEBRAIC[2] is mLss in component INaL (dimensionless).
 * ALGEBRAIC[16] is tmL in component INaL (millisecond).
 * CONSTANTS[41] is thL in component INaL (millisecond).
 * ALGEBRAIC[3] is hLss in component INaL (dimensionless).
 * ALGEBRAIC[4] is hLssp in component INaL (dimensionless).
 * CONSTANTS[116] is thLp in component INaL (millisecond).
 * CONSTANTS[42] is GNaL_b in component INaL (milliS_per_microF).
 * CONSTANTS[117] is GNaL in component INaL (milliS_per_microF).
 * ALGEBRAIC[71] is fINaLp in component INaL (dimensionless).
 * CONSTANTS[43] is Gto_b in component Ito (milliS_per_microF).
 * ALGEBRAIC[5] is ass in component Ito (dimensionless).
 * ALGEBRAIC[17] is ta in component Ito (millisecond).
 * CONSTANTS[44] is EKshift in component Ito (millivolt).
 * ALGEBRAIC[6] is iss in component Ito (dimensionless).
 * ALGEBRAIC[18] is delta_epi in component Ito (dimensionless).
 * ALGEBRAIC[31] is tiF_b in component Ito (millisecond).
 * ALGEBRAIC[39] is tiS_b in component Ito (millisecond).
 * ALGEBRAIC[46] is tiF in component Ito (millisecond).
 * ALGEBRAIC[51] is tiS in component Ito (millisecond).
 * ALGEBRAIC[73] is AiF in component Ito (dimensionless).
 * ALGEBRAIC[74] is AiS in component Ito (dimensionless).
 * ALGEBRAIC[75] is i in component Ito (dimensionless).
 * ALGEBRAIC[32] is assp in component Ito (dimensionless).
 * ALGEBRAIC[54] is dti_develop in component Ito (dimensionless).
 * ALGEBRAIC[57] is dti_recover in component Ito (dimensionless).
 * ALGEBRAIC[60] is tiFp in component Ito (millisecond).
 * ALGEBRAIC[61] is tiSp in component Ito (millisecond).
 * ALGEBRAIC[76] is ip in component Ito (dimensionless).
 * CONSTANTS[118] is Gto in component Ito (milliS_per_microF).
 * ALGEBRAIC[77] is fItop in component Ito (dimensionless).
 * CONSTANTS[45] is Kmn in component ICaL (millimolar).
 * CONSTANTS[46] is k2n in component ICaL (per_millisecond).
 * CONSTANTS[47] is PCa_b in component ICaL (dimensionless).
 * ALGEBRAIC[7] is dss in component ICaL (dimensionless).
 * ALGEBRAIC[8] is fss in component ICaL (dimensionless).
 * CONSTANTS[48] is Aff in component ICaL (dimensionless).
 * CONSTANTS[119] is Afs in component ICaL (dimensionless).
 * ALGEBRAIC[79] is f in component ICaL (dimensionless).
 * ALGEBRAIC[19] is fcass in component ICaL (dimensionless).
 * ALGEBRAIC[9] is jcass in component ICaL (dimensionless).
 * ALGEBRAIC[80] is Afcaf in component ICaL (dimensionless).
 * ALGEBRAIC[81] is Afcas in component ICaL (dimensionless).
 * ALGEBRAIC[82] is fca in component ICaL (dimensionless).
 * ALGEBRAIC[83] is fp in component ICaL (dimensionless).
 * ALGEBRAIC[84] is fcap in component ICaL (dimensionless).
 * ALGEBRAIC[10] is km2n in component ICaL (per_millisecond).
 * ALGEBRAIC[20] is anca_ss in component ICaL (dimensionless).
 * ALGEBRAIC[21] is anca_i in component ICaL (dimensionless).
 * ALGEBRAIC[89] is PhiCaL_ss in component ICaL (dimensionless).
 * ALGEBRAIC[90] is PhiCaNa_ss in component ICaL (dimensionless).
 * ALGEBRAIC[91] is PhiCaK_ss in component ICaL (dimensionless).
 * ALGEBRAIC[108] is PhiCaL_i in component ICaL (dimensionless).
 * ALGEBRAIC[109] is PhiCaNa_i in component ICaL (dimensionless).
 * ALGEBRAIC[110] is PhiCaK_i in component ICaL (dimensionless).
 * CONSTANTS[120] is PCa in component ICaL (dimensionless).
 * CONSTANTS[130] is PCap in component ICaL (dimensionless).
 * CONSTANTS[131] is PCaNa in component ICaL (dimensionless).
 * CONSTANTS[132] is PCaK in component ICaL (dimensionless).
 * CONSTANTS[136] is PCaNap in component ICaL (dimensionless).
 * CONSTANTS[137] is PCaKp in component ICaL (dimensionless).
 * ALGEBRAIC[92] is fICaLp in component ICaL (dimensionless).
 * ALGEBRAIC[22] is td in component ICaL (millisecond).
 * ALGEBRAIC[23] is tff in component ICaL (millisecond).
 * ALGEBRAIC[24] is tfs in component ICaL (millisecond).
 * ALGEBRAIC[33] is tfcaf in component ICaL (millisecond).
 * ALGEBRAIC[34] is tfcas in component ICaL (millisecond).
 * CONSTANTS[49] is tjca in component ICaL (millisecond).
 * ALGEBRAIC[35] is tffp in component ICaL (millisecond).
 * ALGEBRAIC[40] is tfcafp in component ICaL (millisecond).
 * CONSTANTS[50] is vShift in component ICaL (millivolt).
 * CONSTANTS[51] is offset in component ICaL (millisecond).
 * CONSTANTS[121] is Io in component ICaL (dimensionless).
 * ALGEBRAIC[85] is Iss in component ICaL (dimensionless).
 * ALGEBRAIC[100] is Ii in component ICaL (dimensionless).
 * CONSTANTS[52] is dielConstant in component ICaL (per_kelvin).
 * CONSTANTS[133] is constA in component ICaL (dimensionless).
 * CONSTANTS[138] is gamma_cao in component ICaL (dimensionless).
 * ALGEBRAIC[86] is gamma_cass in component ICaL (dimensionless).
 * ALGEBRAIC[103] is gamma_cai in component ICaL (dimensionless).
 * CONSTANTS[139] is gamma_nao in component ICaL (dimensionless).
 * ALGEBRAIC[87] is gamma_nass in component ICaL (dimensionless).
 * ALGEBRAIC[106] is gamma_nai in component ICaL (dimensionless).
 * CONSTANTS[140] is gamma_ko in component ICaL (dimensionless).
 * ALGEBRAIC[88] is gamma_kss in component ICaL (dimensionless).
 * ALGEBRAIC[107] is gamma_ki in component ICaL (dimensionless).
 * CONSTANTS[53] is ICaL_fractionSS in component ICaL (dimensionless).
 * CONSTANTS[54] is GKr_b in component IKr (milliS_per_microF).
 * ALGEBRAIC[41] is alpha in component IKr (per_millisecond).
 * ALGEBRAIC[47] is beta in component IKr (per_millisecond).
 * CONSTANTS[55] is alpha_1 in component IKr (per_millisecond).
 * CONSTANTS[56] is beta_1 in component IKr (per_millisecond).
 * ALGEBRAIC[42] is alpha_2 in component IKr (per_millisecond).
 * ALGEBRAIC[48] is beta_2 in component IKr (per_millisecond).
 * ALGEBRAIC[52] is alpha_i in component IKr (per_millisecond).
 * ALGEBRAIC[55] is beta_i in component IKr (per_millisecond).
 * ALGEBRAIC[58] is alpha_C2ToI in component IKr (per_millisecond).
 * ALGEBRAIC[62] is beta_ItoC2 in component IKr (per_millisecond).
 * CONSTANTS[122] is GKr in component IKr (milliS_per_microF).
 * CONSTANTS[57] is GKs_b in component IKs (milliS_per_microF).
 * CONSTANTS[123] is GKs in component IKs (milliS_per_microF).
 * ALGEBRAIC[12] is xs1ss in component IKs (dimensionless).
 * ALGEBRAIC[26] is xs2ss in component IKs (dimensionless).
 * ALGEBRAIC[27] is txs1 in component IKs (millisecond).
 * ALGEBRAIC[118] is KsCa in component IKs (dimensionless).
 * ALGEBRAIC[36] is txs2 in component IKs (millisecond).
 * CONSTANTS[124] is GK1 in component IK1 (milliS_per_microF).
 * CONSTANTS[58] is GK1_b in component IK1 (milliS_per_microF).
 * ALGEBRAIC[120] is aK1 in component IK1 (dimensionless).
 * ALGEBRAIC[121] is bK1 in component IK1 (dimensionless).
 * ALGEBRAIC[122] is K1ss in component IK1 (dimensionless).
 * CONSTANTS[59] is INaCa_fractionSS in component INaCa (dimensionless).
 * CONSTANTS[60] is kna1 in component INaCa (per_millisecond).
 * CONSTANTS[61] is kna2 in component INaCa (per_millisecond).
 * CONSTANTS[62] is kna3 in component INaCa (per_millisecond).
 * CONSTANTS[63] is kasymm in component INaCa (dimensionless).
 * CONSTANTS[64] is wna in component INaCa (dimensionless).
 * CONSTANTS[65] is wca in component INaCa (dimensionless).
 * CONSTANTS[66] is wnaca in component INaCa (dimensionless).
 * CONSTANTS[67] is kcaon in component INaCa (per_millisecond).
 * CONSTANTS[68] is kcaoff in component INaCa (per_millisecond).
 * CONSTANTS[69] is qna in component INaCa (dimensionless).
 * CONSTANTS[70] is qca in component INaCa (dimensionless).
 * ALGEBRAIC[125] is hna in component INaCa (dimensionless).
 * ALGEBRAIC[124] is hca in component INaCa (dimensionless).
 * CONSTANTS[71] is KmCaAct in component INaCa (millimolar).
 * CONSTANTS[72] is Gncx_b in component INaCa (milliS_per_microF).
 * CONSTANTS[151] is Gncx in component INaCa (milliS_per_microF).
 * ALGEBRAIC[126] is h1_i in component INaCa (dimensionless).
 * ALGEBRAIC[127] is h2_i in component INaCa (dimensionless).
 * ALGEBRAIC[128] is h3_i in component INaCa (dimensionless).
 * ALGEBRAIC[129] is h4_i in component INaCa (dimensionless).
 * ALGEBRAIC[130] is h5_i in component INaCa (dimensionless).
 * ALGEBRAIC[131] is h6_i in component INaCa (dimensionless).
 * ALGEBRAIC[132] is h7_i in component INaCa (dimensionless).
 * ALGEBRAIC[133] is h8_i in component INaCa (dimensionless).
 * ALGEBRAIC[134] is h9_i in component INaCa (dimensionless).
 * CONSTANTS[145] is h10_i in component INaCa (dimensionless).
 * CONSTANTS[146] is h11_i in component INaCa (dimensionless).
 * CONSTANTS[147] is h12_i in component INaCa (dimensionless).
 * CONSTANTS[148] is k1_i in component INaCa (dimensionless).
 * CONSTANTS[149] is k2_i in component INaCa (dimensionless).
 * ALGEBRAIC[135] is k3p_i in component INaCa (dimensionless).
 * ALGEBRAIC[136] is k3pp_i in component INaCa (dimensionless).
 * ALGEBRAIC[137] is k3_i in component INaCa (dimensionless).
 * ALGEBRAIC[140] is k4_i in component INaCa (dimensionless).
 * ALGEBRAIC[138] is k4p_i in component INaCa (dimensionless).
 * ALGEBRAIC[139] is k4pp_i in component INaCa (dimensionless).
 * CONSTANTS[150] is k5_i in component INaCa (dimensionless).
 * ALGEBRAIC[141] is k6_i in component INaCa (dimensionless).
 * ALGEBRAIC[142] is k7_i in component INaCa (dimensionless).
 * ALGEBRAIC[143] is k8_i in component INaCa (dimensionless).
 * ALGEBRAIC[144] is x1_i in component INaCa (dimensionless).
 * ALGEBRAIC[145] is x2_i in component INaCa (dimensionless).
 * ALGEBRAIC[146] is x3_i in component INaCa (dimensionless).
 * ALGEBRAIC[147] is x4_i in component INaCa (dimensionless).
 * ALGEBRAIC[148] is E1_i in component INaCa (dimensionless).
 * ALGEBRAIC[149] is E2_i in component INaCa (dimensionless).
 * ALGEBRAIC[150] is E3_i in component INaCa (dimensionless).
 * ALGEBRAIC[151] is E4_i in component INaCa (dimensionless).
 * ALGEBRAIC[152] is allo_i in component INaCa (dimensionless).
 * ALGEBRAIC[153] is JncxNa_i in component INaCa (millimolar_per_millisecond).
 * ALGEBRAIC[154] is JncxCa_i in component INaCa (millimolar_per_millisecond).
 * ALGEBRAIC[156] is h1_ss in component INaCa (dimensionless).
 * ALGEBRAIC[157] is h2_ss in component INaCa (dimensionless).
 * ALGEBRAIC[158] is h3_ss in component INaCa (dimensionless).
 * ALGEBRAIC[159] is h4_ss in component INaCa (dimensionless).
 * ALGEBRAIC[160] is h5_ss in component INaCa (dimensionless).
 * ALGEBRAIC[161] is h6_ss in component INaCa (dimensionless).
 * ALGEBRAIC[162] is h7_ss in component INaCa (dimensionless).
 * ALGEBRAIC[163] is h8_ss in component INaCa (dimensionless).
 * ALGEBRAIC[164] is h9_ss in component INaCa (dimensionless).
 * CONSTANTS[152] is h10_ss in component INaCa (dimensionless).
 * CONSTANTS[153] is h11_ss in component INaCa (dimensionless).
 * CONSTANTS[154] is h12_ss in component INaCa (dimensionless).
 * CONSTANTS[155] is k1_ss in component INaCa (dimensionless).
 * CONSTANTS[156] is k2_ss in component INaCa (dimensionless).
 * ALGEBRAIC[165] is k3p_ss in component INaCa (dimensionless).
 * ALGEBRAIC[166] is k3pp_ss in component INaCa (dimensionless).
 * ALGEBRAIC[167] is k3_ss in component INaCa (dimensionless).
 * ALGEBRAIC[170] is k4_ss in component INaCa (dimensionless).
 * ALGEBRAIC[168] is k4p_ss in component INaCa (dimensionless).
 * ALGEBRAIC[169] is k4pp_ss in component INaCa (dimensionless).
 * CONSTANTS[157] is k5_ss in component INaCa (dimensionless).
 * ALGEBRAIC[171] is k6_ss in component INaCa (dimensionless).
 * ALGEBRAIC[172] is k7_ss in component INaCa (dimensionless).
 * ALGEBRAIC[173] is k8_ss in component INaCa (dimensionless).
 * ALGEBRAIC[174] is x1_ss in component INaCa (dimensionless).
 * ALGEBRAIC[175] is x2_ss in component INaCa (dimensionless).
 * ALGEBRAIC[176] is x3_ss in component INaCa (dimensionless).
 * ALGEBRAIC[177] is x4_ss in component INaCa (dimensionless).
 * ALGEBRAIC[178] is E1_ss in component INaCa (dimensionless).
 * ALGEBRAIC[179] is E2_ss in component INaCa (dimensionless).
 * ALGEBRAIC[180] is E3_ss in component INaCa (dimensionless).
 * ALGEBRAIC[181] is E4_ss in component INaCa (dimensionless).
 * ALGEBRAIC[182] is allo_ss in component INaCa (dimensionless).
 * ALGEBRAIC[183] is JncxNa_ss in component INaCa (millimolar_per_millisecond).
 * ALGEBRAIC[184] is JncxCa_ss in component INaCa (millimolar_per_millisecond).
 * CONSTANTS[73] is k1p in component INaK (per_millisecond).
 * CONSTANTS[74] is k1m in component INaK (per_millisecond).
 * CONSTANTS[75] is k2p in component INaK (per_millisecond).
 * CONSTANTS[76] is k2m in component INaK (per_millisecond).
 * CONSTANTS[77] is k3p in component INaK (per_millisecond).
 * CONSTANTS[78] is k3m in component INaK (per_millisecond).
 * CONSTANTS[79] is k4p in component INaK (per_millisecond).
 * CONSTANTS[80] is k4m in component INaK (per_millisecond).
 * CONSTANTS[81] is Knai0 in component INaK (millimolar).
 * CONSTANTS[82] is Knao0 in component INaK (millimolar).
 * CONSTANTS[83] is delta in component INaK (millivolt).
 * CONSTANTS[84] is Kki in component INaK (per_millisecond).
 * CONSTANTS[85] is Kko in component INaK (per_millisecond).
 * CONSTANTS[86] is MgADP in component INaK (millimolar).
 * CONSTANTS[87] is MgATP in component INaK (millimolar).
 * CONSTANTS[88] is Kmgatp in component INaK (millimolar).
 * CONSTANTS[89] is H in component INaK (millimolar).
 * CONSTANTS[90] is eP in component INaK (dimensionless).
 * CONSTANTS[91] is Khp in component INaK (millimolar).
 * CONSTANTS[92] is Knap in component INaK (millimolar).
 * CONSTANTS[93] is Kxkur in component INaK (millimolar).
 * CONSTANTS[94] is Pnak_b in component INaK (milliS_per_microF).
 * CONSTANTS[161] is Pnak in component INaK (milliS_per_microF).
 * ALGEBRAIC[186] is Knai in component INaK (millimolar).
 * ALGEBRAIC[187] is Knao in component INaK (millimolar).
 * ALGEBRAIC[188] is P in component INaK (dimensionless).
 * ALGEBRAIC[189] is a1 in component INaK (dimensionless).
 * CONSTANTS[158] is b1 in component INaK (dimensionless).
 * CONSTANTS[159] is a2 in component INaK (dimensionless).
 * ALGEBRAIC[190] is b2 in component INaK (dimensionless).
 * ALGEBRAIC[191] is a3 in component INaK (dimensionless).
 * ALGEBRAIC[192] is b3 in component INaK (dimensionless).
 * CONSTANTS[160] is a4 in component INaK (dimensionless).
 * ALGEBRAIC[193] is b4 in component INaK (dimensionless).
 * ALGEBRAIC[194] is x1 in component INaK (dimensionless).
 * ALGEBRAIC[195] is x2 in component INaK (dimensionless).
 * ALGEBRAIC[196] is x3 in component INaK (dimensionless).
 * ALGEBRAIC[197] is x4 in component INaK (dimensionless).
 * ALGEBRAIC[198] is E1 in component INaK (dimensionless).
 * ALGEBRAIC[199] is E2 in component INaK (dimensionless).
 * ALGEBRAIC[200] is E3 in component INaK (dimensionless).
 * ALGEBRAIC[201] is E4 in component INaK (dimensionless).
 * ALGEBRAIC[202] is JnakNa in component INaK (millimolar_per_millisecond).
 * ALGEBRAIC[203] is JnakK in component INaK (millimolar_per_millisecond).
 * ALGEBRAIC[205] is xkb in component IKb (dimensionless).
 * CONSTANTS[95] is GKb_b in component IKb (milliS_per_microF).
 * CONSTANTS[125] is GKb in component IKb (milliS_per_microF).
 * CONSTANTS[96] is PNab in component INab (milliS_per_microF).
 * CONSTANTS[97] is PCab in component ICab (milliS_per_microF).
 * CONSTANTS[98] is GpCa in component IpCa (milliS_per_microF).
 * CONSTANTS[99] is KmCap in component IpCa (millimolar).
 * CONSTANTS[100] is GClCa in component ICl (milliS_per_microF).
 * CONSTANTS[101] is GClb in component ICl (milliS_per_microF).
 * CONSTANTS[102] is KdClCa in component ICl (millimolar).
 * CONSTANTS[103] is Fjunc in component ICl (dimensionless).
 * CONSTANTS[104] is tauNa in component diff (millisecond).
 * CONSTANTS[105] is tauK in component diff (millisecond).
 * CONSTANTS[106] is tauCa in component diff (millisecond).
 * CONSTANTS[107] is tauCl in component diff (millisecond).
 * CONSTANTS[108] is bt in component ryr (millisecond).
 * CONSTANTS[126] is a_rel in component ryr (millimolar_per_millisecond).
 * ALGEBRAIC[95] is Jrel_inf_b in component ryr (millimolar_per_millisecond).
 * ALGEBRAIC[98] is Jrel_inf in component ryr (millimolar_per_millisecond).
 * ALGEBRAIC[101] is tau_rel_b in component ryr (millisecond).
 * ALGEBRAIC[104] is tau_rel in component ryr (millisecond).
 * CONSTANTS[127] is btp in component ryr (millisecond).
 * CONSTANTS[134] is a_relp in component ryr (millimolar_per_millisecond).
 * ALGEBRAIC[96] is Jrel_infp_b in component ryr (millimolar_per_millisecond).
 * ALGEBRAIC[99] is Jrel_infp in component ryr (millimolar_per_millisecond).
 * ALGEBRAIC[102] is tau_relp_b in component ryr (millisecond).
 * ALGEBRAIC[105] is tau_relp in component ryr (millisecond).
 * CONSTANTS[109] is cajsr_half in component ryr (millimolar).
 * ALGEBRAIC[215] is fJrelp in component ryr (dimensionless).
 * CONSTANTS[110] is Jrel_b in component ryr (dimensionless).
 * CONSTANTS[128] is upScale in component SERCA (dimensionless).
 * ALGEBRAIC[219] is Jupnp in component SERCA (millimolar_per_millisecond).
 * ALGEBRAIC[220] is Jupp in component SERCA (millimolar_per_millisecond).
 * ALGEBRAIC[222] is fJupp in component SERCA (dimensionless).
 * ALGEBRAIC[223] is Jleak in component SERCA (millimolar_per_millisecond).
 * CONSTANTS[111] is Jup_b in component SERCA (dimensionless).
 * RATES[0] is d/dt v in component membrane (millivolt).
 * RATES[1] is d/dt CaMKt in component CaMK (millimolar).
 * RATES[3] is d/dt nai in component intracellular_ions (millimolar).
 * RATES[4] is d/dt nass in component intracellular_ions (millimolar).
 * RATES[5] is d/dt ki in component intracellular_ions (millimolar).
 * RATES[6] is d/dt kss in component intracellular_ions (millimolar).
 * RATES[10] is d/dt cli in component intracellular_ions (millimolar).
 * RATES[11] is d/dt clss in component intracellular_ions (millimolar).
 * RATES[9] is d/dt cai in component intracellular_ions (millimolar).
 * RATES[2] is d/dt cass in component intracellular_ions (millimolar).
 * RATES[7] is d/dt cansr in component intracellular_ions (millimolar).
 * RATES[8] is d/dt cajsr in component intracellular_ions (millimolar).
 * RATES[12] is d/dt m in component INa (dimensionless).
 * RATES[13] is d/dt h in component INa (dimensionless).
 * RATES[14] is d/dt j in component INa (dimensionless).
 * RATES[15] is d/dt hp in component INa (dimensionless).
 * RATES[16] is d/dt jp in component INa (dimensionless).
 * RATES[17] is d/dt mL in component INaL (dimensionless).
 * RATES[18] is d/dt hL in component INaL (dimensionless).
 * RATES[19] is d/dt hLp in component INaL (dimensionless).
 * RATES[20] is d/dt a in component Ito (dimensionless).
 * RATES[21] is d/dt iF in component Ito (dimensionless).
 * RATES[22] is d/dt iS in component Ito (dimensionless).
 * RATES[23] is d/dt ap in component Ito (dimensionless).
 * RATES[24] is d/dt iFp in component Ito (dimensionless).
 * RATES[25] is d/dt iSp in component Ito (dimensionless).
 * RATES[26] is d/dt d in component ICaL (dimensionless).
 * RATES[27] is d/dt ff in component ICaL (dimensionless).
 * RATES[28] is d/dt fs in component ICaL (dimensionless).
 * RATES[29] is d/dt fcaf in component ICaL (dimensionless).
 * RATES[30] is d/dt fcas in component ICaL (dimensionless).
 * RATES[31] is d/dt jca in component ICaL (dimensionless).
 * RATES[32] is d/dt ffp in component ICaL (dimensionless).
 * RATES[33] is d/dt fcafp in component ICaL (dimensionless).
 * RATES[34] is d/dt nca_ss in component ICaL (dimensionless).
 * RATES[35] is d/dt nca_i in component ICaL (dimensionless).
 * RATES[38] is d/dt C3 in component IKr (dimensionless).
 * RATES[37] is d/dt C2 in component IKr (dimensionless).
 * RATES[36] is d/dt C1 in component IKr (dimensionless).
 * RATES[40] is d/dt O in component IKr (dimensionless).
 * RATES[39] is d/dt I in component IKr (dimensionless).
 * RATES[41] is d/dt xs1 in component IKs (dimensionless).
 * RATES[42] is d/dt xs2 in component IKs (dimensionless).
 * RATES[43] is d/dt Jrel_np in component ryr (millimolar_per_millisecond).
 * RATES[44] is d/dt Jrel_p in component ryr (millimolar_per_millisecond).
 */
void
initConsts(double* CONSTANTS, double* RATES, double *STATES)
{
CONSTANTS[0] = 0;
CONSTANTS[1] = 140.0;
CONSTANTS[2] = 1.8;
CONSTANTS[3] = 5.0;
CONSTANTS[4] = 150.0;
CONSTANTS[5] = 8314;
CONSTANTS[6] = 310;
CONSTANTS[7] = 96485;
CONSTANTS[8] = 1;
CONSTANTS[9] = 2;
CONSTANTS[10] = 1;
CONSTANTS[11] = -1;
CONSTANTS[12] = 0.01;
CONSTANTS[13] = 0.0011;
CONSTANTS[14] = 0;
CONSTANTS[15] = 100000000000000000;
CONSTANTS[16] = -53;
CONSTANTS[17] = 1000;
CONSTANTS[18] = 1.0;
CONSTANTS[19] = 0.15;
CONSTANTS[20] = 0.05;
CONSTANTS[21] = 0.00068;
CONSTANTS[22] = 0.05;
CONSTANTS[23] = 0.0015;
CONSTANTS[24] = 0.05;
CONSTANTS[25] = 0.00238;
CONSTANTS[26] = 0.07;
CONSTANTS[27] = 0.0005;
CONSTANTS[28] = 0.047;
CONSTANTS[29] = 0.00087;
CONSTANTS[30] = 1.124;
CONSTANTS[31] = 0.0087;
CONSTANTS[32] = 10;
CONSTANTS[33] = 0.8;
CONSTANTS[34] = 0.01833;
CONSTANTS[35] = 4.3195;
CONSTANTS[36] = 0.0;
CONSTANTS[37] = 5;
CONSTANTS[38] = 2;
CONSTANTS[39] = 0.25;
CONSTANTS[40] = 11.7802;
CONSTANTS[41] = 200;
CONSTANTS[42] = 0.0279;
CONSTANTS[43] = 0.16;
CONSTANTS[44] = 0;
CONSTANTS[45] = 0.002;
CONSTANTS[46] = 500;
CONSTANTS[47] = 8.3757e-05;
CONSTANTS[48] = 0.6;
CONSTANTS[49] = 72.5;
CONSTANTS[50] = 0;
CONSTANTS[51] = 0;
CONSTANTS[52] = 74;
CONSTANTS[53] = 0.8;
CONSTANTS[54] = 0.0321;
CONSTANTS[55] = 0.154375;
CONSTANTS[56] = 0.1911;
CONSTANTS[57] = 0.0011;
CONSTANTS[58] = 0.6992;
CONSTANTS[59] = 0.35;
CONSTANTS[60] = 15;
CONSTANTS[61] = 5;
CONSTANTS[62] = 88.12;
CONSTANTS[63] = 12.5;
CONSTANTS[64] = 6e4;
CONSTANTS[65] = 6e4;
CONSTANTS[66] = 5e3;
CONSTANTS[67] = 1.5e6;
CONSTANTS[68] = 5e3;
CONSTANTS[69] = 0.5224;
CONSTANTS[70] = 0.167;
CONSTANTS[71] = 150e-6;
CONSTANTS[72] = 0.0034;
CONSTANTS[73] = 949.5;
CONSTANTS[74] = 182.4;
CONSTANTS[75] = 687.2;
CONSTANTS[76] = 39.4;
CONSTANTS[77] = 1899;
CONSTANTS[78] = 79300;
CONSTANTS[79] = 639;
CONSTANTS[80] = 40;
CONSTANTS[81] = 9.073;
CONSTANTS[82] = 27.78;
CONSTANTS[83] = -0.155;
CONSTANTS[84] = 0.5;
CONSTANTS[85] = 0.3582;
CONSTANTS[86] = 0.05;
CONSTANTS[87] = 9.8;
CONSTANTS[88] = 1.698e-7;
CONSTANTS[89] = 1e-7;
CONSTANTS[90] = 4.2;
CONSTANTS[91] = 1.698e-7;
CONSTANTS[92] = 224;
CONSTANTS[93] = 292;
CONSTANTS[94] = 15.4509;
CONSTANTS[95] = 0.0189;
CONSTANTS[96] = 1.9239e-09;
CONSTANTS[97] = 5.9194e-08;
CONSTANTS[98] = 5e-04;
CONSTANTS[99] = 0.0005;
CONSTANTS[100] = 0.2843;
CONSTANTS[101] = 1.98e-3;
CONSTANTS[102] = 0.1;
CONSTANTS[103] = 1;
CONSTANTS[104] = 2.0;
CONSTANTS[105] = 2.0;
CONSTANTS[106] = 0.2;
CONSTANTS[107] = 2.0;
CONSTANTS[108] = 4.75;
CONSTANTS[109] = 1.7;
CONSTANTS[110] = 1.5378;
CONSTANTS[111] = 1.0;
CONSTANTS[112] =  1000.00*3.14000*CONSTANTS[13]*CONSTANTS[13]*CONSTANTS[12];
CONSTANTS[113] = (CONSTANTS[0]==1.00000 ?  CONSTANTS[24]*1.30000 : CONSTANTS[24]);
CONSTANTS[114] = pow(CONSTANTS[3]/CONSTANTS[37], 0.240000);
CONSTANTS[115] = 1.00000/(1.00000+pow(CONSTANTS[38]/CONSTANTS[39], 2.00000));
CONSTANTS[116] =  3.00000*CONSTANTS[41];
CONSTANTS[117] = (CONSTANTS[0]==1.00000 ?  CONSTANTS[42]*0.600000 : CONSTANTS[42]);
CONSTANTS[118] = (CONSTANTS[0]==1.00000 ?  CONSTANTS[43]*2.00000 : CONSTANTS[0]==2.00000 ?  CONSTANTS[43]*2.00000 : CONSTANTS[43]);
CONSTANTS[119] = 1.00000 - CONSTANTS[48];
CONSTANTS[120] = (CONSTANTS[0]==1.00000 ?  CONSTANTS[47]*1.20000 : CONSTANTS[0]==2.00000 ?  CONSTANTS[47]*2.00000 : CONSTANTS[47]);
CONSTANTS[121] = ( 0.500000*(CONSTANTS[1]+CONSTANTS[3]+CONSTANTS[4]+ 4.00000*CONSTANTS[2]))/1000.00;
CONSTANTS[122] = (CONSTANTS[0]==1.00000 ?  CONSTANTS[54]*1.30000 : CONSTANTS[0]==2.00000 ?  CONSTANTS[54]*0.800000 : CONSTANTS[54]);
CONSTANTS[123] = (CONSTANTS[0]==1.00000 ?  CONSTANTS[57]*1.40000 : CONSTANTS[57]);
CONSTANTS[124] = (CONSTANTS[0]==1.00000 ?  CONSTANTS[58]*1.20000 : CONSTANTS[0]==2.00000 ?  CONSTANTS[58]*1.30000 : CONSTANTS[58]);
CONSTANTS[125] = (CONSTANTS[0]==1.00000 ?  CONSTANTS[95]*0.600000 : CONSTANTS[95]);
CONSTANTS[126] = ( 0.500000*CONSTANTS[108])/1.00000;
CONSTANTS[127] =  1.25000*CONSTANTS[108];
CONSTANTS[128] = (CONSTANTS[0]==1.00000 ? 1.30000 : 1.00000);
CONSTANTS[129] =  2.00000*3.14000*CONSTANTS[13]*CONSTANTS[13]+ 2.00000*3.14000*CONSTANTS[13]*CONSTANTS[12];
CONSTANTS[130] =  1.10000*CONSTANTS[120];
CONSTANTS[131] =  0.00125000*CONSTANTS[120];
CONSTANTS[132] =  0.000357400*CONSTANTS[120];
CONSTANTS[133] =  1.82000e+06*pow( CONSTANTS[52]*CONSTANTS[6], - 1.50000);
CONSTANTS[134] = ( 0.500000*CONSTANTS[127])/1.00000;
CONSTANTS[135] =  2.00000*CONSTANTS[129];
CONSTANTS[136] =  0.00125000*CONSTANTS[130];
CONSTANTS[137] =  0.000357400*CONSTANTS[130];
CONSTANTS[138] = exp( - CONSTANTS[133]*4.00000*( pow(CONSTANTS[121], 1.0 / 2)/(1.00000+ pow(CONSTANTS[121], 1.0 / 2)) -  0.300000*CONSTANTS[121]));
CONSTANTS[139] = exp( - CONSTANTS[133]*1.00000*( pow(CONSTANTS[121], 1.0 / 2)/(1.00000+ pow(CONSTANTS[121], 1.0 / 2)) -  0.300000*CONSTANTS[121]));
CONSTANTS[140] = exp( - CONSTANTS[133]*1.00000*( pow(CONSTANTS[121], 1.0 / 2)/(1.00000+ pow(CONSTANTS[121], 1.0 / 2)) -  0.300000*CONSTANTS[121]));
CONSTANTS[141] =  0.680000*CONSTANTS[112];
CONSTANTS[142] =  0.0552000*CONSTANTS[112];
CONSTANTS[143] =  0.00480000*CONSTANTS[112];
CONSTANTS[144] =  0.0200000*CONSTANTS[112];
CONSTANTS[145] = CONSTANTS[63]+1.00000+ (CONSTANTS[1]/CONSTANTS[60])*(1.00000+CONSTANTS[1]/CONSTANTS[61]);
CONSTANTS[146] = ( CONSTANTS[1]*CONSTANTS[1])/( CONSTANTS[145]*CONSTANTS[60]*CONSTANTS[61]);
CONSTANTS[147] = 1.00000/CONSTANTS[145];
CONSTANTS[148] =  CONSTANTS[147]*CONSTANTS[2]*CONSTANTS[67];
CONSTANTS[149] = CONSTANTS[68];
CONSTANTS[150] = CONSTANTS[68];
CONSTANTS[151] = (CONSTANTS[0]==1.00000 ?  CONSTANTS[72]*1.10000 : CONSTANTS[0]==2.00000 ?  CONSTANTS[72]*1.40000 : CONSTANTS[72]);
CONSTANTS[152] = CONSTANTS[63]+1.00000+ (CONSTANTS[1]/CONSTANTS[60])*(1.00000+CONSTANTS[1]/CONSTANTS[61]);
CONSTANTS[153] = ( CONSTANTS[1]*CONSTANTS[1])/( CONSTANTS[152]*CONSTANTS[60]*CONSTANTS[61]);
CONSTANTS[154] = 1.00000/CONSTANTS[152];
CONSTANTS[155] =  CONSTANTS[154]*CONSTANTS[2]*CONSTANTS[67];
CONSTANTS[156] = CONSTANTS[68];
CONSTANTS[157] = CONSTANTS[68];
CONSTANTS[158] =  CONSTANTS[74]*CONSTANTS[86];
CONSTANTS[159] = CONSTANTS[75];
CONSTANTS[160] = (( CONSTANTS[79]*CONSTANTS[87])/CONSTANTS[88])/(1.00000+CONSTANTS[87]/CONSTANTS[88]);
CONSTANTS[161] = (CONSTANTS[0]==1.00000 ?  CONSTANTS[94]*0.900000 : CONSTANTS[0]==2.00000 ?  CONSTANTS[94]*0.700000 : CONSTANTS[94]);

STATES[0] = -8.974808e+01;
STATES[1] = 1.095026e-02;
STATES[2] = 6.497341e-05;
STATES[3] = 1.239736e+01;
STATES[4] = 1.239770e+01;
STATES[5] = 1.477115e+02;
STATES[6] = 1.477114e+02;
STATES[7] = 1.528001e+00;
STATES[8] = 1.525693e+00;
STATES[9] = 7.453481e-05;
STATES[10] = 2.920698e+01;
STATES[11] = 2.920696e+01;
STATES[12] = 6.517154e-04;
STATES[13] = 8.473267e-01;
STATES[14] = 8.471657e-01;
STATES[15] = 7.018454e-01;
STATES[16] = 8.469014e-01;
STATES[17] = 1.351203e-04;
STATES[18] = 5.566017e-01;
STATES[19] = 3.115491e-01;
STATES[20] = 8.899259e-04;
STATES[21] = 9.996716e-01;
STATES[22] = 5.988908e-01;
STATES[23] = 4.534165e-04;
STATES[24] = 9.996716e-01;
STATES[25] = 6.620692e-01;
STATES[26] = 1.588841e-31;
STATES[27] = 1.000000e+00;
STATES[28] = 9.401791e-01;
STATES[29] = 1.000000e+00;
STATES[30] = 9.999014e-01;
STATES[31] = 9.999846e-01;
STATES[32] = 1.000000e+00;
STATES[33] = 1.000000e+00;
STATES[34] = 4.899378e-04;
STATES[35] = 8.326009e-04;
STATES[36] = 9.982511e-01;
STATES[37] = 7.936020e-04;
STATES[38] = 6.532143e-04;
STATES[39] = 9.804083e-06;
STATES[40] = 2.922449e-04;
STATES[41] = 2.439590e-01;
STATES[42] = 1.586167e-04;
STATES[43] = 1.808248e-22;
STATES[44] = 4.358608e-21;
}
void
computeRates(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC)
{
RATES[0] = - (ALGEBRAIC[70]+ALGEBRAIC[72]+ALGEBRAIC[78]+ALGEBRAIC[114]+ALGEBRAIC[115]+ALGEBRAIC[116]+ALGEBRAIC[117]+ALGEBRAIC[119]+ALGEBRAIC[123]+ALGEBRAIC[155]+ALGEBRAIC[185]+ALGEBRAIC[204]+ALGEBRAIC[207]+ALGEBRAIC[206]+ALGEBRAIC[211]+ALGEBRAIC[209]+ALGEBRAIC[216]+ALGEBRAIC[218]+ALGEBRAIC[68]+ALGEBRAIC[11]);
RATES[1] =  CONSTANTS[20]*ALGEBRAIC[43]*(ALGEBRAIC[43]+STATES[1]) -  CONSTANTS[21]*STATES[1];
RATES[2] =  ALGEBRAIC[56]*((( - (ALGEBRAIC[93] -  2.00000*ALGEBRAIC[185])*CONSTANTS[135])/( 2.00000*CONSTANTS[7]*CONSTANTS[144])+( ALGEBRAIC[217]*CONSTANTS[143])/CONSTANTS[144]) - ALGEBRAIC[213]);
RATES[3] = ( - (ALGEBRAIC[70]+ALGEBRAIC[72]+ 3.00000*ALGEBRAIC[155]+ALGEBRAIC[112]+ 3.00000*ALGEBRAIC[204]+ALGEBRAIC[207])*CONSTANTS[135])/( CONSTANTS[7]*CONSTANTS[141])+( ALGEBRAIC[210]*CONSTANTS[144])/CONSTANTS[141];
RATES[4] = ( - (ALGEBRAIC[94]+ 3.00000*ALGEBRAIC[185])*CONSTANTS[135])/( CONSTANTS[7]*CONSTANTS[144]) - ALGEBRAIC[210];
RATES[5] = ( - (((ALGEBRAIC[78]+ALGEBRAIC[117]+ALGEBRAIC[119]+ALGEBRAIC[123]+ALGEBRAIC[206]+ALGEBRAIC[68]+ALGEBRAIC[11]) -  2.00000*ALGEBRAIC[204])+ALGEBRAIC[113])*CONSTANTS[135])/( CONSTANTS[7]*CONSTANTS[141])+( ALGEBRAIC[208]*CONSTANTS[144])/CONSTANTS[141];
RATES[6] = ( - ALGEBRAIC[97]*CONSTANTS[135])/( CONSTANTS[7]*CONSTANTS[144]) - ALGEBRAIC[208];
RATES[7] = ALGEBRAIC[224] - ( ALGEBRAIC[225]*CONSTANTS[143])/CONSTANTS[142];
RATES[8] =  ALGEBRAIC[59]*(ALGEBRAIC[225] - ALGEBRAIC[217]);
RATES[9] =  ALGEBRAIC[53]*((( - ((ALGEBRAIC[111]+ALGEBRAIC[211]+ALGEBRAIC[209]) -  2.00000*ALGEBRAIC[155])*CONSTANTS[135])/( 2.00000*CONSTANTS[7]*CONSTANTS[141]) - ( ALGEBRAIC[224]*CONSTANTS[142])/CONSTANTS[141])+( ALGEBRAIC[213]*CONSTANTS[144])/CONSTANTS[141]);
RATES[10] = ( (ALGEBRAIC[218]+ALGEBRAIC[214])*CONSTANTS[135])/( CONSTANTS[7]*CONSTANTS[141])+( ALGEBRAIC[221]*CONSTANTS[144])/CONSTANTS[141];
RATES[11] = ( ALGEBRAIC[212]*CONSTANTS[135])/( CONSTANTS[7]*CONSTANTS[144]) - ALGEBRAIC[221];
RATES[12] = (ALGEBRAIC[0] - STATES[12])/ALGEBRAIC[13];
RATES[13] = (ALGEBRAIC[1] - STATES[13])/ALGEBRAIC[37];
RATES[14] = (ALGEBRAIC[38] - STATES[14])/ALGEBRAIC[44];
RATES[15] = (ALGEBRAIC[45] - STATES[15])/ALGEBRAIC[37];
RATES[16] = (ALGEBRAIC[38] - STATES[16])/ALGEBRAIC[50];
RATES[17] = (ALGEBRAIC[2] - STATES[17])/ALGEBRAIC[16];
RATES[18] = (ALGEBRAIC[3] - STATES[18])/CONSTANTS[41];
RATES[19] = (ALGEBRAIC[4] - STATES[19])/CONSTANTS[116];
RATES[20] = (ALGEBRAIC[5] - STATES[20])/ALGEBRAIC[17];
RATES[21] = (ALGEBRAIC[6] - STATES[21])/ALGEBRAIC[46];
RATES[22] = (ALGEBRAIC[6] - STATES[22])/ALGEBRAIC[51];
RATES[23] = (ALGEBRAIC[32] - STATES[23])/ALGEBRAIC[17];
RATES[24] = (ALGEBRAIC[6] - STATES[24])/ALGEBRAIC[60];
RATES[25] = (ALGEBRAIC[6] - STATES[25])/ALGEBRAIC[61];
RATES[26] = (ALGEBRAIC[7] - STATES[26])/ALGEBRAIC[22];
RATES[27] = (ALGEBRAIC[8] - STATES[27])/ALGEBRAIC[23];
RATES[28] = (ALGEBRAIC[8] - STATES[28])/ALGEBRAIC[24];
RATES[29] = (ALGEBRAIC[19] - STATES[29])/ALGEBRAIC[33];
RATES[30] = (ALGEBRAIC[19] - STATES[30])/ALGEBRAIC[34];
RATES[31] = (ALGEBRAIC[9] - STATES[31])/CONSTANTS[49];
RATES[32] = (ALGEBRAIC[8] - STATES[32])/ALGEBRAIC[35];
RATES[33] = (ALGEBRAIC[19] - STATES[33])/ALGEBRAIC[40];
RATES[34] =  ALGEBRAIC[20]*CONSTANTS[46] -  STATES[34]*ALGEBRAIC[10];
RATES[35] =  ALGEBRAIC[21]*CONSTANTS[46] -  STATES[35]*ALGEBRAIC[10];
RATES[36] = ( CONSTANTS[55]*STATES[37]+ ALGEBRAIC[48]*STATES[40]+ ALGEBRAIC[62]*STATES[39]) -  (CONSTANTS[56]+ALGEBRAIC[42]+ALGEBRAIC[58])*STATES[36];
RATES[37] = ( ALGEBRAIC[41]*STATES[38]+ CONSTANTS[56]*STATES[36]) -  (ALGEBRAIC[47]+CONSTANTS[55])*STATES[37];
RATES[38] =  ALGEBRAIC[47]*STATES[37] -  ALGEBRAIC[41]*STATES[38];
RATES[39] = ( ALGEBRAIC[58]*STATES[36]+ ALGEBRAIC[52]*STATES[40]) -  (ALGEBRAIC[62]+ALGEBRAIC[55])*STATES[39];
RATES[40] = ( ALGEBRAIC[42]*STATES[36]+ ALGEBRAIC[55]*STATES[39]) -  (ALGEBRAIC[48]+ALGEBRAIC[52])*STATES[40];
RATES[41] = (ALGEBRAIC[12] - STATES[41])/ALGEBRAIC[27];
RATES[42] = (ALGEBRAIC[26] - STATES[42])/ALGEBRAIC[36];
RATES[43] = (ALGEBRAIC[98] - STATES[43])/ALGEBRAIC[104];
RATES[44] = (ALGEBRAIC[99] - STATES[44])/ALGEBRAIC[105];

ALGEBRAIC[0] = 1.00000/pow(1.00000+exp(- (STATES[0]+56.8600)/9.03000), 2.00000);
ALGEBRAIC[1] = 1.00000/pow(1.00000+exp((STATES[0]+71.5500)/7.43000), 2.00000);
ALGEBRAIC[2] = 1.00000/(1.00000+exp(- (STATES[0]+42.8500)/5.26400));
ALGEBRAIC[3] = 1.00000/(1.00000+exp((STATES[0]+87.6100)/7.48800));
ALGEBRAIC[4] = 1.00000/(1.00000+exp((STATES[0]+93.8100)/7.48800));
ALGEBRAIC[5] = 1.00000/(1.00000+exp(- ((STATES[0]+CONSTANTS[44]) - 14.3400)/14.8200));
ALGEBRAIC[6] = 1.00000/(1.00000+exp((STATES[0]+CONSTANTS[44]+43.9400)/5.71100));
ALGEBRAIC[7] = (STATES[0]>=31.4978 ? 1.00000 :  1.07630*exp( - 1.00700*exp( - 0.0829000*STATES[0])));
ALGEBRAIC[8] = 1.00000/(1.00000+exp((STATES[0]+19.5800)/3.69600));
ALGEBRAIC[9] = 1.00000/(1.00000+exp((STATES[0]+18.0800)/2.79160));
ALGEBRAIC[10] =  STATES[31]*1.00000;
ALGEBRAIC[11] = (VOI>=CONSTANTS[14]&&(VOI - CONSTANTS[14]) -  floor((VOI - CONSTANTS[14])/CONSTANTS[17])*CONSTANTS[17]<=CONSTANTS[18] ? CONSTANTS[16] : 0.00000);
ALGEBRAIC[12] = 1.00000/(1.00000+exp(- (STATES[0]+11.6000)/8.93200));
ALGEBRAIC[13] =  0.129200*exp(- pow((STATES[0]+45.7900)/15.5400, 2.00000))+ 0.0648700*exp(- pow((STATES[0] - 4.82300)/51.1200, 2.00000));
ALGEBRAIC[14] = (STATES[0]>=- 40.0000 ? 0.00000 :  0.0570000*exp(- (STATES[0]+80.0000)/6.80000));
ALGEBRAIC[15] = (STATES[0]>=- 40.0000 ? 0.00000 : ( ( - 25428.0*exp( 0.244400*STATES[0]) -  6.94800e-06*exp( - 0.0439100*STATES[0]))*(STATES[0]+37.7800))/(1.00000+exp( 0.311000*(STATES[0]+79.2300))));
ALGEBRAIC[16] =  0.129200*exp(- pow((STATES[0]+45.7900)/15.5400, 2.00000))+ 0.0648700*exp(- pow((STATES[0] - 4.82300)/51.1200, 2.00000));
ALGEBRAIC[17] = 1.05150/(1.00000/( 1.20890*(1.00000+exp(- ((STATES[0]+CONSTANTS[44]) - 18.4099)/29.3814)))+3.50000/(1.00000+exp((STATES[0]+CONSTANTS[44]+100.000)/29.3814)));
ALGEBRAIC[18] = (CONSTANTS[0]==1.00000 ? 1.00000 - 0.950000/(1.00000+exp((STATES[0]+CONSTANTS[44]+70.0000)/5.00000)) : 1.00000);
ALGEBRAIC[19] = ALGEBRAIC[8];
ALGEBRAIC[20] = 1.00000/(CONSTANTS[46]/ALGEBRAIC[10]+pow(1.00000+CONSTANTS[45]/STATES[2], 4.00000));
ALGEBRAIC[21] = 1.00000/(CONSTANTS[46]/ALGEBRAIC[10]+pow(1.00000+CONSTANTS[45]/STATES[9], 4.00000));
ALGEBRAIC[22] = CONSTANTS[51]+0.600000+1.00000/(exp( - 0.0500000*(STATES[0]+CONSTANTS[50]+6.00000))+exp( 0.0900000*(STATES[0]+CONSTANTS[50]+14.0000)));
ALGEBRAIC[23] = 7.00000+1.00000/( 0.00450000*exp(- (STATES[0]+20.0000)/10.0000)+ 0.00450000*exp((STATES[0]+20.0000)/10.0000));
ALGEBRAIC[24] = 1000.00+1.00000/( 3.50000e-05*exp(- (STATES[0]+5.00000)/4.00000)+ 3.50000e-05*exp((STATES[0]+5.00000)/6.00000));
ALGEBRAIC[25] = ( STATES[0]*CONSTANTS[7]*CONSTANTS[7])/( CONSTANTS[5]*CONSTANTS[6]);
ALGEBRAIC[26] = ALGEBRAIC[12];
ALGEBRAIC[27] = 817.300+1.00000/( 0.000232600*exp((STATES[0]+48.2800)/17.8000)+ 0.00129200*exp(- (STATES[0]+210.000)/230.000));
ALGEBRAIC[28] = ( STATES[0]*CONSTANTS[7])/( CONSTANTS[5]*CONSTANTS[6]);
ALGEBRAIC[29] = (STATES[0]>=- 40.0000 ? 0.770000/( 0.130000*(1.00000+exp(- (STATES[0]+10.6600)/11.1000))) :  2.70000*exp( 0.0790000*STATES[0])+ 310000.*exp( 0.348500*STATES[0]));
ALGEBRAIC[30] = (STATES[0]>=- 40.0000 ? ( 0.600000*exp( 0.0570000*STATES[0]))/(1.00000+exp( - 0.100000*(STATES[0]+32.0000))) : ( 0.0242400*exp( - 0.0105200*STATES[0]))/(1.00000+exp( - 0.137800*(STATES[0]+40.1400))));
ALGEBRAIC[31] = 4.56200+1.00000/( 0.393300*exp(- (STATES[0]+CONSTANTS[44]+100.000)/100.000)+ 0.0800400*exp((STATES[0]+CONSTANTS[44]+50.0000)/16.5900));
ALGEBRAIC[32] = 1.00000/(1.00000+exp(- ((STATES[0]+CONSTANTS[44]) - 24.3400)/14.8200));
ALGEBRAIC[33] = 7.00000+1.00000/( 0.0400000*exp(- (STATES[0] - 4.00000)/7.00000)+ 0.0400000*exp((STATES[0] - 4.00000)/7.00000));
ALGEBRAIC[34] = 100.000+1.00000/( 0.000120000*exp(- STATES[0]/3.00000)+ 0.000120000*exp(STATES[0]/7.00000));
ALGEBRAIC[35] =  2.50000*ALGEBRAIC[23];
ALGEBRAIC[36] = 1.00000/( 0.0100000*exp((STATES[0] - 50.0000)/20.0000)+ 0.0193000*exp(- (STATES[0]+66.5400)/31.0000));
ALGEBRAIC[37] = 1.00000/(ALGEBRAIC[14]+ALGEBRAIC[29]);
ALGEBRAIC[38] = ALGEBRAIC[1];
ALGEBRAIC[39] = 23.6200+1.00000/( 0.00141600*exp(- (STATES[0]+CONSTANTS[44]+96.5200)/59.0500)+ 1.78000e-08*exp((STATES[0]+CONSTANTS[44]+114.100)/8.07900));
ALGEBRAIC[40] =  2.50000*ALGEBRAIC[33];
ALGEBRAIC[41] =  0.116100*exp( 0.299000*ALGEBRAIC[28]);
ALGEBRAIC[42] =  0.0578000*exp( 0.971000*ALGEBRAIC[28]);
ALGEBRAIC[43] = ( CONSTANTS[22]*(1.00000 - STATES[1]))/(1.00000+CONSTANTS[23]/STATES[2]);
ALGEBRAIC[44] = 1.00000/(ALGEBRAIC[15]+ALGEBRAIC[30]);
ALGEBRAIC[45] = 1.00000/pow(1.00000+exp((STATES[0]+77.5500)/7.43000), 2.00000);
ALGEBRAIC[46] =  ALGEBRAIC[31]*ALGEBRAIC[18];
ALGEBRAIC[47] =  0.244200*exp( - 1.60400*ALGEBRAIC[28]);
ALGEBRAIC[48] =  0.000349000*exp( - 1.06200*ALGEBRAIC[28]);
ALGEBRAIC[49] = ALGEBRAIC[43]+STATES[1];
ALGEBRAIC[50] =  1.46000*ALGEBRAIC[44];
ALGEBRAIC[51] =  ALGEBRAIC[39]*ALGEBRAIC[18];
ALGEBRAIC[52] =  0.253300*exp( 0.595300*ALGEBRAIC[28]);
ALGEBRAIC[53] = 1.00000/(1.00000+( CONSTANTS[113]*CONSTANTS[25])/pow(CONSTANTS[25]+STATES[9], 2.00000)+( CONSTANTS[26]*CONSTANTS[27])/pow(CONSTANTS[27]+STATES[9], 2.00000));
ALGEBRAIC[54] = 1.35400+0.000100000/(exp(((STATES[0]+CONSTANTS[44]) - 167.400)/15.8900)+exp(- ((STATES[0]+CONSTANTS[44]) - 12.2300)/0.215400));
ALGEBRAIC[55] =  0.0652500*exp( - 0.820900*ALGEBRAIC[28]);
ALGEBRAIC[56] = 1.00000/(1.00000+( CONSTANTS[28]*CONSTANTS[29])/pow(CONSTANTS[29]+STATES[2], 2.00000)+( CONSTANTS[30]*CONSTANTS[31])/pow(CONSTANTS[31]+STATES[2], 2.00000));
ALGEBRAIC[57] = 1.00000 - 0.500000/(1.00000+exp((STATES[0]+CONSTANTS[44]+70.0000)/20.0000));
ALGEBRAIC[58] =  5.20000e-05*exp( 1.52500*ALGEBRAIC[28]);
ALGEBRAIC[59] = 1.00000/(1.00000+( CONSTANTS[32]*CONSTANTS[33])/pow(CONSTANTS[33]+STATES[8], 2.00000));
ALGEBRAIC[60] =  ALGEBRAIC[54]*ALGEBRAIC[57]*ALGEBRAIC[46];
ALGEBRAIC[61] =  ALGEBRAIC[54]*ALGEBRAIC[57]*ALGEBRAIC[51];
ALGEBRAIC[62] = ( ALGEBRAIC[48]*ALGEBRAIC[55]*ALGEBRAIC[58])/( ALGEBRAIC[42]*ALGEBRAIC[52]);
ALGEBRAIC[63] =  (( CONSTANTS[5]*CONSTANTS[6])/( CONSTANTS[8]*CONSTANTS[7]))*log(CONSTANTS[1]/STATES[3]);
ALGEBRAIC[64] =  (( CONSTANTS[5]*CONSTANTS[6])/( CONSTANTS[10]*CONSTANTS[7]))*log(CONSTANTS[3]/STATES[5]);
ALGEBRAIC[65] =  (( CONSTANTS[5]*CONSTANTS[6])/( CONSTANTS[10]*CONSTANTS[7]))*log((CONSTANTS[3]+ CONSTANTS[34]*CONSTANTS[1])/(STATES[5]+ CONSTANTS[34]*STATES[3]));
ALGEBRAIC[66] =  (( CONSTANTS[5]*CONSTANTS[6])/( CONSTANTS[11]*CONSTANTS[7]))*log(CONSTANTS[4]/STATES[10]);
ALGEBRAIC[67] =  (( CONSTANTS[5]*CONSTANTS[6])/( CONSTANTS[11]*CONSTANTS[7]))*log(CONSTANTS[4]/STATES[11]);
ALGEBRAIC[68] =  CONSTANTS[36]*CONSTANTS[35]*CONSTANTS[114]*CONSTANTS[115]*(STATES[0] - ALGEBRAIC[64]);
ALGEBRAIC[69] = 1.00000/(1.00000+CONSTANTS[19]/ALGEBRAIC[49]);
ALGEBRAIC[70] =  CONSTANTS[40]*(STATES[0] - ALGEBRAIC[63])*pow(STATES[12], 3.00000)*( (1.00000 - ALGEBRAIC[69])*STATES[13]*STATES[14]+ ALGEBRAIC[69]*STATES[15]*STATES[16]);
ALGEBRAIC[71] = 1.00000/(1.00000+CONSTANTS[19]/ALGEBRAIC[49]);
ALGEBRAIC[72] =  CONSTANTS[117]*(STATES[0] - ALGEBRAIC[63])*STATES[17]*( (1.00000 - ALGEBRAIC[71])*STATES[18]+ ALGEBRAIC[71]*STATES[19]);
ALGEBRAIC[73] = 1.00000/(1.00000+exp(((STATES[0]+CONSTANTS[44]) - 213.600)/151.200));
ALGEBRAIC[74] = 1.00000 - ALGEBRAIC[73];
ALGEBRAIC[75] =  ALGEBRAIC[73]*STATES[21]+ ALGEBRAIC[74]*STATES[22];
ALGEBRAIC[76] =  ALGEBRAIC[73]*STATES[24]+ ALGEBRAIC[74]*STATES[25];
ALGEBRAIC[77] = 1.00000/(1.00000+CONSTANTS[19]/ALGEBRAIC[49]);
ALGEBRAIC[78] =  CONSTANTS[118]*(STATES[0] - ALGEBRAIC[64])*( (1.00000 - ALGEBRAIC[77])*STATES[20]*ALGEBRAIC[75]+ ALGEBRAIC[77]*STATES[23]*ALGEBRAIC[76]);
ALGEBRAIC[79] =  CONSTANTS[48]*STATES[27]+ CONSTANTS[119]*STATES[28];
ALGEBRAIC[80] = 0.300000+0.600000/(1.00000+exp((STATES[0] - 10.0000)/10.0000));
ALGEBRAIC[81] = 1.00000 - ALGEBRAIC[80];
ALGEBRAIC[82] =  ALGEBRAIC[80]*STATES[29]+ ALGEBRAIC[81]*STATES[30];
ALGEBRAIC[83] =  CONSTANTS[48]*STATES[32]+ CONSTANTS[119]*STATES[28];
ALGEBRAIC[84] =  ALGEBRAIC[80]*STATES[33]+ ALGEBRAIC[81]*STATES[30];
ALGEBRAIC[85] = ( 0.500000*(STATES[4]+STATES[6]+STATES[11]+ 4.00000*STATES[2]))/1000.00;
ALGEBRAIC[86] = exp( - CONSTANTS[133]*4.00000*( pow(ALGEBRAIC[85], 1.0 / 2)/(1.00000+ pow(ALGEBRAIC[85], 1.0 / 2)) -  0.300000*ALGEBRAIC[85]));
ALGEBRAIC[87] = exp( - CONSTANTS[133]*1.00000*( pow(ALGEBRAIC[85], 1.0 / 2)/(1.00000+ pow(ALGEBRAIC[85], 1.0 / 2)) -  0.300000*ALGEBRAIC[85]));
ALGEBRAIC[88] = exp( - CONSTANTS[133]*1.00000*( pow(ALGEBRAIC[85], 1.0 / 2)/(1.00000+ pow(ALGEBRAIC[85], 1.0 / 2)) -  0.300000*ALGEBRAIC[85]));
ALGEBRAIC[89] = ( 4.00000*ALGEBRAIC[25]*( ALGEBRAIC[86]*STATES[2]*exp( 2.00000*ALGEBRAIC[28]) -  CONSTANTS[138]*CONSTANTS[2]))/(exp( 2.00000*ALGEBRAIC[28]) - 1.00000);
ALGEBRAIC[90] = ( 1.00000*ALGEBRAIC[25]*( ALGEBRAIC[87]*STATES[4]*exp( 1.00000*ALGEBRAIC[28]) -  CONSTANTS[139]*CONSTANTS[1]))/(exp( 1.00000*ALGEBRAIC[28]) - 1.00000);
ALGEBRAIC[91] = ( 1.00000*ALGEBRAIC[25]*( ALGEBRAIC[88]*STATES[6]*exp( 1.00000*ALGEBRAIC[28]) -  CONSTANTS[140]*CONSTANTS[3]))/(exp( 1.00000*ALGEBRAIC[28]) - 1.00000);
ALGEBRAIC[92] = 1.00000/(1.00000+CONSTANTS[19]/ALGEBRAIC[49]);
ALGEBRAIC[93] =  CONSTANTS[53]*( (1.00000 - ALGEBRAIC[92])*CONSTANTS[120]*ALGEBRAIC[89]*STATES[26]*( ALGEBRAIC[79]*(1.00000 - STATES[34])+ STATES[31]*ALGEBRAIC[82]*STATES[34])+ ALGEBRAIC[92]*CONSTANTS[130]*ALGEBRAIC[89]*STATES[26]*( ALGEBRAIC[83]*(1.00000 - STATES[34])+ STATES[31]*ALGEBRAIC[84]*STATES[34]));
ALGEBRAIC[94] =  CONSTANTS[53]*( (1.00000 - ALGEBRAIC[92])*CONSTANTS[131]*ALGEBRAIC[90]*STATES[26]*( ALGEBRAIC[79]*(1.00000 - STATES[34])+ STATES[31]*ALGEBRAIC[82]*STATES[34])+ ALGEBRAIC[92]*CONSTANTS[136]*ALGEBRAIC[90]*STATES[26]*( ALGEBRAIC[83]*(1.00000 - STATES[34])+ STATES[31]*ALGEBRAIC[84]*STATES[34]));
ALGEBRAIC[95] = (( - CONSTANTS[126]*ALGEBRAIC[93])/1.00000)/(1.00000+pow(CONSTANTS[109]/STATES[8], 8.00000));
ALGEBRAIC[96] = (( - CONSTANTS[134]*ALGEBRAIC[93])/1.00000)/(1.00000+pow(CONSTANTS[109]/STATES[8], 8.00000));
ALGEBRAIC[97] =  CONSTANTS[53]*( (1.00000 - ALGEBRAIC[92])*CONSTANTS[132]*ALGEBRAIC[91]*STATES[26]*( ALGEBRAIC[79]*(1.00000 - STATES[34])+ STATES[31]*ALGEBRAIC[82]*STATES[34])+ ALGEBRAIC[92]*CONSTANTS[137]*ALGEBRAIC[91]*STATES[26]*( ALGEBRAIC[83]*(1.00000 - STATES[34])+ STATES[31]*ALGEBRAIC[84]*STATES[34]));
ALGEBRAIC[98] = (CONSTANTS[0]==2.00000 ?  ALGEBRAIC[95]*1.70000 : ALGEBRAIC[95]);
ALGEBRAIC[99] = (CONSTANTS[0]==2.00000 ?  ALGEBRAIC[96]*1.70000 : ALGEBRAIC[96]);
ALGEBRAIC[100] = ( 0.500000*(STATES[3]+STATES[5]+STATES[10]+ 4.00000*STATES[9]))/1000.00;
ALGEBRAIC[101] = CONSTANTS[108]/(1.00000+0.0123000/STATES[8]);
ALGEBRAIC[102] = CONSTANTS[127]/(1.00000+0.0123000/STATES[8]);
ALGEBRAIC[103] = exp( - CONSTANTS[133]*4.00000*( pow(ALGEBRAIC[100], 1.0 / 2)/(1.00000+ pow(ALGEBRAIC[100], 1.0 / 2)) -  0.300000*ALGEBRAIC[100]));
ALGEBRAIC[104] = (ALGEBRAIC[101]<0.00100000 ? 0.00100000 : ALGEBRAIC[101]);
ALGEBRAIC[105] = (ALGEBRAIC[102]<0.00100000 ? 0.00100000 : ALGEBRAIC[102]);
ALGEBRAIC[106] = exp( - CONSTANTS[133]*1.00000*( pow(ALGEBRAIC[100], 1.0 / 2)/(1.00000+ pow(ALGEBRAIC[100], 1.0 / 2)) -  0.300000*ALGEBRAIC[100]));
ALGEBRAIC[107] = exp( - CONSTANTS[133]*1.00000*( pow(ALGEBRAIC[100], 1.0 / 2)/(1.00000+ pow(ALGEBRAIC[100], 1.0 / 2)) -  0.300000*ALGEBRAIC[100]));
ALGEBRAIC[108] = ( 4.00000*ALGEBRAIC[25]*( ALGEBRAIC[103]*STATES[9]*exp( 2.00000*ALGEBRAIC[28]) -  CONSTANTS[138]*CONSTANTS[2]))/(exp( 2.00000*ALGEBRAIC[28]) - 1.00000);
ALGEBRAIC[109] = ( 1.00000*ALGEBRAIC[25]*( ALGEBRAIC[106]*STATES[3]*exp( 1.00000*ALGEBRAIC[28]) -  CONSTANTS[139]*CONSTANTS[1]))/(exp( 1.00000*ALGEBRAIC[28]) - 1.00000);
ALGEBRAIC[110] = ( 1.00000*ALGEBRAIC[25]*( ALGEBRAIC[107]*STATES[5]*exp( 1.00000*ALGEBRAIC[28]) -  CONSTANTS[140]*CONSTANTS[3]))/(exp( 1.00000*ALGEBRAIC[28]) - 1.00000);
ALGEBRAIC[111] =  (1.00000 - CONSTANTS[53])*( (1.00000 - ALGEBRAIC[92])*CONSTANTS[120]*ALGEBRAIC[108]*STATES[26]*( ALGEBRAIC[79]*(1.00000 - STATES[35])+ STATES[31]*ALGEBRAIC[82]*STATES[35])+ ALGEBRAIC[92]*CONSTANTS[130]*ALGEBRAIC[108]*STATES[26]*( ALGEBRAIC[83]*(1.00000 - STATES[35])+ STATES[31]*ALGEBRAIC[84]*STATES[35]));
ALGEBRAIC[112] =  (1.00000 - CONSTANTS[53])*( (1.00000 - ALGEBRAIC[92])*CONSTANTS[131]*ALGEBRAIC[109]*STATES[26]*( ALGEBRAIC[79]*(1.00000 - STATES[35])+ STATES[31]*ALGEBRAIC[82]*STATES[35])+ ALGEBRAIC[92]*CONSTANTS[136]*ALGEBRAIC[109]*STATES[26]*( ALGEBRAIC[83]*(1.00000 - STATES[35])+ STATES[31]*ALGEBRAIC[84]*STATES[35]));
ALGEBRAIC[113] =  (1.00000 - CONSTANTS[53])*( (1.00000 - ALGEBRAIC[92])*CONSTANTS[132]*ALGEBRAIC[110]*STATES[26]*( ALGEBRAIC[79]*(1.00000 - STATES[35])+ STATES[31]*ALGEBRAIC[82]*STATES[35])+ ALGEBRAIC[92]*CONSTANTS[137]*ALGEBRAIC[110]*STATES[26]*( ALGEBRAIC[83]*(1.00000 - STATES[35])+ STATES[31]*ALGEBRAIC[84]*STATES[35]));
ALGEBRAIC[114] = ALGEBRAIC[93]+ALGEBRAIC[111];
ALGEBRAIC[115] = ALGEBRAIC[94]+ALGEBRAIC[112];
ALGEBRAIC[116] = ALGEBRAIC[97]+ALGEBRAIC[113];
ALGEBRAIC[117] =  CONSTANTS[122]* pow((CONSTANTS[3]/5.00000), 1.0 / 2)*STATES[40]*(STATES[0] - ALGEBRAIC[64]);
ALGEBRAIC[118] = 1.00000+0.600000/(1.00000+pow(3.80000e-05/STATES[9], 1.40000));
ALGEBRAIC[119] =  CONSTANTS[123]*ALGEBRAIC[118]*STATES[41]*STATES[42]*(STATES[0] - ALGEBRAIC[65]);
ALGEBRAIC[120] = 4.09400/(1.00000+exp( 0.121700*((STATES[0] - ALGEBRAIC[64]) - 49.9340)));
ALGEBRAIC[121] = ( 15.7200*exp( 0.0674000*((STATES[0] - ALGEBRAIC[64]) - 3.25700))+exp( 0.0618000*((STATES[0] - ALGEBRAIC[64]) - 594.310)))/(1.00000+exp( - 0.162900*((STATES[0] - ALGEBRAIC[64])+14.2070)));
ALGEBRAIC[122] = ALGEBRAIC[120]/(ALGEBRAIC[120]+ALGEBRAIC[121]);
ALGEBRAIC[123] =  CONSTANTS[124]* pow((CONSTANTS[3]/5.00000), 1.0 / 2)*ALGEBRAIC[122]*(STATES[0] - ALGEBRAIC[64]);
ALGEBRAIC[124] = exp( CONSTANTS[70]*ALGEBRAIC[28]);
ALGEBRAIC[125] = exp( CONSTANTS[69]*ALGEBRAIC[28]);
ALGEBRAIC[126] = 1.00000+ (STATES[3]/CONSTANTS[62])*(1.00000+ALGEBRAIC[125]);
ALGEBRAIC[127] = ( STATES[3]*ALGEBRAIC[125])/( CONSTANTS[62]*ALGEBRAIC[126]);
ALGEBRAIC[128] = 1.00000/ALGEBRAIC[126];
ALGEBRAIC[129] = 1.00000+ (STATES[3]/CONSTANTS[60])*(1.00000+STATES[3]/CONSTANTS[61]);
ALGEBRAIC[130] = ( STATES[3]*STATES[3])/( ALGEBRAIC[129]*CONSTANTS[60]*CONSTANTS[61]);
ALGEBRAIC[131] = 1.00000/ALGEBRAIC[129];
ALGEBRAIC[132] = 1.00000+ (CONSTANTS[1]/CONSTANTS[62])*(1.00000+1.00000/ALGEBRAIC[125]);
ALGEBRAIC[133] = CONSTANTS[1]/( CONSTANTS[62]*ALGEBRAIC[125]*ALGEBRAIC[132]);
ALGEBRAIC[134] = 1.00000/ALGEBRAIC[132];
ALGEBRAIC[135] =  ALGEBRAIC[134]*CONSTANTS[65];
ALGEBRAIC[136] =  ALGEBRAIC[133]*CONSTANTS[66];
ALGEBRAIC[137] = ALGEBRAIC[135]+ALGEBRAIC[136];
ALGEBRAIC[138] = ( ALGEBRAIC[128]*CONSTANTS[65])/ALGEBRAIC[124];
ALGEBRAIC[139] =  ALGEBRAIC[127]*CONSTANTS[66];
ALGEBRAIC[140] = ALGEBRAIC[138]+ALGEBRAIC[139];
ALGEBRAIC[141] =  ALGEBRAIC[131]*STATES[9]*CONSTANTS[67];
ALGEBRAIC[142] =  ALGEBRAIC[130]*ALGEBRAIC[127]*CONSTANTS[64];
ALGEBRAIC[143] =  ALGEBRAIC[133]*CONSTANTS[146]*CONSTANTS[64];
ALGEBRAIC[144] =  CONSTANTS[149]*ALGEBRAIC[140]*(ALGEBRAIC[142]+ALGEBRAIC[141])+ CONSTANTS[150]*ALGEBRAIC[142]*(CONSTANTS[149]+ALGEBRAIC[137]);
ALGEBRAIC[145] =  CONSTANTS[148]*ALGEBRAIC[142]*(ALGEBRAIC[140]+CONSTANTS[150])+ ALGEBRAIC[140]*ALGEBRAIC[141]*(CONSTANTS[148]+ALGEBRAIC[143]);
ALGEBRAIC[146] =  CONSTANTS[148]*ALGEBRAIC[137]*(ALGEBRAIC[142]+ALGEBRAIC[141])+ ALGEBRAIC[143]*ALGEBRAIC[141]*(CONSTANTS[149]+ALGEBRAIC[137]);
ALGEBRAIC[147] =  CONSTANTS[149]*ALGEBRAIC[143]*(ALGEBRAIC[140]+CONSTANTS[150])+ ALGEBRAIC[137]*CONSTANTS[150]*(CONSTANTS[148]+ALGEBRAIC[143]);
ALGEBRAIC[148] = ALGEBRAIC[144]/(ALGEBRAIC[144]+ALGEBRAIC[145]+ALGEBRAIC[146]+ALGEBRAIC[147]);
ALGEBRAIC[149] = ALGEBRAIC[145]/(ALGEBRAIC[144]+ALGEBRAIC[145]+ALGEBRAIC[146]+ALGEBRAIC[147]);
ALGEBRAIC[150] = ALGEBRAIC[146]/(ALGEBRAIC[144]+ALGEBRAIC[145]+ALGEBRAIC[146]+ALGEBRAIC[147]);
ALGEBRAIC[151] = ALGEBRAIC[147]/(ALGEBRAIC[144]+ALGEBRAIC[145]+ALGEBRAIC[146]+ALGEBRAIC[147]);
ALGEBRAIC[152] = 1.00000/(1.00000+pow(CONSTANTS[71]/STATES[9], 2.00000));
ALGEBRAIC[153] = ( 3.00000*( ALGEBRAIC[151]*ALGEBRAIC[142] -  ALGEBRAIC[148]*ALGEBRAIC[143])+ ALGEBRAIC[150]*ALGEBRAIC[139]) -  ALGEBRAIC[149]*ALGEBRAIC[136];
ALGEBRAIC[154] =  ALGEBRAIC[149]*CONSTANTS[149] -  ALGEBRAIC[148]*CONSTANTS[148];
ALGEBRAIC[155] =  (1.00000 - CONSTANTS[59])*CONSTANTS[151]*ALGEBRAIC[152]*( CONSTANTS[8]*ALGEBRAIC[153]+ CONSTANTS[9]*ALGEBRAIC[154]);
ALGEBRAIC[156] = 1.00000+ (STATES[4]/CONSTANTS[62])*(1.00000+ALGEBRAIC[125]);
ALGEBRAIC[157] = ( STATES[4]*ALGEBRAIC[125])/( CONSTANTS[62]*ALGEBRAIC[156]);
ALGEBRAIC[158] = 1.00000/ALGEBRAIC[156];
ALGEBRAIC[159] = 1.00000+ (STATES[4]/CONSTANTS[60])*(1.00000+STATES[4]/CONSTANTS[61]);
ALGEBRAIC[160] = ( STATES[4]*STATES[4])/( ALGEBRAIC[159]*CONSTANTS[60]*CONSTANTS[61]);
ALGEBRAIC[161] = 1.00000/ALGEBRAIC[159];
ALGEBRAIC[162] = 1.00000+ (CONSTANTS[1]/CONSTANTS[62])*(1.00000+1.00000/ALGEBRAIC[125]);
ALGEBRAIC[163] = CONSTANTS[1]/( CONSTANTS[62]*ALGEBRAIC[125]*ALGEBRAIC[162]);
ALGEBRAIC[164] = 1.00000/ALGEBRAIC[162];
ALGEBRAIC[165] =  ALGEBRAIC[164]*CONSTANTS[65];
ALGEBRAIC[166] =  ALGEBRAIC[163]*CONSTANTS[66];
ALGEBRAIC[167] = ALGEBRAIC[165]+ALGEBRAIC[166];
ALGEBRAIC[168] = ( ALGEBRAIC[158]*CONSTANTS[65])/ALGEBRAIC[124];
ALGEBRAIC[169] =  ALGEBRAIC[157]*CONSTANTS[66];
ALGEBRAIC[170] = ALGEBRAIC[168]+ALGEBRAIC[169];
ALGEBRAIC[171] =  ALGEBRAIC[161]*STATES[2]*CONSTANTS[67];
ALGEBRAIC[172] =  ALGEBRAIC[160]*ALGEBRAIC[157]*CONSTANTS[64];
ALGEBRAIC[173] =  ALGEBRAIC[163]*CONSTANTS[153]*CONSTANTS[64];
ALGEBRAIC[174] =  CONSTANTS[156]*ALGEBRAIC[170]*(ALGEBRAIC[172]+ALGEBRAIC[171])+ CONSTANTS[157]*ALGEBRAIC[172]*(CONSTANTS[156]+ALGEBRAIC[167]);
ALGEBRAIC[175] =  CONSTANTS[155]*ALGEBRAIC[172]*(ALGEBRAIC[170]+CONSTANTS[157])+ ALGEBRAIC[170]*ALGEBRAIC[171]*(CONSTANTS[155]+ALGEBRAIC[173]);
ALGEBRAIC[176] =  CONSTANTS[155]*ALGEBRAIC[167]*(ALGEBRAIC[172]+ALGEBRAIC[171])+ ALGEBRAIC[173]*ALGEBRAIC[171]*(CONSTANTS[156]+ALGEBRAIC[167]);
ALGEBRAIC[177] =  CONSTANTS[156]*ALGEBRAIC[173]*(ALGEBRAIC[170]+CONSTANTS[157])+ ALGEBRAIC[167]*CONSTANTS[157]*(CONSTANTS[155]+ALGEBRAIC[173]);
ALGEBRAIC[178] = ALGEBRAIC[174]/(ALGEBRAIC[174]+ALGEBRAIC[175]+ALGEBRAIC[176]+ALGEBRAIC[177]);
ALGEBRAIC[179] = ALGEBRAIC[175]/(ALGEBRAIC[174]+ALGEBRAIC[175]+ALGEBRAIC[176]+ALGEBRAIC[177]);
ALGEBRAIC[180] = ALGEBRAIC[176]/(ALGEBRAIC[174]+ALGEBRAIC[175]+ALGEBRAIC[176]+ALGEBRAIC[177]);
ALGEBRAIC[181] = ALGEBRAIC[177]/(ALGEBRAIC[174]+ALGEBRAIC[175]+ALGEBRAIC[176]+ALGEBRAIC[177]);
ALGEBRAIC[182] = 1.00000/(1.00000+pow(CONSTANTS[71]/STATES[2], 2.00000));
ALGEBRAIC[183] = ( 3.00000*( ALGEBRAIC[181]*ALGEBRAIC[172] -  ALGEBRAIC[178]*ALGEBRAIC[173])+ ALGEBRAIC[180]*ALGEBRAIC[169]) -  ALGEBRAIC[179]*ALGEBRAIC[166];
ALGEBRAIC[184] =  ALGEBRAIC[179]*CONSTANTS[156] -  ALGEBRAIC[178]*CONSTANTS[155];
ALGEBRAIC[185] =  CONSTANTS[59]*CONSTANTS[151]*ALGEBRAIC[182]*( CONSTANTS[8]*ALGEBRAIC[183]+ CONSTANTS[9]*ALGEBRAIC[184]);
ALGEBRAIC[186] =  CONSTANTS[81]*exp(( CONSTANTS[83]*ALGEBRAIC[28])/3.00000);
ALGEBRAIC[187] =  CONSTANTS[82]*exp(( (1.00000 - CONSTANTS[83])*ALGEBRAIC[28])/3.00000);
ALGEBRAIC[188] = CONSTANTS[90]/(1.00000+CONSTANTS[89]/CONSTANTS[91]+STATES[3]/CONSTANTS[92]+STATES[5]/CONSTANTS[93]);
ALGEBRAIC[189] = ( CONSTANTS[73]*pow(STATES[3]/ALGEBRAIC[186], 3.00000))/((pow(1.00000+STATES[3]/ALGEBRAIC[186], 3.00000)+pow(1.00000+STATES[5]/CONSTANTS[84], 2.00000)) - 1.00000);
ALGEBRAIC[190] = ( CONSTANTS[76]*pow(CONSTANTS[1]/ALGEBRAIC[187], 3.00000))/((pow(1.00000+CONSTANTS[1]/ALGEBRAIC[187], 3.00000)+pow(1.00000+CONSTANTS[3]/CONSTANTS[85], 2.00000)) - 1.00000);
ALGEBRAIC[191] = ( CONSTANTS[77]*pow(CONSTANTS[3]/CONSTANTS[85], 2.00000))/((pow(1.00000+CONSTANTS[1]/ALGEBRAIC[187], 3.00000)+pow(1.00000+CONSTANTS[3]/CONSTANTS[85], 2.00000)) - 1.00000);
ALGEBRAIC[192] = ( CONSTANTS[78]*ALGEBRAIC[188]*CONSTANTS[89])/(1.00000+CONSTANTS[87]/CONSTANTS[88]);
ALGEBRAIC[193] = ( CONSTANTS[80]*pow(STATES[5]/CONSTANTS[84], 2.00000))/((pow(1.00000+STATES[3]/ALGEBRAIC[186], 3.00000)+pow(1.00000+STATES[5]/CONSTANTS[84], 2.00000)) - 1.00000);
ALGEBRAIC[194] =  CONSTANTS[160]*ALGEBRAIC[189]*CONSTANTS[159]+ ALGEBRAIC[190]*ALGEBRAIC[193]*ALGEBRAIC[192]+ CONSTANTS[159]*ALGEBRAIC[193]*ALGEBRAIC[192]+ ALGEBRAIC[192]*ALGEBRAIC[189]*CONSTANTS[159];
ALGEBRAIC[195] =  ALGEBRAIC[190]*CONSTANTS[158]*ALGEBRAIC[193]+ ALGEBRAIC[189]*CONSTANTS[159]*ALGEBRAIC[191]+ ALGEBRAIC[191]*CONSTANTS[158]*ALGEBRAIC[193]+ CONSTANTS[159]*ALGEBRAIC[191]*ALGEBRAIC[193];
ALGEBRAIC[196] =  CONSTANTS[159]*ALGEBRAIC[191]*CONSTANTS[160]+ ALGEBRAIC[192]*ALGEBRAIC[190]*CONSTANTS[158]+ ALGEBRAIC[190]*CONSTANTS[158]*CONSTANTS[160]+ ALGEBRAIC[191]*CONSTANTS[160]*CONSTANTS[158];
ALGEBRAIC[197] =  ALGEBRAIC[193]*ALGEBRAIC[192]*ALGEBRAIC[190]+ ALGEBRAIC[191]*CONSTANTS[160]*ALGEBRAIC[189]+ ALGEBRAIC[190]*CONSTANTS[160]*ALGEBRAIC[189]+ ALGEBRAIC[192]*ALGEBRAIC[190]*ALGEBRAIC[189];
ALGEBRAIC[198] = ALGEBRAIC[194]/(ALGEBRAIC[194]+ALGEBRAIC[195]+ALGEBRAIC[196]+ALGEBRAIC[197]);
ALGEBRAIC[199] = ALGEBRAIC[195]/(ALGEBRAIC[194]+ALGEBRAIC[195]+ALGEBRAIC[196]+ALGEBRAIC[197]);
ALGEBRAIC[200] = ALGEBRAIC[196]/(ALGEBRAIC[194]+ALGEBRAIC[195]+ALGEBRAIC[196]+ALGEBRAIC[197]);
ALGEBRAIC[201] = ALGEBRAIC[197]/(ALGEBRAIC[194]+ALGEBRAIC[195]+ALGEBRAIC[196]+ALGEBRAIC[197]);
ALGEBRAIC[202] =  3.00000*( ALGEBRAIC[198]*ALGEBRAIC[191] -  ALGEBRAIC[199]*ALGEBRAIC[192]);
ALGEBRAIC[203] =  2.00000*( ALGEBRAIC[201]*CONSTANTS[158] -  ALGEBRAIC[200]*ALGEBRAIC[189]);
ALGEBRAIC[204] =  CONSTANTS[161]*( CONSTANTS[8]*ALGEBRAIC[202]+ CONSTANTS[10]*ALGEBRAIC[203]);
ALGEBRAIC[205] = 1.00000/(1.00000+exp(- (STATES[0] - 10.8968)/23.9871));
ALGEBRAIC[206] =  CONSTANTS[125]*ALGEBRAIC[205]*(STATES[0] - ALGEBRAIC[64]);
ALGEBRAIC[207] = ( CONSTANTS[96]*ALGEBRAIC[25]*( STATES[3]*exp(ALGEBRAIC[28]) - CONSTANTS[1]))/(exp(ALGEBRAIC[28]) - 1.00000);
ALGEBRAIC[208] = (STATES[6] - STATES[5])/CONSTANTS[105];
ALGEBRAIC[209] = ( CONSTANTS[97]*4.00000*ALGEBRAIC[25]*( ALGEBRAIC[103]*STATES[9]*exp( 2.00000*ALGEBRAIC[28]) -  CONSTANTS[138]*CONSTANTS[2]))/(exp( 2.00000*ALGEBRAIC[28]) - 1.00000);ALGEBRAIC[214] =  (( (1.00000 - CONSTANTS[103])*CONSTANTS[100])/(1.00000+CONSTANTS[102]/STATES[9]))*(STATES[0] - ALGEBRAIC[66]);
ALGEBRAIC[210] = (STATES[4] - STATES[3])/CONSTANTS[104];
ALGEBRAIC[211] = ( CONSTANTS[98]*STATES[9])/(CONSTANTS[99]+STATES[9]);
ALGEBRAIC[212] =  (( CONSTANTS[103]*CONSTANTS[100])/(1.00000+CONSTANTS[102]/STATES[2]))*(STATES[0] - ALGEBRAIC[67]);
ALGEBRAIC[213] = (STATES[2] - STATES[9])/CONSTANTS[106];
ALGEBRAIC[214] =  (( (1.00000 - CONSTANTS[103])*CONSTANTS[100])/(1.00000+CONSTANTS[102]/STATES[9]))*(STATES[0] - ALGEBRAIC[66]);
ALGEBRAIC[215] = 1.00000/(1.00000+CONSTANTS[19]/ALGEBRAIC[49]);
ALGEBRAIC[216] = ALGEBRAIC[212]+ALGEBRAIC[214];
ALGEBRAIC[217] =  CONSTANTS[110]*( (1.00000 - ALGEBRAIC[215])*STATES[43]+ ALGEBRAIC[215]*STATES[44]);
ALGEBRAIC[218] =  CONSTANTS[101]*(STATES[0] - ALGEBRAIC[66]);
ALGEBRAIC[219] = ( CONSTANTS[128]*0.00542500*STATES[9])/(STATES[9]+0.000920000);
ALGEBRAIC[220] = ( CONSTANTS[128]*2.75000*0.00542500*STATES[9])/((STATES[9]+0.000920000) - 0.000170000);
ALGEBRAIC[221] = (STATES[11] - STATES[10])/CONSTANTS[104];
ALGEBRAIC[222] = 1.00000/(1.00000+CONSTANTS[19]/ALGEBRAIC[49]);
ALGEBRAIC[223] = ( 0.00488250*STATES[7])/15.0000;
ALGEBRAIC[224] =  CONSTANTS[111]*(( (1.00000 - ALGEBRAIC[222])*ALGEBRAIC[219]+ ALGEBRAIC[222]*ALGEBRAIC[220]) - ALGEBRAIC[223]);
ALGEBRAIC[225] = (STATES[7] - STATES[8])/60.0000;
}
void
computeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC)
{
ALGEBRAIC[3] = 1.00000/(1.00000+exp((STATES[0]+87.6100)/7.48800));
ALGEBRAIC[4] = 1.00000/(1.00000+exp((STATES[0]+93.8100)/7.48800));
ALGEBRAIC[9] = 1.00000/(1.00000+exp((STATES[0]+18.0800)/2.79160));
ALGEBRAIC[0] = 1.00000/pow(1.00000+exp(- (STATES[0]+56.8600)/9.03000), 2.00000);
ALGEBRAIC[13] =  0.129200*exp(- pow((STATES[0]+45.7900)/15.5400, 2.00000))+ 0.0648700*exp(- pow((STATES[0] - 4.82300)/51.1200, 2.00000));
ALGEBRAIC[2] = 1.00000/(1.00000+exp(- (STATES[0]+42.8500)/5.26400));
ALGEBRAIC[16] =  0.129200*exp(- pow((STATES[0]+45.7900)/15.5400, 2.00000))+ 0.0648700*exp(- pow((STATES[0] - 4.82300)/51.1200, 2.00000));
ALGEBRAIC[5] = 1.00000/(1.00000+exp(- ((STATES[0]+CONSTANTS[44]) - 14.3400)/14.8200));
ALGEBRAIC[17] = 1.05150/(1.00000/( 1.20890*(1.00000+exp(- ((STATES[0]+CONSTANTS[44]) - 18.4099)/29.3814)))+3.50000/(1.00000+exp((STATES[0]+CONSTANTS[44]+100.000)/29.3814)));
ALGEBRAIC[7] = (STATES[0]>=31.4978 ? 1.00000 :  1.07630*exp( - 1.00700*exp( - 0.0829000*STATES[0])));
ALGEBRAIC[22] = CONSTANTS[51]+0.600000+1.00000/(exp( - 0.0500000*(STATES[0]+CONSTANTS[50]+6.00000))+exp( 0.0900000*(STATES[0]+CONSTANTS[50]+14.0000)));
ALGEBRAIC[8] = 1.00000/(1.00000+exp((STATES[0]+19.5800)/3.69600));
ALGEBRAIC[23] = 7.00000+1.00000/( 0.00450000*exp(- (STATES[0]+20.0000)/10.0000)+ 0.00450000*exp((STATES[0]+20.0000)/10.0000));
ALGEBRAIC[24] = 1000.00+1.00000/( 3.50000e-05*exp(- (STATES[0]+5.00000)/4.00000)+ 3.50000e-05*exp((STATES[0]+5.00000)/6.00000));
ALGEBRAIC[10] =  STATES[31]*1.00000;
ALGEBRAIC[20] = 1.00000/(CONSTANTS[46]/ALGEBRAIC[10]+pow(1.00000+CONSTANTS[45]/STATES[2], 4.00000));
ALGEBRAIC[21] = 1.00000/(CONSTANTS[46]/ALGEBRAIC[10]+pow(1.00000+CONSTANTS[45]/STATES[9], 4.00000));
ALGEBRAIC[12] = 1.00000/(1.00000+exp(- (STATES[0]+11.6000)/8.93200));
ALGEBRAIC[27] = 817.300+1.00000/( 0.000232600*exp((STATES[0]+48.2800)/17.8000)+ 0.00129200*exp(- (STATES[0]+210.000)/230.000));
ALGEBRAIC[32] = 1.00000/(1.00000+exp(- ((STATES[0]+CONSTANTS[44]) - 24.3400)/14.8200));
ALGEBRAIC[19] = ALGEBRAIC[8];
ALGEBRAIC[33] = 7.00000+1.00000/( 0.0400000*exp(- (STATES[0] - 4.00000)/7.00000)+ 0.0400000*exp((STATES[0] - 4.00000)/7.00000));
ALGEBRAIC[34] = 100.000+1.00000/( 0.000120000*exp(- STATES[0]/3.00000)+ 0.000120000*exp(STATES[0]/7.00000));
ALGEBRAIC[35] =  2.50000*ALGEBRAIC[23];
ALGEBRAIC[26] = ALGEBRAIC[12];
ALGEBRAIC[36] = 1.00000/( 0.0100000*exp((STATES[0] - 50.0000)/20.0000)+ 0.0193000*exp(- (STATES[0]+66.5400)/31.0000));
ALGEBRAIC[43] = ( CONSTANTS[22]*(1.00000 - STATES[1]))/(1.00000+CONSTANTS[23]/STATES[2]);
ALGEBRAIC[1] = 1.00000/pow(1.00000+exp((STATES[0]+71.5500)/7.43000), 2.00000);
ALGEBRAIC[14] = (STATES[0]>=- 40.0000 ? 0.00000 :  0.0570000*exp(- (STATES[0]+80.0000)/6.80000));
ALGEBRAIC[29] = (STATES[0]>=- 40.0000 ? 0.770000/( 0.130000*(1.00000+exp(- (STATES[0]+10.6600)/11.1000))) :  2.70000*exp( 0.0790000*STATES[0])+ 310000.*exp( 0.348500*STATES[0]));
ALGEBRAIC[37] = 1.00000/(ALGEBRAIC[14]+ALGEBRAIC[29]);
ALGEBRAIC[40] =  2.50000*ALGEBRAIC[33];
ALGEBRAIC[38] = ALGEBRAIC[1];
ALGEBRAIC[15] = (STATES[0]>=- 40.0000 ? 0.00000 : ( ( - 25428.0*exp( 0.244400*STATES[0]) -  6.94800e-06*exp( - 0.0439100*STATES[0]))*(STATES[0]+37.7800))/(1.00000+exp( 0.311000*(STATES[0]+79.2300))));
ALGEBRAIC[30] = (STATES[0]>=- 40.0000 ? ( 0.600000*exp( 0.0570000*STATES[0]))/(1.00000+exp( - 0.100000*(STATES[0]+32.0000))) : ( 0.0242400*exp( - 0.0105200*STATES[0]))/(1.00000+exp( - 0.137800*(STATES[0]+40.1400))));
ALGEBRAIC[44] = 1.00000/(ALGEBRAIC[15]+ALGEBRAIC[30]);
ALGEBRAIC[45] = 1.00000/pow(1.00000+exp((STATES[0]+77.5500)/7.43000), 2.00000);
ALGEBRAIC[6] = 1.00000/(1.00000+exp((STATES[0]+CONSTANTS[44]+43.9400)/5.71100));
ALGEBRAIC[18] = (CONSTANTS[0]==1.00000 ? 1.00000 - 0.950000/(1.00000+exp((STATES[0]+CONSTANTS[44]+70.0000)/5.00000)) : 1.00000);
ALGEBRAIC[31] = 4.56200+1.00000/( 0.393300*exp(- (STATES[0]+CONSTANTS[44]+100.000)/100.000)+ 0.0800400*exp((STATES[0]+CONSTANTS[44]+50.0000)/16.5900));
ALGEBRAIC[46] =  ALGEBRAIC[31]*ALGEBRAIC[18];
ALGEBRAIC[28] = ( STATES[0]*CONSTANTS[7])/( CONSTANTS[5]*CONSTANTS[6]);
ALGEBRAIC[41] =  0.116100*exp( 0.299000*ALGEBRAIC[28]);
ALGEBRAIC[47] =  0.244200*exp( - 1.60400*ALGEBRAIC[28]);
ALGEBRAIC[50] =  1.46000*ALGEBRAIC[44];
ALGEBRAIC[39] = 23.6200+1.00000/( 0.00141600*exp(- (STATES[0]+CONSTANTS[44]+96.5200)/59.0500)+ 1.78000e-08*exp((STATES[0]+CONSTANTS[44]+114.100)/8.07900));
ALGEBRAIC[51] =  ALGEBRAIC[39]*ALGEBRAIC[18];
ALGEBRAIC[42] =  0.0578000*exp( 0.971000*ALGEBRAIC[28]);
ALGEBRAIC[48] =  0.000349000*exp( - 1.06200*ALGEBRAIC[28]);
ALGEBRAIC[52] =  0.253300*exp( 0.595300*ALGEBRAIC[28]);
ALGEBRAIC[55] =  0.0652500*exp( - 0.820900*ALGEBRAIC[28]);
ALGEBRAIC[54] = 1.35400+0.000100000/(exp(((STATES[0]+CONSTANTS[44]) - 167.400)/15.8900)+exp(- ((STATES[0]+CONSTANTS[44]) - 12.2300)/0.215400));
ALGEBRAIC[57] = 1.00000 - 0.500000/(1.00000+exp((STATES[0]+CONSTANTS[44]+70.0000)/20.0000));
ALGEBRAIC[60] =  ALGEBRAIC[54]*ALGEBRAIC[57]*ALGEBRAIC[46];
ALGEBRAIC[61] =  ALGEBRAIC[54]*ALGEBRAIC[57]*ALGEBRAIC[51];
ALGEBRAIC[58] =  5.20000e-05*exp( 1.52500*ALGEBRAIC[28]);
ALGEBRAIC[62] = ( ALGEBRAIC[48]*ALGEBRAIC[55]*ALGEBRAIC[58])/( ALGEBRAIC[42]*ALGEBRAIC[52]);
ALGEBRAIC[79] =  CONSTANTS[48]*STATES[27]+ CONSTANTS[119]*STATES[28];
ALGEBRAIC[80] = 0.300000+0.600000/(1.00000+exp((STATES[0] - 10.0000)/10.0000));
ALGEBRAIC[81] = 1.00000 - ALGEBRAIC[80];
ALGEBRAIC[82] =  ALGEBRAIC[80]*STATES[29]+ ALGEBRAIC[81]*STATES[30];
ALGEBRAIC[83] =  CONSTANTS[48]*STATES[32]+ CONSTANTS[119]*STATES[28];
ALGEBRAIC[84] =  ALGEBRAIC[80]*STATES[33]+ ALGEBRAIC[81]*STATES[30];
ALGEBRAIC[25] = ( STATES[0]*CONSTANTS[7]*CONSTANTS[7])/( CONSTANTS[5]*CONSTANTS[6]);
ALGEBRAIC[85] = ( 0.500000*(STATES[4]+STATES[6]+STATES[11]+ 4.00000*STATES[2]))/1000.00;
ALGEBRAIC[86] = exp( - CONSTANTS[133]*4.00000*( pow(ALGEBRAIC[85], 1.0 / 2)/(1.00000+ pow(ALGEBRAIC[85], 1.0 / 2)) -  0.300000*ALGEBRAIC[85]));
ALGEBRAIC[89] = ( 4.00000*ALGEBRAIC[25]*( ALGEBRAIC[86]*STATES[2]*exp( 2.00000*ALGEBRAIC[28]) -  CONSTANTS[138]*CONSTANTS[2]))/(exp( 2.00000*ALGEBRAIC[28]) - 1.00000);
ALGEBRAIC[49] = ALGEBRAIC[43]+STATES[1];
ALGEBRAIC[92] = 1.00000/(1.00000+CONSTANTS[19]/ALGEBRAIC[49]);
ALGEBRAIC[93] =  CONSTANTS[53]*( (1.00000 - ALGEBRAIC[92])*CONSTANTS[120]*ALGEBRAIC[89]*STATES[26]*( ALGEBRAIC[79]*(1.00000 - STATES[34])+ STATES[31]*ALGEBRAIC[82]*STATES[34])+ ALGEBRAIC[92]*CONSTANTS[130]*ALGEBRAIC[89]*STATES[26]*( ALGEBRAIC[83]*(1.00000 - STATES[34])+ STATES[31]*ALGEBRAIC[84]*STATES[34]));
ALGEBRAIC[95] = (( - CONSTANTS[126]*ALGEBRAIC[93])/1.00000)/(1.00000+pow(CONSTANTS[109]/STATES[8], 8.00000));
ALGEBRAIC[98] = (CONSTANTS[0]==2.00000 ?  ALGEBRAIC[95]*1.70000 : ALGEBRAIC[95]);
ALGEBRAIC[101] = CONSTANTS[108]/(1.00000+0.0123000/STATES[8]);
ALGEBRAIC[104] = (ALGEBRAIC[101]<0.00100000 ? 0.00100000 : ALGEBRAIC[101]);
ALGEBRAIC[96] = (( - CONSTANTS[134]*ALGEBRAIC[93])/1.00000)/(1.00000+pow(CONSTANTS[109]/STATES[8], 8.00000));
ALGEBRAIC[99] = (CONSTANTS[0]==2.00000 ?  ALGEBRAIC[96]*1.70000 : ALGEBRAIC[96]);
ALGEBRAIC[102] = CONSTANTS[127]/(1.00000+0.0123000/STATES[8]);
ALGEBRAIC[105] = (ALGEBRAIC[102]<0.00100000 ? 0.00100000 : ALGEBRAIC[102]);
ALGEBRAIC[64] =  (( CONSTANTS[5]*CONSTANTS[6])/( CONSTANTS[10]*CONSTANTS[7]))*log(CONSTANTS[3]/STATES[5]);
ALGEBRAIC[73] = 1.00000/(1.00000+exp(((STATES[0]+CONSTANTS[44]) - 213.600)/151.200));
ALGEBRAIC[74] = 1.00000 - ALGEBRAIC[73];
ALGEBRAIC[75] =  ALGEBRAIC[73]*STATES[21]+ ALGEBRAIC[74]*STATES[22];
ALGEBRAIC[76] =  ALGEBRAIC[73]*STATES[24]+ ALGEBRAIC[74]*STATES[25];
ALGEBRAIC[77] = 1.00000/(1.00000+CONSTANTS[19]/ALGEBRAIC[49]);
ALGEBRAIC[78] =  CONSTANTS[118]*(STATES[0] - ALGEBRAIC[64])*( (1.00000 - ALGEBRAIC[77])*STATES[20]*ALGEBRAIC[75]+ ALGEBRAIC[77]*STATES[23]*ALGEBRAIC[76]);
ALGEBRAIC[117] =  CONSTANTS[122]* pow((CONSTANTS[3]/5.00000), 1.0 / 2)*STATES[40]*(STATES[0] - ALGEBRAIC[64]);
ALGEBRAIC[65] =  (( CONSTANTS[5]*CONSTANTS[6])/( CONSTANTS[10]*CONSTANTS[7]))*log((CONSTANTS[3]+ CONSTANTS[34]*CONSTANTS[1])/(STATES[5]+ CONSTANTS[34]*STATES[3]));
ALGEBRAIC[118] = 1.00000+0.600000/(1.00000+pow(3.80000e-05/STATES[9], 1.40000));
ALGEBRAIC[119] =  CONSTANTS[123]*ALGEBRAIC[118]*STATES[41]*STATES[42]*(STATES[0] - ALGEBRAIC[65]);
ALGEBRAIC[120] = 4.09400/(1.00000+exp( 0.121700*((STATES[0] - ALGEBRAIC[64]) - 49.9340)));
ALGEBRAIC[121] = ( 15.7200*exp( 0.0674000*((STATES[0] - ALGEBRAIC[64]) - 3.25700))+exp( 0.0618000*((STATES[0] - ALGEBRAIC[64]) - 594.310)))/(1.00000+exp( - 0.162900*((STATES[0] - ALGEBRAIC[64])+14.2070)));
ALGEBRAIC[122] = ALGEBRAIC[120]/(ALGEBRAIC[120]+ALGEBRAIC[121]);
ALGEBRAIC[123] =  CONSTANTS[124]* pow((CONSTANTS[3]/5.00000), 1.0 / 2)*ALGEBRAIC[122]*(STATES[0] - ALGEBRAIC[64]);
ALGEBRAIC[187] =  CONSTANTS[82]*exp(( (1.00000 - CONSTANTS[83])*ALGEBRAIC[28])/3.00000);
ALGEBRAIC[191] = ( CONSTANTS[77]*pow(CONSTANTS[3]/CONSTANTS[85], 2.00000))/((pow(1.00000+CONSTANTS[1]/ALGEBRAIC[187], 3.00000)+pow(1.00000+CONSTANTS[3]/CONSTANTS[85], 2.00000)) - 1.00000);
ALGEBRAIC[188] = CONSTANTS[90]/(1.00000+CONSTANTS[89]/CONSTANTS[91]+STATES[3]/CONSTANTS[92]+STATES[5]/CONSTANTS[93]);
ALGEBRAIC[192] = ( CONSTANTS[78]*ALGEBRAIC[188]*CONSTANTS[89])/(1.00000+CONSTANTS[87]/CONSTANTS[88]);
ALGEBRAIC[186] =  CONSTANTS[81]*exp(( CONSTANTS[83]*ALGEBRAIC[28])/3.00000);
ALGEBRAIC[189] = ( CONSTANTS[73]*pow(STATES[3]/ALGEBRAIC[186], 3.00000))/((pow(1.00000+STATES[3]/ALGEBRAIC[186], 3.00000)+pow(1.00000+STATES[5]/CONSTANTS[84], 2.00000)) - 1.00000);
ALGEBRAIC[190] = ( CONSTANTS[76]*pow(CONSTANTS[1]/ALGEBRAIC[187], 3.00000))/((pow(1.00000+CONSTANTS[1]/ALGEBRAIC[187], 3.00000)+pow(1.00000+CONSTANTS[3]/CONSTANTS[85], 2.00000)) - 1.00000);
ALGEBRAIC[193] = ( CONSTANTS[80]*pow(STATES[5]/CONSTANTS[84], 2.00000))/((pow(1.00000+STATES[3]/ALGEBRAIC[186], 3.00000)+pow(1.00000+STATES[5]/CONSTANTS[84], 2.00000)) - 1.00000);
ALGEBRAIC[194] =  CONSTANTS[160]*ALGEBRAIC[189]*CONSTANTS[159]+ ALGEBRAIC[190]*ALGEBRAIC[193]*ALGEBRAIC[192]+ CONSTANTS[159]*ALGEBRAIC[193]*ALGEBRAIC[192]+ ALGEBRAIC[192]*ALGEBRAIC[189]*CONSTANTS[159];
ALGEBRAIC[195] =  ALGEBRAIC[190]*CONSTANTS[158]*ALGEBRAIC[193]+ ALGEBRAIC[189]*CONSTANTS[159]*ALGEBRAIC[191]+ ALGEBRAIC[191]*CONSTANTS[158]*ALGEBRAIC[193]+ CONSTANTS[159]*ALGEBRAIC[191]*ALGEBRAIC[193];
ALGEBRAIC[196] =  CONSTANTS[159]*ALGEBRAIC[191]*CONSTANTS[160]+ ALGEBRAIC[192]*ALGEBRAIC[190]*CONSTANTS[158]+ ALGEBRAIC[190]*CONSTANTS[158]*CONSTANTS[160]+ ALGEBRAIC[191]*CONSTANTS[160]*CONSTANTS[158];
ALGEBRAIC[197] =  ALGEBRAIC[193]*ALGEBRAIC[192]*ALGEBRAIC[190]+ ALGEBRAIC[191]*CONSTANTS[160]*ALGEBRAIC[189]+ ALGEBRAIC[190]*CONSTANTS[160]*ALGEBRAIC[189]+ ALGEBRAIC[192]*ALGEBRAIC[190]*ALGEBRAIC[189];
ALGEBRAIC[198] = ALGEBRAIC[194]/(ALGEBRAIC[194]+ALGEBRAIC[195]+ALGEBRAIC[196]+ALGEBRAIC[197]);
ALGEBRAIC[199] = ALGEBRAIC[195]/(ALGEBRAIC[194]+ALGEBRAIC[195]+ALGEBRAIC[196]+ALGEBRAIC[197]);
ALGEBRAIC[202] =  3.00000*( ALGEBRAIC[198]*ALGEBRAIC[191] -  ALGEBRAIC[199]*ALGEBRAIC[192]);
ALGEBRAIC[200] = ALGEBRAIC[196]/(ALGEBRAIC[194]+ALGEBRAIC[195]+ALGEBRAIC[196]+ALGEBRAIC[197]);
ALGEBRAIC[201] = ALGEBRAIC[197]/(ALGEBRAIC[194]+ALGEBRAIC[195]+ALGEBRAIC[196]+ALGEBRAIC[197]);
ALGEBRAIC[203] =  2.00000*( ALGEBRAIC[201]*CONSTANTS[158] -  ALGEBRAIC[200]*ALGEBRAIC[189]);
ALGEBRAIC[204] =  CONSTANTS[161]*( CONSTANTS[8]*ALGEBRAIC[202]+ CONSTANTS[10]*ALGEBRAIC[203]);
ALGEBRAIC[205] = 1.00000/(1.00000+exp(- (STATES[0] - 10.8968)/23.9871));
ALGEBRAIC[206] =  CONSTANTS[125]*ALGEBRAIC[205]*(STATES[0] - ALGEBRAIC[64]);
ALGEBRAIC[68] =  CONSTANTS[36]*CONSTANTS[35]*CONSTANTS[114]*CONSTANTS[115]*(STATES[0] - ALGEBRAIC[64]);
ALGEBRAIC[11] = (VOI>=CONSTANTS[14]&&(VOI - CONSTANTS[14]) -  floor((VOI - CONSTANTS[14])/CONSTANTS[17])*CONSTANTS[17]<=CONSTANTS[18] ? CONSTANTS[16] : 0.00000);
ALGEBRAIC[100] = ( 0.500000*(STATES[3]+STATES[5]+STATES[10]+ 4.00000*STATES[9]))/1000.00;
ALGEBRAIC[107] = exp( - CONSTANTS[133]*1.00000*( pow(ALGEBRAIC[100], 1.0 / 2)/(1.00000+ pow(ALGEBRAIC[100], 1.0 / 2)) -  0.300000*ALGEBRAIC[100]));
ALGEBRAIC[110] = ( 1.00000*ALGEBRAIC[25]*( ALGEBRAIC[107]*STATES[5]*exp( 1.00000*ALGEBRAIC[28]) -  CONSTANTS[140]*CONSTANTS[3]))/(exp( 1.00000*ALGEBRAIC[28]) - 1.00000);
ALGEBRAIC[113] =  (1.00000 - CONSTANTS[53])*( (1.00000 - ALGEBRAIC[92])*CONSTANTS[132]*ALGEBRAIC[110]*STATES[26]*( ALGEBRAIC[79]*(1.00000 - STATES[35])+ STATES[31]*ALGEBRAIC[82]*STATES[35])+ ALGEBRAIC[92]*CONSTANTS[137]*ALGEBRAIC[110]*STATES[26]*( ALGEBRAIC[83]*(1.00000 - STATES[35])+ STATES[31]*ALGEBRAIC[84]*STATES[35]));
ALGEBRAIC[208] = (STATES[6] - STATES[5])/CONSTANTS[105];
ALGEBRAIC[88] = exp( - CONSTANTS[133]*1.00000*( pow(ALGEBRAIC[85], 1.0 / 2)/(1.00000+ pow(ALGEBRAIC[85], 1.0 / 2)) -  0.300000*ALGEBRAIC[85]));
ALGEBRAIC[91] = ( 1.00000*ALGEBRAIC[25]*( ALGEBRAIC[88]*STATES[6]*exp( 1.00000*ALGEBRAIC[28]) -  CONSTANTS[140]*CONSTANTS[3]))/(exp( 1.00000*ALGEBRAIC[28]) - 1.00000);
ALGEBRAIC[97] =  CONSTANTS[53]*( (1.00000 - ALGEBRAIC[92])*CONSTANTS[132]*ALGEBRAIC[91]*STATES[26]*( ALGEBRAIC[79]*(1.00000 - STATES[34])+ STATES[31]*ALGEBRAIC[82]*STATES[34])+ ALGEBRAIC[92]*CONSTANTS[137]*ALGEBRAIC[91]*STATES[26]*( ALGEBRAIC[83]*(1.00000 - STATES[34])+ STATES[31]*ALGEBRAIC[84]*STATES[34]));
ALGEBRAIC[63] =  (( CONSTANTS[5]*CONSTANTS[6])/( CONSTANTS[8]*CONSTANTS[7]))*log(CONSTANTS[1]/STATES[3]);
ALGEBRAIC[69] = 1.00000/(1.00000+CONSTANTS[19]/ALGEBRAIC[49]);
ALGEBRAIC[70] =  CONSTANTS[40]*(STATES[0] - ALGEBRAIC[63])*pow(STATES[12], 3.00000)*( (1.00000 - ALGEBRAIC[69])*STATES[13]*STATES[14]+ ALGEBRAIC[69]*STATES[15]*STATES[16]);
ALGEBRAIC[71] = 1.00000/(1.00000+CONSTANTS[19]/ALGEBRAIC[49]);
ALGEBRAIC[72] =  CONSTANTS[117]*(STATES[0] - ALGEBRAIC[63])*STATES[17]*( (1.00000 - ALGEBRAIC[71])*STATES[18]+ ALGEBRAIC[71]*STATES[19]);
ALGEBRAIC[152] = 1.00000/(1.00000+pow(CONSTANTS[71]/STATES[9], 2.00000));
ALGEBRAIC[125] = exp( CONSTANTS[69]*ALGEBRAIC[28]);
ALGEBRAIC[132] = 1.00000+ (CONSTANTS[1]/CONSTANTS[62])*(1.00000+1.00000/ALGEBRAIC[125]);
ALGEBRAIC[133] = CONSTANTS[1]/( CONSTANTS[62]*ALGEBRAIC[125]*ALGEBRAIC[132]);
ALGEBRAIC[136] =  ALGEBRAIC[133]*CONSTANTS[66];
ALGEBRAIC[126] = 1.00000+ (STATES[3]/CONSTANTS[62])*(1.00000+ALGEBRAIC[125]);
ALGEBRAIC[127] = ( STATES[3]*ALGEBRAIC[125])/( CONSTANTS[62]*ALGEBRAIC[126]);
ALGEBRAIC[139] =  ALGEBRAIC[127]*CONSTANTS[66];
ALGEBRAIC[129] = 1.00000+ (STATES[3]/CONSTANTS[60])*(1.00000+STATES[3]/CONSTANTS[61]);
ALGEBRAIC[130] = ( STATES[3]*STATES[3])/( ALGEBRAIC[129]*CONSTANTS[60]*CONSTANTS[61]);
ALGEBRAIC[142] =  ALGEBRAIC[130]*ALGEBRAIC[127]*CONSTANTS[64];
ALGEBRAIC[143] =  ALGEBRAIC[133]*CONSTANTS[146]*CONSTANTS[64];
ALGEBRAIC[134] = 1.00000/ALGEBRAIC[132];
ALGEBRAIC[135] =  ALGEBRAIC[134]*CONSTANTS[65];
ALGEBRAIC[137] = ALGEBRAIC[135]+ALGEBRAIC[136];
ALGEBRAIC[124] = exp( CONSTANTS[70]*ALGEBRAIC[28]);
ALGEBRAIC[128] = 1.00000/ALGEBRAIC[126];
ALGEBRAIC[138] = ( ALGEBRAIC[128]*CONSTANTS[65])/ALGEBRAIC[124];
ALGEBRAIC[140] = ALGEBRAIC[138]+ALGEBRAIC[139];
ALGEBRAIC[131] = 1.00000/ALGEBRAIC[129];
ALGEBRAIC[141] =  ALGEBRAIC[131]*STATES[9]*CONSTANTS[67];
ALGEBRAIC[144] =  CONSTANTS[149]*ALGEBRAIC[140]*(ALGEBRAIC[142]+ALGEBRAIC[141])+ CONSTANTS[150]*ALGEBRAIC[142]*(CONSTANTS[149]+ALGEBRAIC[137]);
ALGEBRAIC[145] =  CONSTANTS[148]*ALGEBRAIC[142]*(ALGEBRAIC[140]+CONSTANTS[150])+ ALGEBRAIC[140]*ALGEBRAIC[141]*(CONSTANTS[148]+ALGEBRAIC[143]);
ALGEBRAIC[146] =  CONSTANTS[148]*ALGEBRAIC[137]*(ALGEBRAIC[142]+ALGEBRAIC[141])+ ALGEBRAIC[143]*ALGEBRAIC[141]*(CONSTANTS[149]+ALGEBRAIC[137]);
ALGEBRAIC[147] =  CONSTANTS[149]*ALGEBRAIC[143]*(ALGEBRAIC[140]+CONSTANTS[150])+ ALGEBRAIC[137]*CONSTANTS[150]*(CONSTANTS[148]+ALGEBRAIC[143]);
ALGEBRAIC[148] = ALGEBRAIC[144]/(ALGEBRAIC[144]+ALGEBRAIC[145]+ALGEBRAIC[146]+ALGEBRAIC[147]);
ALGEBRAIC[149] = ALGEBRAIC[145]/(ALGEBRAIC[144]+ALGEBRAIC[145]+ALGEBRAIC[146]+ALGEBRAIC[147]);
ALGEBRAIC[150] = ALGEBRAIC[146]/(ALGEBRAIC[144]+ALGEBRAIC[145]+ALGEBRAIC[146]+ALGEBRAIC[147]);
ALGEBRAIC[151] = ALGEBRAIC[147]/(ALGEBRAIC[144]+ALGEBRAIC[145]+ALGEBRAIC[146]+ALGEBRAIC[147]);
ALGEBRAIC[153] = ( 3.00000*( ALGEBRAIC[151]*ALGEBRAIC[142] -  ALGEBRAIC[148]*ALGEBRAIC[143])+ ALGEBRAIC[150]*ALGEBRAIC[139]) -  ALGEBRAIC[149]*ALGEBRAIC[136];
ALGEBRAIC[154] =  ALGEBRAIC[149]*CONSTANTS[149] -  ALGEBRAIC[148]*CONSTANTS[148];
ALGEBRAIC[155] =  (1.00000 - CONSTANTS[59])*CONSTANTS[151]*ALGEBRAIC[152]*( CONSTANTS[8]*ALGEBRAIC[153]+ CONSTANTS[9]*ALGEBRAIC[154]);
ALGEBRAIC[207] = ( CONSTANTS[96]*ALGEBRAIC[25]*( STATES[3]*exp(ALGEBRAIC[28]) - CONSTANTS[1]))/(exp(ALGEBRAIC[28]) - 1.00000);
ALGEBRAIC[106] = exp( - CONSTANTS[133]*1.00000*( pow(ALGEBRAIC[100], 1.0 / 2)/(1.00000+ pow(ALGEBRAIC[100], 1.0 / 2)) -  0.300000*ALGEBRAIC[100]));
ALGEBRAIC[109] = ( 1.00000*ALGEBRAIC[25]*( ALGEBRAIC[106]*STATES[3]*exp( 1.00000*ALGEBRAIC[28]) -  CONSTANTS[139]*CONSTANTS[1]))/(exp( 1.00000*ALGEBRAIC[28]) - 1.00000);
ALGEBRAIC[112] =  (1.00000 - CONSTANTS[53])*( (1.00000 - ALGEBRAIC[92])*CONSTANTS[131]*ALGEBRAIC[109]*STATES[26]*( ALGEBRAIC[79]*(1.00000 - STATES[35])+ STATES[31]*ALGEBRAIC[82]*STATES[35])+ ALGEBRAIC[92]*CONSTANTS[136]*ALGEBRAIC[109]*STATES[26]*( ALGEBRAIC[83]*(1.00000 - STATES[35])+ STATES[31]*ALGEBRAIC[84]*STATES[35]));
ALGEBRAIC[210] = (STATES[4] - STATES[3])/CONSTANTS[104];
ALGEBRAIC[182] = 1.00000/(1.00000+pow(CONSTANTS[71]/STATES[2], 2.00000));
ALGEBRAIC[162] = 1.00000+ (CONSTANTS[1]/CONSTANTS[62])*(1.00000+1.00000/ALGEBRAIC[125]);
ALGEBRAIC[163] = CONSTANTS[1]/( CONSTANTS[62]*ALGEBRAIC[125]*ALGEBRAIC[162]);
ALGEBRAIC[166] =  ALGEBRAIC[163]*CONSTANTS[66];
ALGEBRAIC[156] = 1.00000+ (STATES[4]/CONSTANTS[62])*(1.00000+ALGEBRAIC[125]);
ALGEBRAIC[157] = ( STATES[4]*ALGEBRAIC[125])/( CONSTANTS[62]*ALGEBRAIC[156]);
ALGEBRAIC[169] =  ALGEBRAIC[157]*CONSTANTS[66];
ALGEBRAIC[159] = 1.00000+ (STATES[4]/CONSTANTS[60])*(1.00000+STATES[4]/CONSTANTS[61]);
ALGEBRAIC[160] = ( STATES[4]*STATES[4])/( ALGEBRAIC[159]*CONSTANTS[60]*CONSTANTS[61]);
ALGEBRAIC[172] =  ALGEBRAIC[160]*ALGEBRAIC[157]*CONSTANTS[64];
ALGEBRAIC[173] =  ALGEBRAIC[163]*CONSTANTS[153]*CONSTANTS[64];
ALGEBRAIC[164] = 1.00000/ALGEBRAIC[162];
ALGEBRAIC[165] =  ALGEBRAIC[164]*CONSTANTS[65];
ALGEBRAIC[167] = ALGEBRAIC[165]+ALGEBRAIC[166];
ALGEBRAIC[158] = 1.00000/ALGEBRAIC[156];
ALGEBRAIC[168] = ( ALGEBRAIC[158]*CONSTANTS[65])/ALGEBRAIC[124];
ALGEBRAIC[170] = ALGEBRAIC[168]+ALGEBRAIC[169];
ALGEBRAIC[161] = 1.00000/ALGEBRAIC[159];
ALGEBRAIC[171] =  ALGEBRAIC[161]*STATES[2]*CONSTANTS[67];
ALGEBRAIC[174] =  CONSTANTS[156]*ALGEBRAIC[170]*(ALGEBRAIC[172]+ALGEBRAIC[171])+ CONSTANTS[157]*ALGEBRAIC[172]*(CONSTANTS[156]+ALGEBRAIC[167]);
ALGEBRAIC[175] =  CONSTANTS[155]*ALGEBRAIC[172]*(ALGEBRAIC[170]+CONSTANTS[157])+ ALGEBRAIC[170]*ALGEBRAIC[171]*(CONSTANTS[155]+ALGEBRAIC[173]);
ALGEBRAIC[176] =  CONSTANTS[155]*ALGEBRAIC[167]*(ALGEBRAIC[172]+ALGEBRAIC[171])+ ALGEBRAIC[173]*ALGEBRAIC[171]*(CONSTANTS[156]+ALGEBRAIC[167]);
ALGEBRAIC[177] =  CONSTANTS[156]*ALGEBRAIC[173]*(ALGEBRAIC[170]+CONSTANTS[157])+ ALGEBRAIC[167]*CONSTANTS[157]*(CONSTANTS[155]+ALGEBRAIC[173]);
ALGEBRAIC[178] = ALGEBRAIC[174]/(ALGEBRAIC[174]+ALGEBRAIC[175]+ALGEBRAIC[176]+ALGEBRAIC[177]);
ALGEBRAIC[179] = ALGEBRAIC[175]/(ALGEBRAIC[174]+ALGEBRAIC[175]+ALGEBRAIC[176]+ALGEBRAIC[177]);
ALGEBRAIC[180] = ALGEBRAIC[176]/(ALGEBRAIC[174]+ALGEBRAIC[175]+ALGEBRAIC[176]+ALGEBRAIC[177]);
ALGEBRAIC[181] = ALGEBRAIC[177]/(ALGEBRAIC[174]+ALGEBRAIC[175]+ALGEBRAIC[176]+ALGEBRAIC[177]);
ALGEBRAIC[183] = ( 3.00000*( ALGEBRAIC[181]*ALGEBRAIC[172] -  ALGEBRAIC[178]*ALGEBRAIC[173])+ ALGEBRAIC[180]*ALGEBRAIC[169]) -  ALGEBRAIC[179]*ALGEBRAIC[166];
ALGEBRAIC[184] =  ALGEBRAIC[179]*CONSTANTS[156] -  ALGEBRAIC[178]*CONSTANTS[155];
ALGEBRAIC[185] =  CONSTANTS[59]*CONSTANTS[151]*ALGEBRAIC[182]*( CONSTANTS[8]*ALGEBRAIC[183]+ CONSTANTS[9]*ALGEBRAIC[184]);
ALGEBRAIC[87] = exp( - CONSTANTS[133]*1.00000*( pow(ALGEBRAIC[85], 1.0 / 2)/(1.00000+ pow(ALGEBRAIC[85], 1.0 / 2)) -  0.300000*ALGEBRAIC[85]));
ALGEBRAIC[90] = ( 1.00000*ALGEBRAIC[25]*( ALGEBRAIC[87]*STATES[4]*exp( 1.00000*ALGEBRAIC[28]) -  CONSTANTS[139]*CONSTANTS[1]))/(exp( 1.00000*ALGEBRAIC[28]) - 1.00000);
ALGEBRAIC[94] =  CONSTANTS[53]*( (1.00000 - ALGEBRAIC[92])*CONSTANTS[131]*ALGEBRAIC[90]*STATES[26]*( ALGEBRAIC[79]*(1.00000 - STATES[34])+ STATES[31]*ALGEBRAIC[82]*STATES[34])+ ALGEBRAIC[92]*CONSTANTS[136]*ALGEBRAIC[90]*STATES[26]*( ALGEBRAIC[83]*(1.00000 - STATES[34])+ STATES[31]*ALGEBRAIC[84]*STATES[34]));
ALGEBRAIC[213] = (STATES[2] - STATES[9])/CONSTANTS[106];
ALGEBRAIC[215] = 1.00000/(1.00000+CONSTANTS[19]/ALGEBRAIC[49]);
ALGEBRAIC[217] =  CONSTANTS[110]*( (1.00000 - ALGEBRAIC[215])*STATES[43]+ ALGEBRAIC[215]*STATES[44]);
ALGEBRAIC[56] = 1.00000/(1.00000+( CONSTANTS[28]*CONSTANTS[29])/pow(CONSTANTS[29]+STATES[2], 2.00000)+( CONSTANTS[30]*CONSTANTS[31])/pow(CONSTANTS[31]+STATES[2], 2.00000));
ALGEBRAIC[103] = exp( - CONSTANTS[133]*4.00000*( pow(ALGEBRAIC[100], 1.0 / 2)/(1.00000+ pow(ALGEBRAIC[100], 1.0 / 2)) -  0.300000*ALGEBRAIC[100]));
ALGEBRAIC[108] = ( 4.00000*ALGEBRAIC[25]*( ALGEBRAIC[103]*STATES[9]*exp( 2.00000*ALGEBRAIC[28]) -  CONSTANTS[138]*CONSTANTS[2]))/(exp( 2.00000*ALGEBRAIC[28]) - 1.00000);
ALGEBRAIC[111] =  (1.00000 - CONSTANTS[53])*( (1.00000 - ALGEBRAIC[92])*CONSTANTS[120]*ALGEBRAIC[108]*STATES[26]*( ALGEBRAIC[79]*(1.00000 - STATES[35])+ STATES[31]*ALGEBRAIC[82]*STATES[35])+ ALGEBRAIC[92]*CONSTANTS[130]*ALGEBRAIC[108]*STATES[26]*( ALGEBRAIC[83]*(1.00000 - STATES[35])+ STATES[31]*ALGEBRAIC[84]*STATES[35]));
ALGEBRAIC[114] = ALGEBRAIC[93]+ALGEBRAIC[111];
ALGEBRAIC[115] = ALGEBRAIC[94]+ALGEBRAIC[112];
ALGEBRAIC[116] = ALGEBRAIC[97]+ALGEBRAIC[113];
ALGEBRAIC[211] = ( CONSTANTS[98]*STATES[9])/(CONSTANTS[99]+STATES[9]);
ALGEBRAIC[209] = ( CONSTANTS[97]*4.00000*ALGEBRAIC[25]*( ALGEBRAIC[103]*STATES[9]*exp( 2.00000*ALGEBRAIC[28]) -  CONSTANTS[138]*CONSTANTS[2]))/(exp( 2.00000*ALGEBRAIC[28]) - 1.00000);
ALGEBRAIC[66] =  (( CONSTANTS[5]*CONSTANTS[6])/( CONSTANTS[11]*CONSTANTS[7]))*log(CONSTANTS[4]/STATES[10]);
ALGEBRAIC[214] =  (( (1.00000 - CONSTANTS[103])*CONSTANTS[100])/(1.00000+CONSTANTS[102]/STATES[9]))*(STATES[0] - ALGEBRAIC[66]);
ALGEBRAIC[67] =  (( CONSTANTS[5]*CONSTANTS[6])/( CONSTANTS[11]*CONSTANTS[7]))*log(CONSTANTS[4]/STATES[11]);
ALGEBRAIC[212] =  (( CONSTANTS[103]*CONSTANTS[100])/(1.00000+CONSTANTS[102]/STATES[2]))*(STATES[0] - ALGEBRAIC[67]);
ALGEBRAIC[216] = ALGEBRAIC[212]+ALGEBRAIC[214];
ALGEBRAIC[218] =  CONSTANTS[101]*(STATES[0] - ALGEBRAIC[66]);
ALGEBRAIC[221] = (STATES[11] - STATES[10])/CONSTANTS[104];
ALGEBRAIC[219] = ( CONSTANTS[128]*0.00542500*STATES[9])/(STATES[9]+0.000920000);
ALGEBRAIC[220] = ( CONSTANTS[128]*2.75000*0.00542500*STATES[9])/((STATES[9]+0.000920000) - 0.000170000);
ALGEBRAIC[222] = 1.00000/(1.00000+CONSTANTS[19]/ALGEBRAIC[49]);
ALGEBRAIC[223] = ( 0.00488250*STATES[7])/15.0000;
ALGEBRAIC[224] =  CONSTANTS[111]*(( (1.00000 - ALGEBRAIC[222])*ALGEBRAIC[219]+ ALGEBRAIC[222]*ALGEBRAIC[220]) - ALGEBRAIC[223]);
ALGEBRAIC[53] = 1.00000/(1.00000+( CONSTANTS[113]*CONSTANTS[25])/pow(CONSTANTS[25]+STATES[9], 2.00000)+( CONSTANTS[26]*CONSTANTS[27])/pow(CONSTANTS[27]+STATES[9], 2.00000));
ALGEBRAIC[225] = (STATES[7] - STATES[8])/60.0000;
ALGEBRAIC[59] = 1.00000/(1.00000+( CONSTANTS[32]*CONSTANTS[33])/pow(CONSTANTS[33]+STATES[8], 2.00000));
}
