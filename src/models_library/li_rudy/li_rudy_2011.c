#include "li_rudy_2011.h"

// TODO: Move this to a function and remember that each cell must have this variables ...
// Do the same for the GPU code ...
//real ical, camkactive;
//real  qtr1,qtr2;
//real qrel1, qrel2, qup2;

GET_CELL_MODEL_DATA(init_cell_model_data) {

    if(get_initial_v)
        cell_model->initial_v = INITIAL_V;
    if(get_neq)
        cell_model->number_of_ode_equations = NEQ;
}

SET_ODE_INITIAL_CONDITIONS_CPU(set_model_initial_conditions_cpu) {

    log_info("Using Li & Rudy 2011 CPU model\n");

    uint32_t num_cells = solver->original_num_cells;
    solver->sv = (real*)malloc(NEQ*num_cells*sizeof(real));

    OMP(parallel for)
    for(uint32_t i = 0; i < num_cells; i++) {
        real *sv = &solver->sv[i * NEQ];
        sv[0] = -84.058830;  // V millivolt
        sv[1] = 0.000821;    // m dimensionless
        sv[2] = 0.995741;    // h dimensionless
        sv[3] = 0.999872;    // j dimensionless
        sv[4] = 0.000016;    // d dimensionless
        sv[5] = 0.999193;    // f dimensionless
        sv[6] = 0.988692;    // f2 dimensionless
        sv[7] = 0.965405;    // fca dimensionless
        sv[8] = 0.739378;    // fca2 dimensionless
        sv[9] = 0.001114;    // xs1 dimensionless
        sv[10] = 0.042234;   // xs2 dimensionless
        sv[11] = 0.069808;   // xr dimensionless
        sv[12] = 0.000119;   // a dimensionless
        sv[13] = 0.992541;   // i dimensionless
        sv[14] = 0.745628;   // i2 dimensionless
        sv[15] = 0.000329;   // ml dimensionless
        sv[16] = 0.046538;   // ml3 dimensionless
        sv[17] = 0.984170;   // hl dimensionless
        sv[18] = 0.853893;   // hl3 dimensionless
        sv[19] = 0.912569;   // jl dimensionless
        sv[20] = 0.827885;   // jl3 dimensionless
        sv[21] = 0.000135;   // casss dimensionless
        sv[22] = 1.510741;   // cajsr dimensionless
        sv[23] = 1.537577;   // cacsr dimensionless
        sv[24] = 1.538668;   // cansr dimensionless
        sv[25] = 0.000130;   // cassl dimensionless
        sv[26] = 11.501546;  // nai dimensionless
        sv[27] = 11.501230;  // nassl dimensionless
        sv[28] = 11.501240;  // nasss dimensionless
        sv[29] = 136.422946; // ki dimensionless
        sv[30] = 0.000053;   // cai millimolar
        sv[31] = 0.000437;   // b dimensionless
        sv[32] = 0.990384;   // g dimensionless
        sv[33] = 0.535627;   // u dimensionless
        sv[34] = 0.182859;   // y dimensionless
        sv[35] = 0.010600;   // camktrap dimensionless

        // Additional parameters
        sv[36] = 0.0; // ical
        sv[37] = 0.0; // camkactive
        sv[38] = 0.0; // qrel1
        sv[39] = 0.0; // qrel2
        sv[40] = 0.0; // qup2
    }

}

SOLVE_MODEL_ODES(solve_model_odes_cpu) {

    uint32_t sv_id;

    size_t num_cells_to_solve = ode_solver->num_cells_to_solve;
    uint32_t * cells_to_solve = ode_solver->cells_to_solve;
    real *sv = ode_solver->sv;
    real dt = ode_solver->min_dt;
    uint32_t num_steps = ode_solver->num_steps;

    OMP(parallel for private(sv_id))
    for (uint32_t i = 0; i < num_cells_to_solve; i++) {

        if(cells_to_solve)
            sv_id = cells_to_solve[i];
        else
            sv_id = i;

        for (int j = 0; j < num_steps; ++j) {
            solve_model_ode_cpu(dt, sv + (sv_id * NEQ), stim_currents[i]);
        }
    }
}

void solve_model_ode_cpu(real dt, real *sv, real stim_current)  {

    real rY[NEQ], rDY[NEQ];

    // Save odl value of the state vector
    for(int i = 0; i < NEQ; i++)
        rY[i] = sv[i];

    // Compute Right-hand-side of the ODE's
    RHS_cpu(rY, rDY, dt, stim_current);

    // Solve the ODE's using a mix between Forward Euler and Rush-Larsen
    for(int i = 0; i < NEQ; i++)
        sv[i] = rDY[i];
}

void RHS_cpu(const real *sv, real *rDY_, real dt, real stim_current) {

    //const double dtmin = 0.001;
    //const double dtmed = 0.005;
    //const double dtmax = 0.1;

    real v_old = sv[0];
    real m_old = sv[1];
    real h_old = sv[2];
    real j_old = sv[3];
    real d_old = sv[4];
    real f_old = sv[5];
    real f2_old = sv[6];
    real fca_old = sv[7];
    real fca2_old = sv[8];
    real xs1_old = sv[9];
    real xs2_old = sv[10];
    real xr_old = sv[11];
    real a_old = sv[12];
    real i_old = sv[13];
    real i2_old = sv[14];
    real ml_old = sv[15];
    real ml3_old = sv[16];
    real hl_old = sv[17];
    real hl3_old = sv[18];
    real jl_old = sv[19];
    real jl3_old = sv[20];
    real casss_old = sv[21];
    real cajsr_old = sv[22];
    real cacsr_old = sv[23];
    real cansr_old = sv[24];
    real cassl_old = sv[25];
    real nai_old = sv[26];
    real nassl_old = sv[27];
    real nasss_old = sv[28];
    real ki_old = sv[29];
    real cai_old = sv[30];
    real b_old = sv[31];
    real g_old = sv[32];
    real u_old = sv[33];
    real y_old = sv[34];
    real camktrap_old = sv[35];

    real ical       = sv[36];
    real camkactive = sv[37];
    real qrel1      = sv[38];
    real qrel2      = sv[39];
    real qup2       = sv[40];

    // Parameters
    //  CELL GEOMETRY
    const real pi = 3.14;
    const real radius = 0.00175;
    const real length = 0.0164;
    const real rcg = 1.54;

    const real vcell    = 1000*pi*radius*radius*length;
    const real ageo = 2*pi*radius*radius + 2*pi*radius*length;
    const real acap = rcg*ageo;
    const real vmyo = vcell * 0.60;
    const real vnsr = vcell * 0.04;
    //const real vmito   = vcell * 0.18;
    const real vjsr = vcell * 0.002;
    const real vcsr = vcell * 0.008;
    const real vsss = vcell * 0.02;
    const real vssl    = vcell * 0.15;

    // PHYSICAL CONSTANTS
    const real frdy = 96485;
    const real R = 8314;
    const real temp = 310;

    const real nao = 140;
    const real cao = 1.8;
    const real ko  = 5.4;
    //const real clo = 100;

    const real zna = 1;
    const real zk  = 1;
    //const real zcl = -1;
    const real zca = 2;
    //const real ganai = 0.75;
    //const real ganao = 0.75;
    //const real gaki  = 0.75;
    //const real gako  = 0.75;
    const real gacai = 1.0;
    const real gacao = 0.341;

    // CAMKII DYNAMICS
    const real camk0 = 0.05;
    const real alphacamk = 0.05;
    const real betacamk = 0.00068;
    const real kmcam = 0.0015;
    const real kmcamk = 0.15;
    //const real fca_dtaucamkbar = 10.0;

    // MEMBRANE IONIC CURRENTS
    const real gna = 18;
    const real gnal2 = 0.052;
    const real gnal3 = 0.018;
    const real pca = 1.9926e-4;
    //const real powtau = 10;
    const real gcat = 0.07875;
    const real gtos = 0.1414;
    const real gtof = 0.042;
    const real prnak = 0.014;
    //const real gnab = 0.0025;
    const real pcab = 3.99e-8;
    const real pnab = 0.64e-8;
    const real inacamax = 2.52;
    const real kmcaact = 0.000125;
    const real kmnai1 = 12.3;
    const real kmnao = 87.5;
    const real kmcai = 0.0036;
    const real kmcao = 1.3;
    const real nu = 0.35;
    const real ksat = 0.27;
    const real ibarnak = 1.1004;
    const real ipcabar = 0.0115;
    const real kmpca = 0.0005;

    // CALCIUM FLUXES AND CONCENTRATIONS
    const real IP3 = 0.0001;
    const real k1 = 150000;
    const real k1a = 16.5;
    const real k0 = 96000;
    const real k0a = 9.6;
    const real k2 = 1800;
    const real k2a = 0.21;
    const real tauip3r = 3.7;

    const real  dqupcamkbar = 0.75;
    const real  dkmplbbar = 0.00017;
    const real kmup   = 0.00028;
    const real nsrbar = 15.0;

    const real bsrbar = 0.019975;
    const real kmbsr = 0.00087;
    const real bslbar = 0.4777;
    const real kmbsl = 0.0087;

    const real csqnbar = 2.88;
    const real kmcsqn = 0.8;

    const real cmdnbar = 0.1125;
    const real kmcmdn = 2.38e-3;
    const real trpnbar = 3.15e-2;
    const real kmtrpn = 0.5e-3;
    const real trpnbar1 = 3.5e-3;
    const real cmdnbar1 = 1.25e-2;
    const real csqnbar1 = 1.2;

    // CALCIUM FLUXES RATE CONSTANTS
    const real tautr1 = 120;
    const real tautr2 = 120;
    const real gaptau = 12;
    const real sstau = 0.2;

    // comp_revs()
    real eca    = (R*temp/(zca*frdy))*log(cao/cassl_old);
    real ena    = (R*temp/frdy)*log(nao/nassl_old);
    real ek = (R*temp/frdy)*log(ko/ki_old);

    // comp_ina()
    real ma = 0.64*(v_old+37.13)/(1-exp(-0.1*(v_old+37.13)));
    real mb = 0.16*exp(-v_old/11);
    real ha, hb, ja, jb;
    if (v_old<-40)
    {
        ha = 0.135*exp((70+v_old)/-6.8);
        hb = 3.56*exp(0.079*v_old)+310000*exp(0.35*v_old);
        ja = (-127140*exp(0.2444*v_old)-0.003474*exp(-0.04391*v_old))*(v_old+37.78)/(1+exp(0.311*(v_old+79.23)));
        jb = 0.1212*exp(-0.01052*v_old)/(1+exp(-0.1378*(v_old+40.14)));
    }
    else
    {
        ha = 0.0;
        hb = 1/(0.13*(1+exp((v_old+10.66)/-11.1)));
        ja = 0.0;
        jb = 0.3*exp(-0.0000002535*v_old)/(1+exp(-0.1*(v_old+32)));
    }
    real mtau   = 1/(ma+mb);
    real htau   = 1/(ha + hb);
    real jtau   = 1/(ja+jb);
    real mss    = ma*mtau;
    real hss    = ha*htau;
    real jss    = 1*ja*jtau;

    // Rush-Larsen
    m_old   = mss-(mss-m_old)*exp(-dt/mtau);
    h_old   = hss-(hss-h_old)*exp(-dt/htau);
    j_old   = jss-(jss-j_old)*exp(-dt/jtau);

    real ina    = gna*pow(m_old,3)*h_old*j_old*(v_old-ena);

    // comp_inal()
    real mltau  = 1/(0.64*(v_old+37.13)/(1-exp(-0.1*(v_old+37.13))) + 0.16*exp(-v_old/11));
    real ml3tau  = mltau;
    real mlss   = 1/(1+exp(-(v_old+28)/7));
    real ml3ss   = 1/(1+exp(-(v_old+63)/7));
    real hltau   = 162+132/(1+exp(-(v_old+28)/5.5));
    real hl3tau  = 0.5*hltau;
    real hlss   = 1/(1+exp((v_old+28)/12));
    real hl3ss  = 1/(1+exp((v_old+63)/12));
    real jltau   = 411;
    real jl3tau  = 0.5*jltau;
    real jlss   = hlss;
    real jl3ss  = hl3ss;

    // Rush-Larsen
    ml_old      = mlss-(mlss-ml_old)*exp(-dt/mltau);
    ml3_old     = ml3ss-(ml3ss-ml3_old)*exp(-dt/ml3tau);
    hl_old      = hlss-(hlss-hl_old)*exp(-dt/hltau);
    hl3_old     = hl3ss-(hl3ss-hl3_old)*exp(-dt/hl3tau);
    jl_old      = jlss-(jlss-jl_old)*exp(-dt/jltau);
    jl3_old     = jl3ss-(jl3ss-jl3_old)*exp(-dt/jl3tau);

    real inal2   = gnal2*ml_old*hl_old*jl_old*(v_old-ena);
    real inal3   = gnal3*ml3_old*hl3_old*jl3_old*(v_old-ena);
    real inal    = inal2 + inal3;

    // comp_inab()
    real inab    = pnab*frdy*((frdy*v_old)/(R*temp))*(nassl_old*exp((frdy*v_old)/(R*temp)) - nao)/(exp((frdy*v_old)/(R*temp))-1);

    // comp_ical()
    real ibarca     = pca*zca*zca*(((v_old-15)*frdy*frdy)/(R*temp))*((gacai*casss_old*exp((zca*(v_old-15)*frdy)/(R*temp))-gacao*cao)/(exp((zca*(v_old-15)*frdy)/(R*temp))-1));
    real dss            = (1/(1.0+exp(-(v_old-2.0)/7.8)));
    real dtau       = (0.59+0.8*exp(0.052*(v_old+13))/(1+exp(0.132*(v_old+13))));
    real fss            = 1/(1.0 + exp((v_old+16.5)/9.5));
    real ftau        = 0.92/(0.125*exp(-(0.058*(v_old-2.5))*(0.045*(v_old-2.5)))+0.1);
    real f2ss        = fss;
    real f2tau       = 0.90/(0.02*exp(-(0.04*(v_old-18.6))*(0.045*(v_old-18.6)))+0.005);
    real fcass      = 0.3/(1 - ical/0.05) + 0.55/(1.0+casss_old/0.003)+0.15;
    real fcatau     = 10*camkactive/(camkactive+kmcam) + 0.5+1/(1.0+casss_old/0.003);
    real fca2ss     = 1.0/(1.0-ical/0.01);
    real fca2tau        = 1*(300.0/(1.0+exp((-ical-0.175)/0.04))+125.0);

    // Rush-Larsen
    d_old           = dss-(dss-d_old)*exp(-dt/dtau);
    f_old           = fss-(fss-f_old)*exp(-dt/ftau);
    f2_old          = f2ss-(f2ss-f2_old)*exp(-dt/f2tau);
    fca_old         = fcass-(fcass-fca_old)*exp(-dt/fcatau);
    fca2_old        = fca2ss-(fca2ss-fca2_old)*exp(-dt/fca2tau);

    ical        = d_old*f_old*f2_old*fca_old*fca2_old*ibarca;

    // comp_icat()
    real bss        = 1/(1+ exp (-(v_old+30)/7));
    real gss        = 1/(1+exp((v_old+61)/5));
    real taub   = 1/(1.068*exp((v_old+16.3)/30)+1.068*exp(-(v_old+16.3)/30));
    real taug    = 1/(0.015*exp(-(v_old+71.7)/83.3)+0.015*exp((v_old+71.7)/15.4));

    // Rush-Larsen
    b_old       = bss-(bss-b_old)*exp(-dt/taub);
    g_old       = gss-(gss-g_old)*exp(-dt/taug);

    real icat   = gcat*b_old*g_old*(v_old-eca);

    // comp_icab()
    real icab   = pcab*zca*zca*((v_old*frdy*frdy)/(R*temp))*((gacai*cassl_old*exp((zca*v_old*frdy)/(R*temp))-gacao*cao)/(exp((zca*v_old*frdy)/(R*temp))-1));

    // comp_itol()
    real atau   = 1/(25*exp((v_old-82)/18)/(1+exp((v_old-82)/18))+25*exp(-(v_old+52)/18)/(1+exp(-(v_old+52)/18)));
    real itau   = 2.86+ 1/(exp(-(v_old+125)/15)*0.1 + 0.1*exp((v_old+2)/26.5));
    real i2tau  = 21.5+ 1/(exp(-(v_old+138.2)/52)*0.005 + 0.003*exp((v_old+18)/12.5));
    real ass        = 1/(1+exp(-(v_old-8.9)/10.3));
    real iss        = 1/(1+exp((v_old+30)/11));
    real i2ss   = iss;

    // Rush-Larsen
    a_old       = ass-(ass-a_old)*exp(-dt/atau);
    i_old       = iss-(iss-i_old)*exp(-dt/itau);
    i2_old      = i2ss-(i2ss-i2_old)*exp(-dt/i2tau);

    real itos    = gtos*a_old*i_old*i2_old*(v_old-ek);
    real itof    = gtof*(v_old-ek)/(1+exp(-(v_old-3)/19.8));
    real ito1   = itos + itof;

    // comp_ikr()
    real gkr        = 0.0326*sqrt(ko/5.4);
    real xrss   = 1/(1+exp(-(v_old)/15));
    real xrtau   = 400.0/(1.0+exp(v_old/10.0)) + 100.0;
    real rkr        = 1/(1+exp((v_old)/35));

    // Rush-Larsen
    xr_old      = xrss-(xrss-xr_old)*exp(-dt/xrtau);

    real ikr        = gkr*xr_old*rkr*(v_old-ek);

    // comp_iks()
    real eks        = (R*temp/frdy)*log((ko+prnak*nao)/(ki_old+prnak*nassl_old));
    real gks        = 0.053*(1+0.6/(1+pow((0.000038/cassl_old),1.4)));
    real xsss   = 1/(1+exp(-(v_old-9)/13.7));
    real xs1tau = 200/(exp(-(v_old+10)/6) + exp((v_old-62)/55));
    real xs2tau = 1500+ 350/(exp(-(v_old+10)/4) + exp((v_old-90)/58));

    // Rush-Larsen
    xs1_old     = xsss-(xsss-xs1_old)*exp(-dt/xs1tau);
    xs2_old     = xsss-(xsss-xs2_old)*exp(-dt/xs2tau);

    real iks        = gks*xs1_old*xs2_old*(v_old-eks);

    // comp_ik1()
    real k1ss      = 1/(1+exp((v_old+103-(2.9+ko*2.175))/10.15));
    real gk1          = 0.12*sqrt(ko);
    real ik1          = gk1*k1ss*(v_old-ek);

    // comp_inaca()
    real allo       = 1/(1+pow((kmcaact/(1.5*casss_old)),2));
    real num            = inacamax*(pow(nasss_old,3)*cao*exp(nu*v_old*frdy/(R*temp))-pow(nao,3)*1.5*casss_old*exp((nu-1)*v_old*frdy/(R*temp)));
    real denommult  = 1+ksat*exp((nu-1)*v_old*frdy/(R*temp));
    real denomterm1 = kmcao*pow(nasss_old,3)+pow(kmnao,3)*1.5*casss_old+pow(kmnai1,3)*cao*(1+1.5*casss_old/kmcai);
    real denomterm2 = kmcai*pow(nao,3)*(1+pow(nasss_old/kmnai1,3))+pow(nasss_old,3)*cao+pow(nao,3)*1.5*casss_old;
    real deltaE     = num/(denommult*(denomterm1+denomterm2));
    real inacass        = 0.2*allo*deltaE;

    allo        = 1/(1+pow((kmcaact/(1.5*cassl_old)),2));
    num         = inacamax*(pow(nassl_old,3)*cao*exp(nu*v_old*frdy/(R*temp))-pow(nao,3)*1.5*cassl_old*exp((nu-1)*v_old*frdy/(R*temp)));
    denommult   = 1+ksat*exp((nu-1)*v_old*frdy/(R*temp));
    denomterm1  = kmcao*pow(nassl_old,3)+pow(kmnao,3)*1.5*cassl_old+pow(kmnai1,3)*cao*(1+1.5*cassl_old/kmcai);
    denomterm2  = kmcai*pow(nao,3)*(1+pow(nassl_old/kmnai1,3))+pow(nassl_old,3)*cao+pow(nao,3)*1.5*cassl_old;
    deltaE      = num/(denommult*(denomterm1+denomterm2));
    real inaca      = 0.8*allo*deltaE;

    // comp_inak()
    real inak   = ibarnak*(1/(1+exp(-1*(v_old+92)*frdy/(R*temp))))*pow((nassl_old/(nassl_old+2.6)),3)*(ko/(ko+0.8));

    // comp_ipca()
    real ipca   = ipcabar/((kmpca/cassl_old)+1);

    // comp_if()
    real yss       = 1/(1+exp((v_old+87)/9.5));
    real ytau      = 2000/(exp(-(v_old+132)/10) + exp((v_old+57)/60));

    // Rush-Larsen
    y_old         = yss - (yss-y_old)*exp(-dt/ytau);

    real ifna     = 0.012*y_old*y_old*(v_old-ena);
    real ifk       = 0.024*y_old*y_old*(v_old-ek);
    //real iftotal   = ifna + ifk;

    // comp_istim()
    real istim = stim_current;

    // comp_itot()
    real icatot = ical+icat+ipca+icab-2*inaca-2*inacass;
    real iktot  = ikr+iks+ik1-2*inak+ito1+ifk+1*istim;
    real inatot = 3*inak+ina+3*inaca+3*inacass+inal+ifna+inab;
    real itot   = icatot+iktot+inatot;

    // comp_ip3()
    // Forward Euler
    real du = dt*(casss_old*k2*(1-u_old) - k2a*u_old);
    u_old += du;
    real POip3 = tauip3r*IP3*casss_old*(1-u_old)/((1+IP3*k0/k0a)*(1+casss_old*k1/k1a));
    real qip3 = 10.920*(cajsr_old-casss_old)*(POip3);

    // comp_qrel1()
    real qdiff  = (casss_old-cassl_old)/sstau;
    real REL  = -((ical)*acap/(vsss*2.0*frdy) - (qrel1 + qip3)*vjsr/vsss + qdiff);
    real ireltau = 2*(1+1*(1/(1+pow((0.28/camkactive),8))))/(1+(0.0123/cajsr_old));
    real irelss;
    if (REL > 0)
        irelss  = 15*(1+1*(1/(1+pow((0.28/camkactive),8))))*REL/(1 + pow((1.0/cajsr_old),8));
    else
        irelss = 0;
    // Forward Euler
    qrel1 += dt*((irelss-qrel1)/ireltau);

    // comp_qrel2()
    real qgap  = (cassl_old-cai_old)/gaptau;
    REL  = (-qup2*vnsr/vmyo + qgap*vssl/vmyo+ (qrel2)*vcsr/vmyo);
    ireltau = 6*(1+1*(1/(1+pow((0.28/camkactive),8))))/(1+(0.0123/cacsr_old));
    if (REL > 0)
        irelss  = 91*(1+1*(1/(1+pow((0.28/camkactive),8))))*(REL)/(1 + pow((1/cacsr_old),8));
    else
        irelss = 0;
    // Forward Euler
    qrel2 += dt*((irelss-qrel2)/ireltau);

    // comp_qup1()
    real dkmplb     = dkmplbbar*camkactive/(kmcamk+camkactive);
    real dqupcamk   = dqupcamkbar*camkactive/(kmcamk+camkactive);
    real qup1       = 0.0002*(dqupcamk+1)/(1+pow((kmup-dkmplb)/cassl_old,1))-0.00105*cansr_old/nsrbar;

    dkmplb      = dkmplbbar*camkactive/(kmcamk+camkactive);
    dqupcamk    = dqupcamkbar*camkactive/(kmcamk+camkactive);
    qup2        = 0.0026*(dqupcamk+1)/(1+pow((kmup-dkmplb)/cai_old,1))-0.0042*cansr_old/nsrbar;

    // comp_qtr1()
    real qtr1       = (cansr_old-cajsr_old)/tautr1;

    // comp_qtr2()
    real qtr2       = (cansr_old-cacsr_old)/tautr2;

    // comp_conc()
    qdiff       = (casss_old-cassl_old)/sstau;
    qgap        = (cassl_old-cai_old)/gaptau;
    real qdiffna     = (nasss_old-nassl_old)/sstau;
    real qgapna      = (nassl_old-nai_old)/gaptau;


    // Forward Euler
    real dcasss     = dt*(-(ical-2*inacass)*acap/(vsss*2.0*frdy)+(qrel1+qip3)*vjsr/vsss-qdiff);
    real bsss        = 1/(1+(bsrbar*kmbsr/pow(kmbsr+casss_old,2))+(bslbar*kmbsl/pow(kmbsl+casss_old,2)));
    casss_old += bsss*dcasss;

    // Forward Euler
    real dcassl     = dt*(-(qup1)*vnsr/vssl+qdiff*vsss/vssl-qgap-(icat+ipca+icab-2*inaca)*acap/(vssl*2.0*frdy));
    real trpn        = trpnbar1*(cassl_old/(cassl_old+kmtrpn));
    real cmdn       = cmdnbar1*(cassl_old/(cassl_old+kmcmdn));
    real catotal        = trpn+cmdn+dcassl+cassl_old;
    real bmyo       = cmdnbar1+trpnbar1-catotal+kmtrpn+kmcmdn;
    real cmyo       = kmcmdn*kmtrpn-catotal*(kmtrpn+kmcmdn)+(trpnbar1*kmcmdn)+cmdnbar1*kmtrpn;
    real dmyo       = -kmtrpn*kmcmdn*catotal;
    cassl_old       = (2.0/3.0)*sqrt(bmyo*bmyo-3.0*cmyo)*cos(acos((9.0*bmyo*cmyo-2*bmyo*bmyo*bmyo-27*dmyo)/(2.0*pow((bmyo*bmyo-3.0*cmyo),1.5)))/3.0)-bmyo/3.0;

    real dcajsr     = dt*(qtr1-qrel1-qip3);
    real csqn1       = csqnbar1*(cajsr_old/(cajsr_old+kmcsqn));
    real bjsr        = csqnbar1 - csqn1-cajsr_old-dcajsr+kmcsqn;
    real cjsr        = kmcsqn*(csqn1+cajsr_old+dcajsr);
    cajsr_old       = (sqrt(bjsr*bjsr+4*cjsr)-bjsr)/2;

    real dcacsr     = dt*(qtr2-qrel2);
    real csqn        = csqnbar*(cacsr_old/(cacsr_old+kmcsqn));
    real bcsr        = csqnbar - csqn-cacsr_old-dcacsr+kmcsqn;
    real ccsr        = kmcsqn*(csqn+cacsr_old+dcacsr);
    cacsr_old       = (sqrt(bcsr*bcsr+4*ccsr)-bcsr)/2;

    // Forward Euler
    real dcansr     = dt*(qup1+qup2-qtr1*vjsr/vnsr-qtr2*vcsr/vnsr);
    cansr_old      += dcansr;

    // Forward Euler
    real dnasss     = dt*((-(3*inacass)*acap)/((vsss)*zna*frdy)-qdiffna);
    nasss_old      += dnasss;

    // Forward Euler
    real dnassl     = dt*((-(3*inak+ina+inal+3*inaca+ifna+inab)*acap)/((vssl)*zna*frdy)+qdiffna*vsss/vssl-qgapna);
    nassl_old      += dnassl;

    // Forward Euler
    real dnai        = dt*(qgapna*vssl/vmyo);
    nai_old        += dnai;

    // Forward Euler
    real dki            = dt*((-iktot*acap)/((vmyo+vssl+vsss)*zk*frdy));
    ki_old         += dki;

    // Forward Euler
    real dcai       = dt*(-(qup2)*vnsr/vmyo+qgap*vssl/vmyo+(qrel2)*vcsr/vmyo);
    trpn        = trpnbar*(cai_old/(cai_old+kmtrpn));
    cmdn        = cmdnbar*(cai_old/(cai_old+kmcmdn));
    catotal     = trpn+cmdn+dcai+cai_old;
    bmyo        = cmdnbar+trpnbar-catotal+kmtrpn+kmcmdn;
    cmyo        = kmcmdn*kmtrpn-catotal*(kmtrpn+kmcmdn)+(trpnbar*kmcmdn)+cmdnbar*kmtrpn;
    dmyo        = -kmtrpn*kmcmdn*catotal;
    cai_old         = (2.0/3.0)*sqrt(bmyo*bmyo-3.0*cmyo)*cos(acos((9.0*bmyo*cmyo-2*bmyo*bmyo*bmyo-27*dmyo)/(2.0*pow((bmyo*bmyo-3.0*cmyo),1.5)))/3.0)-bmyo/3.0;

    // Unused ...
    //real caavg       = (casss_old*vsss+cassl*vssl+cai_old*vmyo)/(vsss+vmyo+vssl);

    real camkbound  = camk0*(1-camktrap_old)*1/(1+(kmcam/casss_old));

    // Forward Euler
    camktrap_old    = dt*(alphacamk*camkbound*(camkbound+camktrap_old)-betacamk*camktrap_old) + camktrap_old;
    camkactive  = camkbound+camktrap_old;

    real dvdt = -itot;
    v_old   += dt*dvdt;

    // Rush-Larsen
    rDY_[1] = m_old;
    rDY_[2] = h_old;
    rDY_[3] = j_old;
    rDY_[4] = d_old;
    rDY_[5] = f_old;
    rDY_[6] = f2_old;
    rDY_[7] = fca_old;
    rDY_[8] = fca2_old;
    rDY_[9] = xs1_old;
    rDY_[10] = xs2_old;
    rDY_[11] = xr_old;
    rDY_[12] = a_old;
    rDY_[13] = i_old;
    rDY_[14] = i2_old;
    rDY_[15] = ml_old;
    rDY_[16] = ml3_old;
    rDY_[17] = hl_old;
    rDY_[18] = hl3_old;
    rDY_[19] = jl_old;
    rDY_[20] = jl3_old;
    rDY_[31] = b_old;
    rDY_[32] = g_old;
    rDY_[34] = y_old;

    // Forward Euler (I already calculated the Forward Euler scheme here ...)
    rDY_[0] = v_old;
    rDY_[21] = casss_old;
    rDY_[22] = cajsr_old;
    rDY_[23] = cacsr_old;
    rDY_[24] = cansr_old;
    rDY_[25] = cassl_old;
    rDY_[26] = nai_old;
    rDY_[27] = nassl_old;
    rDY_[28] = nasss_old;
    rDY_[29] = ki_old;
    rDY_[30] = cai_old;
    rDY_[33] = u_old;
    rDY_[35] = camktrap_old;

    rDY_[36]      = ical;
    rDY_[37]      = camkactive;
    rDY_[38]      = qrel1;
    rDY_[39]      = qrel2;
    rDY_[40]      = qup2;


}
