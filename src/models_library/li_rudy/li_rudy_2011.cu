#include "../../gpu_utils/gpu_utils.h"
#include <stddef.h>
#include <stdint.h>

#include "li_rudy_2011.h"

extern "C" SET_ODE_INITIAL_CONDITIONS_GPU(set_model_initial_conditions_gpu) {

    log_info("Using Li & Rudy 2011 GPU model\n");

    uint32_t num_volumes = solver->original_num_cells;

    // execution configuration
    const int GRID  = (num_volumes + BLOCK_SIZE - 1)/BLOCK_SIZE;

    size_t size = num_volumes*sizeof(real);

    check_cuda_error(cudaMallocPitch((void **) &(solver->sv), &pitch_h, size, (size_t )NEQ));
    check_cuda_error(cudaMemcpyToSymbol(pitch, &pitch_h, sizeof(size_t)));


    kernel_set_model_inital_conditions <<<GRID, BLOCK_SIZE>>>(solver->sv, num_volumes);

    check_cuda_error( cudaPeekAtLastError() );
    cudaDeviceSynchronize();
    return pitch_h;

}

extern "C" SOLVE_MODEL_ODES(solve_model_odes_gpu) {

    size_t num_cells_to_solve = ode_solver->num_cells_to_solve;
    uint32_t * cells_to_solve = ode_solver->cells_to_solve;
    real *sv = ode_solver->sv;
    real dt = ode_solver->min_dt;
    uint32_t num_steps = ode_solver->num_steps;

    // execution configuration
    const int GRID  = ((int)num_cells_to_solve + BLOCK_SIZE - 1)/BLOCK_SIZE;

    size_t stim_currents_size = sizeof(real)*num_cells_to_solve;
    size_t cells_to_solve_size = sizeof(uint32_t)*num_cells_to_solve;

    real *stims_currents_device;
    check_cuda_error(cudaMalloc((void **) &stims_currents_device, stim_currents_size));
    check_cuda_error(cudaMemcpy(stims_currents_device, stim_currents, stim_currents_size, cudaMemcpyHostToDevice));

    //the array cells to solve is passed when we are using and adaptive mesh
    uint32_t *cells_to_solve_device = NULL;
    if(cells_to_solve != NULL) {
        check_cuda_error(cudaMalloc((void **) &cells_to_solve_device, cells_to_solve_size));
        check_cuda_error(cudaMemcpy(cells_to_solve_device, cells_to_solve, cells_to_solve_size, cudaMemcpyHostToDevice));
    }
    solve_gpu <<<GRID, BLOCK_SIZE>>>(dt, sv, stims_currents_device, cells_to_solve_device, num_cells_to_solve, num_steps);

    check_cuda_error( cudaPeekAtLastError() );

    check_cuda_error(cudaFree(stims_currents_device));
    if(cells_to_solve_device) check_cuda_error(cudaFree(cells_to_solve_device));

}

__global__ void kernel_set_model_inital_conditions(real *sv, int num_volumes) {
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;

    if (threadID < num_volumes) {

        *((real * )((char *) sv + pitch * 0) + threadID)		= -84.058830;       // V millivolt
        *((real * )((char *) sv + pitch * 1) + threadID)		= 0.000821;         // m dimensionless
        *((real * )((char *) sv + pitch * 2) + threadID)		= 0.995741;         // h dimensionless
        *((real * )((char *) sv + pitch * 3) + threadID)		= 0.999872;         // j dimensionless
        *((real * )((char *) sv + pitch * 4) + threadID)		= 0.000016;         // d dimensionless
        *((real * )((char *) sv + pitch * 5) + threadID)		= 0.999193;         // f dimensionless
        *((real * )((char *) sv + pitch * 6) + threadID)		= 0.988692;         // f2 dimensionless
        *((real * )((char *) sv + pitch * 7) + threadID)		= 0.965405;         // fca dimensionless
        *((real * )((char *) sv + pitch * 8) + threadID)	    = 0.739378;         // fca2 dimensionless
        *((real * )((char *) sv + pitch * 9) + threadID)		= 0.001114;         // xs1 dimensionless
        *((real * )((char *) sv + pitch * 10) + threadID)		= 0.042234;         // xs2 dimensionless
        *((real * )((char *) sv + pitch * 11) + threadID)		= 0.069808;         // xr dimensionless
        *((real * )((char *) sv + pitch * 12) + threadID)		= 0.000119;         // a dimensionless
        *((real * )((char *) sv + pitch * 13) + threadID)		= 0.992541;         // i dimensionless
        *((real * )((char *) sv + pitch * 14) + threadID)		= 0.745628;         // i2 dimensionless
        *((real * )((char *) sv + pitch * 15) + threadID)		= 0.000329;         // ml dimensionless
        *((real * )((char *) sv + pitch * 16) + threadID)		= 0.046538;         // ml3 dimensionless
        *((real * )((char *) sv + pitch * 17) + threadID)		= 0.984170;         // hl dimensionless
        *((real * )((char *) sv + pitch * 18) + threadID)		= 0.853893;         // hl3 dimensionless
        *((real * )((char *) sv + pitch * 19) + threadID)		= 0.912569;         // jl dimensionless
        *((real * )((char *) sv + pitch * 20) + threadID)		= 0.827885;         // jl3 dimensionless
        *((real * )((char *) sv + pitch * 21) + threadID)	    = 0.000135;         // casss dimensionless
        *((real * )((char *) sv + pitch * 22) + threadID)	    = 1.510741;         // cajsr dimensionless
        *((real * )((char *) sv + pitch * 23) + threadID)	    = 1.537577;         // cacsr dimensionless
        *((real * )((char *) sv + pitch * 24) + threadID)	    = 1.538668;         // cansr dimensionless
        *((real * )((char *) sv + pitch * 25) + threadID)	    = 0.000130;         // cassl dimensionless
        *((real * )((char *) sv + pitch * 26) + threadID)	    = 11.501546;        // nai dimensionless
        *((real * )((char *) sv + pitch * 27) + threadID)	    = 11.501230;        // nassl dimensionless
        *((real * )((char *) sv + pitch * 28) + threadID)	    = 11.501240;        // nasss dimensionless
        *((real * )((char *) sv + pitch * 29) + threadID)		= 136.422946;       // ki dimensionless
        *((real * )((char *) sv + pitch * 30) + threadID)		= 0.000053;         // cai millimolar
        *((real * )((char *) sv + pitch * 31) + threadID)	    = 0.000437;         // b dimensionless
        *((real * )((char *) sv + pitch * 32) + threadID)	    = 0.990384;         // g dimensionless
        *((real * )((char *) sv + pitch * 33) + threadID)       = 0.535627;        // u dimensionless
        *((real * )((char *) sv + pitch * 34) + threadID)       = 0.182859;        // y dimensionless
        *((real * )((char *) sv + pitch * 35) + threadID)       = 0.010600;        // camktrap dimensionless

        // Additional parameters
        *((real * )((char *) sv + pitch * 36) + threadID)      = 0.0;              // ical
        *((real * )((char *) sv + pitch * 37) + threadID)      = 0.0;              // camkactive
        *((real * )((char *) sv + pitch * 38) + threadID)      = 0.0;              // qrel1
        *((real * )((char *) sv + pitch * 39) + threadID)      = 0.0;              // qrel2
        *((real * )((char *) sv + pitch * 40) + threadID)      = 0.0;              // qup2

    }
}

// Solving the model for each cell in the tissue matrix ni x nj
__global__ void solve_gpu(real dt, real *sv, real* stim_currents,
                          uint32_t *cells_to_solve, uint32_t num_cells_to_solve,
                          int num_steps)
{
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;
    int sv_id;

    // Each thread solves one cell model
    if(threadID < num_cells_to_solve) {
        if(cells_to_solve)
            sv_id = cells_to_solve[threadID];
        else
            sv_id = threadID;

        real rDY[NEQ];

        for (int n = 0; n < num_steps; ++n) {

            // Compute Right-hand-side of the ODE's
            RHS_gpu(sv, rDY, stim_currents[threadID], dt, sv_id);

            // Solve the ODE's using a mix between Forward Euler and Rush-Larsen
            for(int i = 0; i < NEQ; i++) {
                *((real*)((char*)sv + pitch * i) + sv_id) = rDY[i];
            }            

        }

    }
}

inline __device__ void RHS_gpu(real *sv_, real *rDY_, real stim_current, real dt, int threadID_) {

    //const double dtmin = 0.001;		     
    //const double dtmed = 0.005;		
    //const double dtmax = 0.1;

    real v_old = *((real*)((char*)sv_ + pitch * 0) + threadID_);
    real m_old = *((real*)((char*)sv_ + pitch * 1) + threadID_);
    real h_old = *((real*)((char*)sv_ + pitch * 2) + threadID_);
    real j_old = *((real*)((char*)sv_ + pitch * 3) + threadID_);
    real d_old = *((real*)((char*)sv_ + pitch * 4) + threadID_);
    real f_old = *((real*)((char*)sv_ + pitch * 5) + threadID_);
    real f2_old = *((real*)((char*)sv_ + pitch * 6) + threadID_);
    real fca_old = *((real*)((char*)sv_ + pitch * 7) + threadID_);
    real fca2_old = *((real*)((char*)sv_ + pitch * 8) + threadID_);
    real xs1_old = *((real*)((char*)sv_ + pitch * 9) + threadID_);
    real xs2_old = *((real*)((char*)sv_ + pitch * 10) + threadID_);
    real xr_old = *((real*)((char*)sv_ + pitch * 11) + threadID_);
    real a_old = *((real*)((char*)sv_ + pitch * 12) + threadID_);
    real i_old = *((real*)((char*)sv_ + pitch * 13) + threadID_);
    real i2_old = *((real*)((char*)sv_ + pitch * 14) + threadID_);
    real ml_old = *((real*)((char*)sv_ + pitch * 15) + threadID_);
    real ml3_old = *((real*)((char*)sv_ + pitch * 16) + threadID_);
    real hl_old = *((real*)((char*)sv_ + pitch * 17) + threadID_);
    real hl3_old = *((real*)((char*)sv_ + pitch * 18) + threadID_);
    real jl_old = *((real*)((char*)sv_ + pitch * 19) + threadID_);
    real jl3_old = *((real*)((char*)sv_ + pitch * 20) + threadID_);
    real casss_old = *((real*)((char*)sv_ + pitch * 21) + threadID_);
    real cajsr_old = *((real*)((char*)sv_ + pitch * 22) + threadID_);
    real cacsr_old = *((real*)((char*)sv_ + pitch * 23) + threadID_);
    real cansr_old = *((real*)((char*)sv_ + pitch * 24) + threadID_);
    real cassl_old = *((real*)((char*)sv_ + pitch * 25) + threadID_);
    real nai_old = *((real*)((char*)sv_ + pitch * 26) + threadID_);
    real nassl_old = *((real*)((char*)sv_ + pitch * 27) + threadID_);
    real nasss_old = *((real*)((char*)sv_ + pitch * 28) + threadID_);
    real ki_old = *((real*)((char*)sv_ + pitch * 29) + threadID_);
    real cai_old = *((real*)((char*)sv_ + pitch * 30) + threadID_);
    real b_old = *((real*)((char*)sv_ + pitch * 31) + threadID_);
    real g_old = *((real*)((char*)sv_ + pitch * 32) + threadID_);
    real u_old = *((real*)((char*)sv_ + pitch * 33) + threadID_);
    real y_old = *((real*)((char*)sv_ + pitch * 34) + threadID_);
    real camktrap_old = *((real*)((char*)sv_ + pitch * 35) + threadID_);

    real ical       = *((real*)((char*)sv_ + pitch * 36) + threadID_);
    real camkactive = *((real*)((char*)sv_ + pitch * 37) + threadID_);
    real qrel1      = *((real*)((char*)sv_ + pitch * 38) + threadID_);      
    real qrel2      = *((real*)((char*)sv_ + pitch * 39) + threadID_);
    real qup2       = *((real*)((char*)sv_ + pitch * 40) + threadID_);

    // Parameters
    //  CELL GEOMETRY
    const real pi = 3.14;			
    const real radius = 0.00175;	
    const real length = 0.0164;	
    const real rcg = 1.54;

	const real vcell	= 1000*pi*radius*radius*length;
	const real ageo	= 2*pi*radius*radius + 2*pi*radius*length;
	const real acap	= rcg*ageo;
	const real vmyo	= vcell * 0.60;
	const real vnsr	= vcell * 0.04;
	//const real vmito   = vcell * 0.18;
	const real vjsr	= vcell * 0.002;
	const real vcsr	= vcell * 0.008;
	const real vsss	= vcell * 0.02;
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
    real eca	= (R*temp/(zca*frdy))*log(cao/cassl_old);
	real ena	= (R*temp/frdy)*log(nao/nassl_old);
	real ek	= (R*temp/frdy)*log(ko/ki_old);

    // comp_ina()
    real ma	= 0.64*(v_old+37.13)/(1-exp(-0.1*(v_old+37.13)));
	real mb	= 0.16*exp(-v_old/11);
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
	real mtau	= 1/(ma+mb);
	real htau	= 1/(ha + hb);
	real jtau	= 1/(ja+jb);
	real mss	= ma*mtau;
	real hss	= ha*htau;
	real jss	= 1*ja*jtau;
	
    // Rush-Larsen
    m_old	= mss-(mss-m_old)*exp(-dt/mtau);
	h_old	= hss-(hss-h_old)*exp(-dt/htau);
	j_old	= jss-(jss-j_old)*exp(-dt/jtau);

	real ina	= gna*pow(m_old,3)*h_old*j_old*(v_old-ena);

    // comp_inal()
    real mltau	= 1/(0.64*(v_old+37.13)/(1-exp(-0.1*(v_old+37.13))) + 0.16*exp(-v_old/11));
	real ml3tau  = mltau;
	real mlss	= 1/(1+exp(-(v_old+28)/7));
	real ml3ss   = 1/(1+exp(-(v_old+63)/7));
	real hltau   = 162+132/(1+exp(-(v_old+28)/5.5));
	real hl3tau  = 0.5*hltau;
	real hlss	= 1/(1+exp((v_old+28)/12));
	real hl3ss	= 1/(1+exp((v_old+63)/12));
	real jltau   = 411;
	real jl3tau  = 0.5*jltau;
	real jlss	= hlss;
	real jl3ss	= hl3ss;
	
    // Rush-Larsen
    ml_old	    = mlss-(mlss-ml_old)*exp(-dt/mltau);
	ml3_old     = ml3ss-(ml3ss-ml3_old)*exp(-dt/ml3tau);
	hl_old	    = hlss-(hlss-hl_old)*exp(-dt/hltau);
	hl3_old     = hl3ss-(hl3ss-hl3_old)*exp(-dt/hl3tau);
	jl_old	    = jlss-(jlss-jl_old)*exp(-dt/jltau);
	jl3_old     = jl3ss-(jl3ss-jl3_old)*exp(-dt/jl3tau);
	
    real inal2   = gnal2*ml_old*hl_old*jl_old*(v_old-ena);
	real inal3   = gnal3*ml3_old*hl3_old*jl3_old*(v_old-ena);
	real inal    = inal2 + inal3;

    // comp_inab()
    real inab    = pnab*frdy*((frdy*v_old)/(R*temp))*(nassl_old*exp((frdy*v_old)/(R*temp)) - nao)/(exp((frdy*v_old)/(R*temp))-1);

    // comp_ical()
    real ibarca		= pca*zca*zca*(((v_old-15)*frdy*frdy)/(R*temp))*((gacai*casss_old*exp((zca*(v_old-15)*frdy)/(R*temp))-gacao*cao)/(exp((zca*(v_old-15)*frdy)/(R*temp))-1));
	real dss		    = (1/(1.0+exp(-(v_old-2.0)/7.8)));
	real dtau		= (0.59+0.8*exp(0.052*(v_old+13))/(1+exp(0.132*(v_old+13))));
	real fss	        = 1/(1.0 + exp((v_old+16.5)/9.5));
	real ftau        = 0.92/(0.125*exp(-(0.058*(v_old-2.5))*(0.045*(v_old-2.5)))+0.1);
	real f2ss        = fss;
	real f2tau       = 0.90/(0.02*exp(-(0.04*(v_old-18.6))*(0.045*(v_old-18.6)))+0.005);
	real fcass		= 0.3/(1 - ical/0.05) + 0.55/(1.0+casss_old/0.003)+0.15;
	real fcatau		= 10*camkactive/(camkactive+kmcam) + 0.5+1/(1.0+casss_old/0.003);
	real fca2ss		= 1.0/(1.0-ical/0.01);
	real fca2tau		= 1*(300.0/(1.0+exp((-ical-0.175)/0.04))+125.0);

    // Rush-Larsen
    d_old		    = dss-(dss-d_old)*exp(-dt/dtau);
	f_old		    = fss-(fss-f_old)*exp(-dt/ftau);
	f2_old		    = f2ss-(f2ss-f2_old)*exp(-dt/f2tau);
	fca_old		    = fcass-(fcass-fca_old)*exp(-dt/fcatau);
	fca2_old		= fca2ss-(fca2ss-fca2_old)*exp(-dt/fca2tau);

	ical		= d_old*f_old*f2_old*fca_old*fca2_old*ibarca;

    // comp_icat()
    real bss	    = 1/(1+ exp (-(v_old+30)/7));
	real gss	    = 1/(1+exp((v_old+61)/5));
	real taub	= 1/(1.068*exp((v_old+16.3)/30)+1.068*exp(-(v_old+16.3)/30));
	real taug    = 1/(0.015*exp(-(v_old+71.7)/83.3)+0.015*exp((v_old+71.7)/15.4));
	
    // Rush-Larsen
    b_old	    = bss-(bss-b_old)*exp(-dt/taub);
	g_old	    = gss-(gss-g_old)*exp(-dt/taug);
	
    real icat	= gcat*b_old*g_old*(v_old-eca);

    // comp_icab()
    real icab	= pcab*zca*zca*((v_old*frdy*frdy)/(R*temp))*((gacai*cassl_old*exp((zca*v_old*frdy)/(R*temp))-gacao*cao)/(exp((zca*v_old*frdy)/(R*temp))-1));

    // comp_itol()
    real atau	= 1/(25*exp((v_old-82)/18)/(1+exp((v_old-82)/18))+25*exp(-(v_old+52)/18)/(1+exp(-(v_old+52)/18)));
	real itau	= 2.86+ 1/(exp(-(v_old+125)/15)*0.1 + 0.1*exp((v_old+2)/26.5));
	real i2tau	= 21.5+ 1/(exp(-(v_old+138.2)/52)*0.005 + 0.003*exp((v_old+18)/12.5));
	real ass	    = 1/(1+exp(-(v_old-8.9)/10.3));
	real iss	    = 1/(1+exp((v_old+30)/11));
	real i2ss	= iss;
	
    // Rush-Larsen
    a_old	    = ass-(ass-a_old)*exp(-dt/atau);
	i_old	    = iss-(iss-i_old)*exp(-dt/itau);
	i2_old	    = i2ss-(i2ss-i2_old)*exp(-dt/i2tau);

	real itos    = gtos*a_old*i_old*i2_old*(v_old-ek);
	real itof    = gtof*(v_old-ek)/(1+exp(-(v_old-3)/19.8));
	real ito1	= itos + itof;

    // comp_ikr()
    real gkr	    = 0.0326*sqrt(ko/5.4);
	real xrss	= 1/(1+exp(-(v_old)/15));
	real xrtau   = 400.0/(1.0+exp(v_old/10.0)) + 100.0;
	real rkr	    = 1/(1+exp((v_old)/35));
	
    // Rush-Larsen
    xr_old	    = xrss-(xrss-xr_old)*exp(-dt/xrtau);
	
    real ikr	    = gkr*xr_old*rkr*(v_old-ek);

    // comp_iks()
    real eks	    = (R*temp/frdy)*log((ko+prnak*nao)/(ki_old+prnak*nassl_old));
	real gks	    = 0.053*(1+0.6/(1+pow((0.000038/cassl_old),1.4)));
	real xsss	= 1/(1+exp(-(v_old-9)/13.7));
	real xs1tau	= 200/(exp(-(v_old+10)/6) + exp((v_old-62)/55));
	real xs2tau	= 1500+ 350/(exp(-(v_old+10)/4) + exp((v_old-90)/58));

    // Rush-Larsen
	xs1_old	    = xsss-(xsss-xs1_old)*exp(-dt/xs1tau);
	xs2_old	    = xsss-(xsss-xs2_old)*exp(-dt/xs2tau);
	
    real iks	    = gks*xs1_old*xs2_old*(v_old-eks);

    // comp_ik1()
    real k1ss      = 1/(1+exp((v_old+103-(2.9+ko*2.175))/10.15));
	real gk1	      = 0.12*sqrt(ko);
	real ik1	      = gk1*k1ss*(v_old-ek);

    // comp_inaca()
    real allo		= 1/(1+pow((kmcaact/(1.5*casss_old)),2));
	real num		    = inacamax*(pow(nasss_old,3)*cao*exp(nu*v_old*frdy/(R*temp))-pow(nao,3)*1.5*casss_old*exp((nu-1)*v_old*frdy/(R*temp)));
	real denommult	= 1+ksat*exp((nu-1)*v_old*frdy/(R*temp));
	real denomterm1	= kmcao*pow(nasss_old,3)+pow(kmnao,3)*1.5*casss_old+pow(kmnai1,3)*cao*(1+1.5*casss_old/kmcai);
	real denomterm2	= kmcai*pow(nao,3)*(1+pow(nasss_old/kmnai1,3))+pow(nasss_old,3)*cao+pow(nao,3)*1.5*casss_old;
	real deltaE		= num/(denommult*(denomterm1+denomterm2));
	real inacass		= 0.2*allo*deltaE;
	
	allo		= 1/(1+pow((kmcaact/(1.5*cassl_old)),2));
	num		    = inacamax*(pow(nassl_old,3)*cao*exp(nu*v_old*frdy/(R*temp))-pow(nao,3)*1.5*cassl_old*exp((nu-1)*v_old*frdy/(R*temp)));
	denommult	= 1+ksat*exp((nu-1)*v_old*frdy/(R*temp));
	denomterm1	= kmcao*pow(nassl_old,3)+pow(kmnao,3)*1.5*cassl_old+pow(kmnai1,3)*cao*(1+1.5*cassl_old/kmcai);
	denomterm2	= kmcai*pow(nao,3)*(1+pow(nassl_old/kmnai1,3))+pow(nassl_old,3)*cao+pow(nao,3)*1.5*cassl_old;
	deltaE		= num/(denommult*(denomterm1+denomterm2));
	real inaca		= 0.8*allo*deltaE;

    // comp_inak()
    real inak	= ibarnak*(1/(1+exp(-1*(v_old+92)*frdy/(R*temp))))*pow((nassl_old/(nassl_old+2.6)),3)*(ko/(ko+0.8));

    // comp_ipca()
    real ipca	= ipcabar/((kmpca/cassl_old)+1);

    // comp_if()
    real yss       = 1/(1+exp((v_old+87)/9.5));
	real ytau      = 2000/(exp(-(v_old+132)/10) + exp((v_old+57)/60));
	
    // Rush-Larsen
    y_old         = yss - (yss-y_old)*exp(-dt/ytau);

	real ifna	  = 0.012*y_old*y_old*(v_old-ena);
	real ifk       = 0.024*y_old*y_old*(v_old-ek);
	//real iftotal   = ifna + ifk;

    // comp_istim()
    real istim = stim_current;
    
    // comp_itot()
    real icatot	= ical+icat+ipca+icab-2*inaca-2*inacass;
    real iktot	= ikr+iks+ik1-2*inak+ito1+ifk+1*istim;
    real inatot	= 3*inak+ina+3*inaca+3*inacass+inal+ifna+inab;
    real itot	= icatot+iktot+inatot;

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
    real dkmplb		= dkmplbbar*camkactive/(kmcamk+camkactive);
	real dqupcamk	= dqupcamkbar*camkactive/(kmcamk+camkactive); 
	real qup1		= 0.0002*(dqupcamk+1)/(1+pow((kmup-dkmplb)/cassl_old,1))-0.00105*cansr_old/nsrbar;

    dkmplb		= dkmplbbar*camkactive/(kmcamk+camkactive);
	dqupcamk	= dqupcamkbar*camkactive/(kmcamk+camkactive); 
	qup2		= 0.0026*(dqupcamk+1)/(1+pow((kmup-dkmplb)/cai_old,1))-0.0042*cansr_old/nsrbar;

    // comp_qtr1()
    real qtr1		= (cansr_old-cajsr_old)/tautr1;

    // comp_qtr2()
    real qtr2		= (cansr_old-cacsr_old)/tautr2;

    // comp_conc()
    qdiff       = (casss_old-cassl_old)/sstau;  
	qgap        = (cassl_old-cai_old)/gaptau;  
    real qdiffna     = (nasss_old-nassl_old)/sstau;
    real qgapna      = (nassl_old-nai_old)/gaptau;
    
    
    // Forward Euler
    real dcasss		= dt*(-(ical-2*inacass)*acap/(vsss*2.0*frdy)+(qrel1+qip3)*vjsr/vsss-qdiff);
    real bsss        = 1/(1+(bsrbar*kmbsr/pow(kmbsr+casss_old,2))+(bslbar*kmbsl/pow(kmbsl+casss_old,2)));
	casss_old += bsss*dcasss;
	
    // Forward Euler
	real dcassl		= dt*(-(qup1)*vnsr/vssl+qdiff*vsss/vssl-qgap-(icat+ipca+icab-2*inaca)*acap/(vssl*2.0*frdy));
	real trpn        = trpnbar1*(cassl_old/(cassl_old+kmtrpn));
	real cmdn		= cmdnbar1*(cassl_old/(cassl_old+kmcmdn));
	real catotal		= trpn+cmdn+dcassl+cassl_old;
	real bmyo		= cmdnbar1+trpnbar1-catotal+kmtrpn+kmcmdn;
	real cmyo		= kmcmdn*kmtrpn-catotal*(kmtrpn+kmcmdn)+(trpnbar1*kmcmdn)+cmdnbar1*kmtrpn;
	real dmyo		= -kmtrpn*kmcmdn*catotal;
	cassl_old		= (2.0/3.0)*sqrt(bmyo*bmyo-3.0*cmyo)*cos(acos((9.0*bmyo*cmyo-2*bmyo*bmyo*bmyo-27*dmyo)/(2.0*pow((bmyo*bmyo-3.0*cmyo),1.5)))/3.0)-bmyo/3.0;   
 	
	real dcajsr		= dt*(qtr1-qrel1-qip3);
	real csqn1       = csqnbar1*(cajsr_old/(cajsr_old+kmcsqn));
	real bjsr        = csqnbar1 - csqn1-cajsr_old-dcajsr+kmcsqn;
	real cjsr        = kmcsqn*(csqn1+cajsr_old+dcajsr);
	cajsr_old       = (sqrt(bjsr*bjsr+4*cjsr)-bjsr)/2;
	
	real dcacsr		= dt*(qtr2-qrel2);
	real csqn        = csqnbar*(cacsr_old/(cacsr_old+kmcsqn));
	real bcsr        = csqnbar - csqn-cacsr_old-dcacsr+kmcsqn;
	real ccsr        = kmcsqn*(csqn+cacsr_old+dcacsr);
	cacsr_old       = (sqrt(bcsr*bcsr+4*ccsr)-bcsr)/2;
	
    // Forward Euler
	real dcansr	    = dt*(qup1+qup2-qtr1*vjsr/vnsr-qtr2*vcsr/vnsr);
 	cansr_old	   += dcansr;
 	
    // Forward Euler 
	real dnasss	    = dt*((-(3*inacass)*acap)/((vsss)*zna*frdy)-qdiffna); 
	nasss_old      += dnasss;
	
    // Forward Euler
	real dnassl	    = dt*((-(3*inak+ina+inal+3*inaca+ifna+inab)*acap)/((vssl)*zna*frdy)+qdiffna*vsss/vssl-qgapna);
	nassl_old	   += dnassl;
	
    // Forward Euler
	real dnai        = dt*(qgapna*vssl/vmyo);
	nai_old        += dnai;
	
    // Forward Euler
	real dki	        = dt*((-iktot*acap)/((vmyo+vssl+vsss)*zk*frdy));
	ki_old         += dki;
	
    // Forward Euler
	real dcai		= dt*(-(qup2)*vnsr/vmyo+qgap*vssl/vmyo+(qrel2)*vcsr/vmyo);
	trpn        = trpnbar*(cai_old/(cai_old+kmtrpn));
	cmdn		= cmdnbar*(cai_old/(cai_old+kmcmdn));
	catotal		= trpn+cmdn+dcai+cai_old;
	bmyo		= cmdnbar+trpnbar-catotal+kmtrpn+kmcmdn;
	cmyo		= kmcmdn*kmtrpn-catotal*(kmtrpn+kmcmdn)+(trpnbar*kmcmdn)+cmdnbar*kmtrpn;
	dmyo		= -kmtrpn*kmcmdn*catotal;
	cai_old		    = (2.0/3.0)*sqrt(bmyo*bmyo-3.0*cmyo)*cos(acos((9.0*bmyo*cmyo-2*bmyo*bmyo*bmyo-27*dmyo)/(2.0*pow((bmyo*bmyo-3.0*cmyo),1.5)))/3.0)-bmyo/3.0;  
	
    // Unused ...
	//real caavg       = (casss_old*vsss+cassl*vssl+cai_old*vmyo)/(vsss+vmyo+vssl);
	
 	real camkbound	= camk0*(1-camktrap_old)*1/(1+(kmcam/casss_old));
	
    // Forward Euler
    camktrap_old	= dt*(alphacamk*camkbound*(camkbound+camktrap_old)-betacamk*camktrap_old) + camktrap_old;
	camkactive	= camkbound+camktrap_old; 

    real dvdt = -itot;
    v_old	+= dt*dvdt;

    // Rush-Larsen
    rDY_[1]	= m_old;
	rDY_[2]	= h_old;
	rDY_[3]	= j_old;
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

