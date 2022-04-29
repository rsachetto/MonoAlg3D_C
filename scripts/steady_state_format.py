import numpy as np

def print_steady_state_format (sv):
    for i in range(len(sv)):
        #print("sv[%d] = %e;" % (i,sv[i]))
        #print("extra_data[%d+offset] = %e;" % (i,sv[i]))
        print("*((real * )((char *) sv + pitch * %d) + threadID) = %e;" % (i,sv[i]))

def get_correct_indexes (sv):
    nlin = np.shape(sv)

    v       =  sv[0] 
    nai     =  sv[1] 
    nass    =  sv[2] 
    ki      =  sv[3] 
    kss     =  sv[4] 
    cai     =  sv[5] 
    cass    =  sv[6] 
    cansr   =  sv[7] 
    cajsr   =  sv[8] 
    m       =  sv[9] 
    hp      =  sv[10]
    h       =  sv[11]
    j       =  sv[12]
    jp      =  sv[13]
    mL      =  sv[14]
    hL      =  sv[15]
    hLp     =  sv[16]
    a       =  sv[17]
    iF      =  sv[18]
    iS      =  sv[19]
    ap      =  sv[20]
    iFp     =  sv[21]
    iSp     =  sv[22]
    d       =  sv[23]
    ff      =  sv[24]
    fs      =  sv[25]
    fcaf    =  sv[26]
    fcas    =  sv[27]
    jca     =  sv[28]
    nca_ss  =  sv[29]
    nca_i   =  sv[30]
    ffp     =  sv[31]
    fcafp   =  sv[32]
    xs1     =  sv[33]
    xs2     =  sv[34]
    Jrel_np =  sv[35]
    CaMKt   =  sv[36]
    C3      =  sv[37]
    C2      =  sv[38]
    C1      =  sv[39]
    O       =  sv[40]
    I       =  sv[41]
    Jrel_p  =  sv[42]

    # Change the columns to the Lucas`s mapping
    sv_new = np.zeros(nlin)

    sv_new[0]  = v        
    sv_new[1]  = CaMKt    
    sv_new[2]  = cass     
    sv_new[3]  = nai      
    sv_new[4]  = nass     
    sv_new[5]  = ki       
    sv_new[6]  = kss      
    sv_new[7]  = cansr    
    sv_new[8]  = cajsr    
    sv_new[9]  = cai      
    sv_new[10] = m        
    sv_new[11] = h        
    sv_new[12] = j        
    sv_new[13] = hp       
    sv_new[14] = jp       
    sv_new[15] = mL       
    sv_new[16] = hL       
    sv_new[17] = hLp      
    sv_new[18] = a        
    sv_new[19] = iF       
    sv_new[20] = iS       
    sv_new[21] = ap       
    sv_new[22] = iFp      
    sv_new[23] = iSp      
    sv_new[24] = d        
    sv_new[25] = ff       
    sv_new[26] = fs       
    sv_new[27] = fcaf     
    sv_new[28] = fcas     
    sv_new[29] = jca      
    sv_new[30] = ffp      
    sv_new[31] = fcafp    
    sv_new[32] = nca_ss   
    sv_new[33] = nca_i    
    sv_new[34] = C1       
    sv_new[35] = C2       
    sv_new[36] = C3       
    sv_new[37] = I        
    sv_new[38] = O        
    sv_new[39] = xs1      
    sv_new[40] = xs2      
    sv_new[41] = Jrel_np  
    sv_new[42] = Jrel_p

    return sv_new

# ToRORd_fkatp steady-state solution after 200 beats using Jakub Tomek Matlab code
torord_baseline_endo_ss = [-8.890585e+01,1.210818e+01,1.210851e+01,1.426206e+02,1.426205e+02,7.455488e-05,6.504164e-05,1.530373e+00,1.528032e+00,7.814592e-04,6.752873e-01,8.313839e-01,8.311938e-01,8.308255e-01,1.585610e-04,5.294475e-01,2.896996e-01,9.419166e-04,9.996194e-01,5.938602e-01,4.799180e-04,9.996194e-01,6.543754e-01,-2.898677e-33,1.000000e+00,9.389659e-01,1.000000e+00,9.999003e-01,9.999773e-01,4.920606e-04,8.337021e-04,1.000000e+00,1.000000e+00,2.471690e-01,1.742987e-04,5.421027e-24,1.107642e-02,9.980807e-01,8.425453e-04,6.962775e-04,3.675442e-04,1.289824e-05,6.407933e-23]
torord_baseline_epi_ss = [-8.917755e+01,1.284260e+01,1.284291e+01,1.429114e+02,1.429113e+02,6.631866e-05,5.767956e-05,1.812268e+00,1.810520e+00,7.370422e-04,6.840260e-01,8.366816e-01,8.366012e-01,8.363958e-01,1.505860e-04,5.412669e-01,3.043382e-01,9.248184e-04,9.996371e-01,9.996342e-01,4.712023e-04,9.996371e-01,9.996366e-01,4.333129e-43,1.000000e+00,9.485160e-01,1.000000e+00,9.999339e-01,9.999822e-01,3.086885e-04,5.303737e-04,1.000000e+00,1.000000e+00,2.308784e-01,1.690386e-04,-1.103286e-23,1.288116e-02,9.982135e-01,8.264829e-04,6.775197e-04,2.730221e-04,9.433146e-06,-6.177055e-22]
torord_baseline_mid_ss = [-8.924177e+01,1.503347e+01,1.503401e+01,1.434407e+02,1.434406e+02,8.177438e-05,6.585066e-05,1.959747e+00,1.963459e+00,7.269124e-04,6.860578e-01,8.379059e-01,8.377164e-01,8.372100e-01,1.487602e-04,5.350003e-01,2.851164e-01,9.208259e-04,9.996411e-01,5.673539e-01,4.691672e-04,9.996412e-01,6.265825e-01,-4.922960e-40,1.000000e+00,9.200354e-01,1.000000e+00,9.997888e-01,9.999665e-01,5.161178e-04,1.189422e-03,1.000000e+00,1.000000e+00,2.650323e-01,1.678628e-04,2.091039e-25,1.922391e-02,9.979358e-01,8.225453e-04,6.917041e-04,5.316232e-04,1.835276e-05,2.438403e-23]

# Trovato steady-state solution after 200 beats using MonoAlg3D
trovato_baseline_ss = [-86.2436,0.00563079,0.000105437,9.46644,9.46626,9.46627,142.526,142.527,142.527,4.51062e-05,0.000105684,1.30622,1.28855,1.30783,1.94256e-05,0,0.00670099,0.776682,0.776592,0.77863,0.562506,0.778611,0.0002629,0.445838,0.222675,0.000282192,0.631896,0.989626,7.73809e-09,1,0.915897,1,1,0.999969,1,1,0.00628123,0.00032388,0.993688,0.000566292,0.586241,0.209664,0.000233827,0.208843,0.997186,0.475622]

#torord_baseline_endo_ss = get_correct_indexes(torord_baseline_endo_ss)
#print_steady_state_format(torord_baseline_endo_ss)
#print()

#torord_baseline_epi_ss = get_correct_indexes(torord_baseline_epi_ss)
#print_steady_state_format(torord_baseline_epi_ss)
#print()

#torord_baseline_mid_ss = get_correct_indexes(torord_baseline_mid_ss)
#print_steady_state_format(torord_baseline_mid_ss)
#print()

print_steady_state_format(trovato_baseline_ss)
print()