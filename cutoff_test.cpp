//macos:g++-11 printer.cpp -o printer -O3 -std=c++14 -I ~/homebrew/include/eigen3/ -L ~/homebrew/lib/
//wahab:g++ printer.cpp -o printer -O3 -std=c++14 -I /cm/shared/applications/eigen/3.3.7/include/eigen3/
//wahab_gpu:nvcc -I /cm/shared/applications/cuda-toolkit/10.0.130/include -I /cm/shared/applications/eigen/3.3.7/include/eigen3/ printer_gpu_momrep.cpp  -O3 -std=c++14 -L /cm/shared/applications/cuda-toolkit/10.0.130/lib64 -lcublas -lcusolver -lcudart -o printer_gpu
//wahab_gpu_omp:nvcc -I /cm/shared/applications/cuda-toolkit/10.0.130/include -I /cm/shared/applications/eigen/3.3.7/include/eigen3/ printer_gpu_momrep.cpp  -O3 -std=c++14 -L /cm/shared/applications/cuda-toolkit/10.0.130/lib64 -lcublas -lcusolver -lcudart -lgomp -Xcompiler -fopenmp -o printer_gpu
//ubuntupc_omp:nvcc117 -I/usr/include/eigen3/ -I/opt/nvidia/hpc_sdk/Linux_x86_64/2022/math_libs/include -I/opt/nvidia/hpc_sdk/Linux_x86_64/2022/cuda/include -L/opt/nvidia/hpc_sdk/Linux_x86_64/2022/math_libs/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/2022/cuda/lib64 printer_gpu_momrep.cpp  -O3 -std=c++14 -lcublas -lcusolver -lcudart -lgomp -Xcompiler -fopenmp -o printer_gpu
//ubuntupc_omp1:nvcc -I/usr/include/eigen3/ -I/opt/nvidia/hpc_sdk/Linux_x86_64/2022/math_libs/include -I/opt/nvidia/hpc_sdk/Linux_x86_64/2022/cuda/include -L/opt/nvidia/hpc_sdk/Linux_x86_64/2022/math_libs/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/2022/cuda/lib64 printer_gpu_momrep.cpp  -O3 -std=c++14 -lcublas -lcusolver -lcudart -lgomp -Xcompiler -fopenmp -o printer_gpu

#include<bits/stdc++.h>
#include "functions_momrep_based.h"
#include "momentum_vector_maker.h"
#include "integralequation_momrep_based.h"
#include "interpolator_momrep.h"
#include "solvers.h"
#include<omp.h>
#include <Eigen/Eigenvalues> //for vertex factor calculations
//#include "LUSolver_modified.h"
#include "fastGL_based_functions.h"
#include "SA_method_functions.h"
#include <sys/time.h>

using namespace std;

typedef complex<double> comp;


comp s0_based_Jfunc_comp_smooth(    double s0, 
                                    comp x  )
{
    if(real(x)<=0.0) return 0.0;
    else if(real(x)>0.0 && real(x)<s0) return exp((-1.0/x)*exp(-1.0/(1.0-x)));
    else return 1.0;
}

comp s0_based_Jfunc_comp_hard(      double s0, 
                                    comp x  )
{
    if(real(x)<=s0) return 0.0;
    //else if(real(x)>0.0 && real(x)<s0) return exp((-1.0/x)*exp(-1.0/(1.0-x)));
    else return 1.0;
}

comp s0_based_Hfunc_comp_smooth(    comp sigmap,
                                    double s0,
                                    double m    )
{
    return s0_based_Jfunc_comp_smooth(s0, sigmap/(4.0*m*m));
}

comp s0_based_Hfunc_comp_hard(      comp sigmap,
                                    double s0,
                                    double m    )
{
    return s0_based_Jfunc_comp_hard(s0, sigmap/(4.0*m*m));
}


comp K2_scattlength(    comp s,
                        double scattering_length)
{
    comp pi = acos(-1.0);

    return -8.0*pi*sqrt(s)*scattering_length;
}

comp G_00_smooth(   comp s,
                    double px,
                    double py,
                    double pz,
                    double kx,
                    double ky,
                    double kz,
                    double s0,
                    double m    )

{
    comp p = sqrt(px*px + py*py + pz*pz);
    comp k = sqrt(kx*kx + ky*ky + kz*kz);
    comp p_plus_k = sqrt(       px*px + py*py + pz*pz
                            +   kx*kx + ky*ky + kz*kz 
                            +  2.0*px*kx + 2.0*py*ky + 2.0*pz*kz    );
    comp s2p = sigma_p(s,p,m);
    comp s2k = sigma_p(s,k,m);

    comp cutoff1 = s0_based_Hfunc_comp_smooth(s2p,s0,m);
    comp cutoff2 = s0_based_Hfunc_comp_smooth(s2k,s0,m);

    comp omegap = omega_comp(p,m);
    comp omegak = omega_comp(k,m);
    comp omegakp = omega_comp(p_plus_k,m);

    comp res = cutoff1*cutoff2/(2.0*omegakp*(sqrt(s) - omegak - omegap - omegakp));

    return res;

}

comp G_00_hard(   comp s,
                    double px,
                    double py,
                    double pz,
                    double kx,
                    double ky,
                    double kz,
                    double s0,
                    double m    )

{
    comp p = sqrt(px*px + py*py + pz*pz);
    comp k = sqrt(kx*kx + ky*ky + kz*kz);
    comp p_plus_k = sqrt(       px*px + py*py + pz*pz
                            +   kx*kx + ky*ky + kz*kz 
                            +  2.0*px*kx + 2.0*py*ky + 2.0*pz*kz    );
    comp s2p = sigma_p(s,p,m);
    comp s2k = sigma_p(s,k,m);

    comp cutoff1 = s0_based_Hfunc_comp_hard(s2p,s0,m);
    comp cutoff2 = s0_based_Hfunc_comp_hard(s2k,s0,m);

    comp omegap = omega_comp(p,m);
    comp omegak = omega_comp(k,m);
    comp omegakp = omega_comp(p_plus_k,m);

    comp res = cutoff1*cutoff2/(2.0*omegakp*(sqrt(s) - omegak - omegap - omegakp));

    return res;

}

comp F_ieps_00_smooth(  comp s,
                        double kx,
                        double ky,
                        double kz, 
                        double s0,
                        double m,
                        double L,
                        double n_max    )
{
    
}