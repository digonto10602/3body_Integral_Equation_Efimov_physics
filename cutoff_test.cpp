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
                        double epsilon,
                        double L,
                        double kpoints,
                        double n_max    )
{
    comp ii = {0.0,1.0};
    comp res_I = {0.0,0.0};
    double pi = acos(-1.0);

    double kmin = 0.0;
    comp kmax = pmom(s,0.0,m);
    double real_kmax = real(kmax);
    //double kpoints = 500.0;
    double delk = abs(real_kmax - kmin)/kpoints;
    for(int nx=0;nx<n_max;++nx)
    {
        for(int ny=0;ny<n_max;++ny)
        {
            for(int nz=0;nz<n_max;++nz)
            {
                double n = sqrt(nx*nx + ny*ny + nz*nz);
                if(n<=0) continue;
                if(n>=n_max) continue;
                
                for(double ax=kmin;ax<=real_kmax;ax=ax+delk)
                {
                    for(double ay=kmin;ay<=real_kmax;ay=ay+delk)
                    {
                        for(double az=kmin;az<=real_kmax;az=az+delk)
                        {
                            double amom = sqrt(ax*ax + ay*ay + az*az);

                            double phasespace = delk*delk*delk/pow(2.0*pi,3.0);

                            double amom_n = ax*nx + ay*ny + az*nz;
                            comp expon = exp(ii*amom_n*L);

                            double k = sqrt(kx*kx + ky*ky + kz*kz);
                            comp sigk = sigma_p(s,k,m);
                            comp siga = sigma_p(s,amom,m);

                            double k_plus_a = sqrt(     kx*kx + ky*ky + kz*kz 
                                                    +   ax*ax + ay*ay + az*az   
                                                    +  2.0*kx*ax + 2.0*ky*ay + 2.0*kz*az    );
                            comp sigka = sigma_p(s,k_plus_a,m);

                            comp cutoffs = s0_based_Hfunc_comp_smooth(sigk,s0,m)
                                            *s0_based_Hfunc_comp_smooth(siga,s0,m)
                                            *s0_based_Hfunc_comp_smooth(sigka,s0,m);

                            comp omegak = omega_comp(k,m);
                            comp omegaa = omega_comp(amom, m);
                            comp omegaka = omega_comp(k_plus_a, m);

                            comp denom = 2.0*omegaa*2.0*omegaka*(sqrt(s) - omegak - omegaa - omegaka + ii*epsilon);
                            
                            /*
                            cout<<"nx = "<<nx<<'\t'<<"ny = "<<ny<<'\t'<<"nz = "<<nz<<endl;
                            cout<<"n = "<<n<<endl;
                            
                            cout<<"ax = "<<ax<<'\t'<<"ay = "<<ay<<'\t'<<"az = "<<az<<endl;
                            cout<<"a = "<<amom<<endl;

                            cout<<"exponent = "<<expon<<'\t'<<" cutoffs = "<<cutoffs<<endl;
                            */
                            comp tot = phasespace * expon * cutoffs / denom; 
                            
                            //cout<<"result = "<<tot<<endl;

                            //cout<<"====================================="<<endl;

                            res_I = res_I + tot;
                        }
                    }
                }
            }
        }
    }

    return res_I*0.5;
}

int main()
{

    // Inputs
    double pi = acos(-1.0);
    double a = 2.0;
    double s = 6.5;
    double kx = 0.2;
    double ky = 0.2;
    double kz = 0.2;

    double s0 = 1.0;
    double m = 1.0;
    double epsilon = 0.001;
    double L = 4.0;
    double kpoints = 100.0;
    double n_max = 10.0;

    for(int L=2;L<20;++L)
    {
        comp F = F_ieps_00_smooth(  s, kx, ky, kz, s0, m, epsilon, L, kpoints, n_max);

        cout<<L<<'\t'<<F<<endl;
    }
    return 0;
}