//
// Created by yi_min on 26.03.18.
//

#ifndef seiler2App_seiler2_u1_H
#define seiler2App_seiler2_u1_H

#include "TimeKernel.h"

class seiler2_u1;

template <>
InputParameters validParams<seiler2_u1>();

class seiler2_u1:public Kernel
{
public:
    seiler2_u1(const InputParameters &parameters);

protected:
    virtual Real computeQpResidual() override ;
    virtual Real computeQpJacobian() override ;
    virtual Real computeQpOffDiagJacobian(unsigned int jvar) override ;
    
        // Input parameter
    const Real &_f1;
    const VariableValue &_couple_u2,&_couple_d;
    const VariableGradient &_grad_couple_u2;
    unsigned int _u2_var,_d_var;
    const MaterialProperty<Real> &_sigma1, &_sigma2, &_sigma3,&_sigmap1,&_sigmap2,&_sigmap3,&_c11,&_c12,&_c13,&_c21,&_c22,&_c23,&_c31,&_c32,&_c33;

    // this is the kernel for c
    // we need gradient of mu, diffusivity
//    const VariableValue &_couple_m2_dot,&_couple_m3_dot,&_couple_m2_dotdu,&_couple_m3_dotdu;

    //const Real Gamm0=1.76086E11,mu0=4*3.141592653E-7;
};



#endif //MFMAPP_MFM_POISSON_H
