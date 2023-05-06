//
//

#ifndef seiler2App_seiler2_d_H
#define seiler2App_seiler2_d_H

#include "TimeKernel.h"

#include "RankTwoTensor.h"
#include "RankFourTensor.h"

class seiler2_d;

template <>
InputParameters validParams<seiler2_d>();

class seiler2_d:public TimeKernel
{
public:
    seiler2_d(const InputParameters &parameters);

protected:
    virtual Real computeQpResidual() override ;
    virtual Real computeQpJacobian() override ;
    virtual Real computeQpOffDiagJacobian(unsigned int jvar) override ;
    
        // Input parameter
    const Real &_G,&_l0;
    const VariableValue &_couple_u1,&_couple_u2;
    unsigned int _u1_var,_u2_var;
    const MaterialProperty<Real> &_pc;
    const MaterialProperty<Real> &_hist,&_sigma1,&_sigma2,&_sigma3,&_sigmap1,&_sigmap2,&_sigmap3,&_sigmapfc1,&_sigmapfc2,&_sigmapfc3,&_alpha;

    // this is the kernel for c
    // we need gradient of mu, diffusivity
//    const VariableValue &_couple_m2_dot,&_couple_m3_dot,&_couple_m2_dotdu,&_couple_m3_dotdu;

    //const Real Gamm0=1.76086E11,mu0=4*3.141592653E-7;
};



#endif //MFMAPP_MFM_POISSON_H
