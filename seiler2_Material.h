//
//

#ifndef seiler2APP_seiler2_MATERIAL_H
#define seiler2APP_seiler2_MATERIAL_H

#include "Material.h"
#include "DerivativeMaterialInterface.h"

//class seiler2_Material;

//template <>
//InputParameters validParams<seiler2_Material>();

class seiler2_Material:public DerivativeMaterialInterface<Material>
{
public:
    static InputParameters validParams();
    
    seiler2_Material(const InputParameters &parameters);
    
protected:
    virtual void computeQpProperties() override;
    
    virtual void initQpStatefulProperties();
    
    const Real &_E,&_nu,&_cm11,&_cm12,&_cm22,&_cm33,&_K,&_n,&_sigmaf,&_epsilonf,&_a0,&_b,&_c,&_alpha0,&_xi,&_sigmau,&_D0,&_alpha1,&_R;
    
    unsigned int _model;
    
    const MaterialProperty<Real> &_cycle,&_Gd,&_dGddc;
    
    const VariableValue &_couple_u1,&_couple_u2,&_couple_d;
    // next line is for gradients of coupled variables
    const VariableGradient &_grad_couple_u1,&_grad_couple_u2;
     
    MaterialProperty<Real> &_hel,&_epsilon1,&_epsilon2,&_epsilon3,&_sigma1,&_sigma2,&_sigma3,&_Pswt1,&_Pswt2,&_Pswt3,&_sigmap1,&_sigmap2,&_sigmap3,&_sigmapfc1,&_sigmapfc2,&_sigmapfc3,&_c11,&_c12,&_c13,&_c21,&_c22,&_c23,&_c31,&_c32,&_c33,&_alpha,&_D,&_Ni,&_sigmareva,&_sigmamax,&_sigma_max,&_hist,&_sigmabar;
    
    const MaterialProperty<Real> &_sigmamax_old,& _hist_old,&_D_old;
    
    const bool _print_convergence_message;
};

#endif //seiler2APP_seiler2_MATERIAL_H
