//
// Created by yi_min on 26.03.18.
//

#include "seiler2_u1.h"
registerMooseObject("seiler2App",seiler2_u1);
template <>
InputParameters validParams<seiler2_u1>()
{
    InputParameters params=validParams<TimeKernel>();

    params.addRequiredParam<Real>("f1","description of f1");

    params.addRequiredCoupledVar("u2","displacement uy"); 
    params.addCoupledVar("d", "phase field value");

    return params;
}

seiler2_u1::seiler2_u1(const InputParameters &parameters)
        :Kernel(parameters),
         _f1(getParam<Real>("f1")),
         _couple_u2(coupledValue("u2")),
         _couple_d(coupledValue("d")),
         _grad_couple_u2(coupledGradient("u2")),
         _u2_var(coupled("u2")),
         _d_var(coupled("d")),
         _sigma1(getMaterialProperty<Real>("sigma1")),
         _sigma2(getMaterialProperty<Real>("sigma2")),
         _sigma3(getMaterialProperty<Real>("sigma3")),
         _sigmap1(getMaterialProperty<Real>("sigmap1")),
         _sigmap2(getMaterialProperty<Real>("sigmap2")),
         _sigmap3(getMaterialProperty<Real>("sigmap3")),
         _c11(getMaterialProperty<Real>("c11")),
         _c12(getMaterialProperty<Real>("c12")),
         _c13(getMaterialProperty<Real>("c13")),
         _c21(getMaterialProperty<Real>("c21")),
         _c22(getMaterialProperty<Real>("c22")),
         _c23(getMaterialProperty<Real>("c23")),
         _c31(getMaterialProperty<Real>("c31")),
         _c32(getMaterialProperty<Real>("c32")),
         _c33(getMaterialProperty<Real>("c33"))

{}

Real seiler2_u1::computeQpResidual()
{
    // compute R_u1
    return 
    _sigma1[_qp]*_grad_test[_i][_qp](0)+_sigma3[_qp]*_grad_test[_i][_qp](1)+_f1*_test[_i][_qp];
    //_grad_u[_qp](0)*_grad_test[_i][_qp](0)+_grad_u[_qp](1)*_grad_test[_i][_qp](1)+_grad_u[_qp](2)*_grad_test[_i][_qp](2);
    // u always represent current kernel's variable, here it represents phi
}

Real seiler2_u1::computeQpJacobian()
{
    // compute K_u1_u1
    return (_c11[_qp]*_grad_test[_i][_qp](0)+_c31[_qp]*_grad_test[_i][_qp](1))*_grad_phi[_j][_qp](0)+(_c13[_qp]*_grad_test[_i][_qp](0)+_c33[_qp]*_grad_test[_i][_qp](1))*_grad_phi[_j][_qp](1);
}

// Now we calculate the coupling term K(phi,m3), if we have coupling between phi and m3
Real seiler2_u1::computeQpOffDiagJacobian(unsigned int jvar)
{

    if(jvar==_u2_var)
    {
         // K_u1_u2
        return (_c13[_qp]*_grad_test[_i][_qp](0)+_c33[_qp]*_grad_test[_i][_qp](1))*_grad_phi[_j][_qp](0)+(_c12[_qp]*_grad_test[_i][_qp](0)+_c32[_qp]*_grad_test[_i][_qp](1))*_grad_phi[_j][_qp](1);

    }
    
    if(jvar==_d_var)
    {
         // K_u1_d
        return -2.0*(1.0-_couple_d[_qp])*(_sigmap1[_qp]*_grad_test[_i][_qp](0)+_sigmap3[_qp]*_grad_test[_i][_qp](1))*_phi[_j][_qp];
//                        
    }

    return 0.0; // if not the Kcmu ,return zero contribution
}
    
