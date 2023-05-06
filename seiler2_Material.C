//
//

#include "seiler2_Material.h"

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <Eigen/Dense>
// #include <algorithm>
using namespace std;

registerMooseObject("seiler2App",seiler2_Material);

//template <>
InputParameters seiler2_Material::validParams()
{
    InputParameters params = Material::validParams();

    params.addRequiredParam<Real>("E","Young's modulu");
    params.addRequiredParam<Real>("nu","Possion's ratio");
    params.addRequiredParam<Real>("K","cyclic hardening coefficient");
    params.addRequiredParam<Real>("n","cyclic hardening exponent");
    params.addRequiredParam<Real>("sigmaf","SWT parameter");
    params.addRequiredParam<Real>("epsilonf","SWT parameter");
    params.addRequiredParam<Real>("a0","degradation function constant");
    params.addRequiredParam<Real>("b","SWT parameter");
    params.addRequiredParam<Real>("c","SWT parameter");
    params.addRequiredParam<Real>("alpha0","lifetime parameter");
    params.addRequiredParam<Real>("D0","damage0");
    params.addRequiredParam<Real>("alpha1","lifetime parameter1");
    params.addRequiredParam<Real>("xi","lifetime parameter");
    params.addRequiredParam<Real>("sigmau","yield strength");
    params.addRequiredParam<Real>("cm11","Elasticity modulu");
    params.addRequiredParam<Real>("cm12","Elasticity modulu");
    params.addRequiredParam<Real>("cm22","Elasticity modulu");
    params.addRequiredParam<Real>("cm33","Elasticity modulu");
    params.addRequiredParam<Real>("R","Stress ratio");
    params.addRequiredParam<unsigned int>("model","model=1 (strain method), 2 (stresstrain method),3 (stress method)");
    
    params.addCoupledVar("u1","displacement ux"); 
    params.addCoupledVar("u2","displacement uy"); 
    params.addCoupledVar("d","fracture phase field value");
    params.addRequiredParam<std::string>("pro_cycle","get the cycle number");
//     params.addRequiredParam<std::string>("pro_rsigma1","get the residual stress_xx");
//     params.addRequiredParam<std::string>("pro_rsigma2","get the residual stress_yy");
//     params.addRequiredParam<std::string>("pro_rsigma3","get the residual stress_xy");
    params.addParam<MaterialPropertyName>("g_name", "g_d", "degradation function g(d)");
    params.addParam<bool>(
      "print_state_variable_convergence_error_messages",
      false,
      "Whether or not to print warning messages from the crystal plasticity specific convergence "
      "checks on the stress measure and general constitutive model quantinties.");

    return params;
}

seiler2_Material::seiler2_Material(const InputParameters &parameters)
    : DerivativeMaterialInterface<Material>(parameters),        
//:Material(parameters),
         // Input parameters
         _E(getParam<Real>("E")),
         _nu(getParam<Real>("nu")),
         _cm11(getParam<Real>("cm11")),
         _cm12(getParam<Real>("cm12")),
         _cm22(getParam<Real>("cm22")),
         _cm33(getParam<Real>("cm33")),
         _K(getParam<Real>("K")),
         _n(getParam<Real>("n")),
         _sigmaf(getParam<Real>("sigmaf")),
         _epsilonf(getParam<Real>("epsilonf")),
         _a0(getParam<Real>("a0")),
         _b(getParam<Real>("b")),
         _c(getParam<Real>("c")),
         _alpha0(getParam<Real>("alpha0")),
         _xi(getParam<Real>("xi")),
         _sigmau(getParam<Real>("sigmau")),
         _D0(getParam<Real>("D0")),
         _alpha1(getParam<Real>("alpha1")),
         _R(getParam<Real>("R")),
         _model(getParam<unsigned int>("model")),
         
         _cycle(getMaterialProperty<Real>(parameters.get<std::string>("pro_cycle"))),

         _Gd(declareProperty<Real>(getParam<MaterialPropertyName>("g_name"))),
         _dGddc(getMaterialPropertyDerivative<Real>("g_name", getVar("d", 0)->name())),
         
         _couple_u1(coupledValue("u1")),
         _couple_u2(coupledValue("u2")),
         _couple_d(coupledValue("d")),
         
         _grad_couple_u1(coupledGradient("u1")),
         _grad_couple_u2(coupledGradient("u2")),
         // Parameters to be calculated with this Materials block
         //_hen(declareProperty<RealVectorValue>("hen")),
         _hel(declareProperty<Real>("hel")),
         
         _epsilon1(declareProperty<Real>("epsilon1")),
         _epsilon2(declareProperty<Real>("epsilon2")),
         _epsilon3(declareProperty<Real>("epsilon3")),
         
         _sigma1(declareProperty<Real>("sigma1")),
         _sigma2(declareProperty<Real>("sigma2")),
         _sigma3(declareProperty<Real>("sigma3")),
         
         _Pswt1(declareProperty<Real>("Pswt1")),
         _Pswt2(declareProperty<Real>("Pswt2")),
         _Pswt3(declareProperty<Real>("Pswt3")),
         
         _sigmap1(declareProperty<Real>("sigmap1")),
         _sigmap2(declareProperty<Real>("sigmap2")),
         _sigmap3(declareProperty<Real>("sigmap3")),
         
//          _sigmapd1(declareProperty<Real>("sigmapd1")),
//          _sigmapd2(declareProperty<Real>("sigmapd2")),
//          _sigmapd3(declareProperty<Real>("sigmapd3")),
         
         _sigmapfc1(declareProperty<Real>("sigmapfc1")),
         _sigmapfc2(declareProperty<Real>("sigmapfc2")),
         _sigmapfc3(declareProperty<Real>("sigmapfc3")),
         
         _c11(declareProperty<Real>("c11")),
         _c12(declareProperty<Real>("c12")),
         _c13(declareProperty<Real>("c13")),
         _c21(declareProperty<Real>("c21")),
         _c22(declareProperty<Real>("c22")),
         _c23(declareProperty<Real>("c23")),
         _c31(declareProperty<Real>("c31")),
         _c32(declareProperty<Real>("c32")),
         _c33(declareProperty<Real>("c33")),
        
         
         _alpha(declareProperty<Real>("alpha")),
         _D(declareProperty<Real>("D")),
         _Ni(declareProperty<Real>("Ni")),
         
         
         _sigmareva(declareProperty<Real>("sigmareva")),
         _sigmamax(declareProperty<Real>("sigmamax")), //sigmamax=max{sigmareva}
         _sigma_max(declareProperty<Real>("sigma_max")), //sigma_max=2/(1-R)*sigmareva
         _hist(declareProperty<Real>("hist")),
         _sigmabar(declareProperty<Real>("sigmabar")),
         _sigmamax_old(getMaterialPropertyOld<Real>("sigmamax")),
         _hist_old(getMaterialPropertyOld<Real>("hist")),
         _D_old(getMaterialPropertyOld<Real>("D")),
         
         _print_convergence_message(getParam<bool>("print_state_variable_convergence_error_messages"))

{}

void seiler2_Material::initQpStatefulProperties()
{
    _sigmamax[_qp] =0.0; 
    _hist[_qp] =0.0;
    _D[_qp]=0.0;

}

void seiler2_Material::computeQpProperties()
{    
    _epsilon1[_qp]=_grad_couple_u1[_qp](0);
    _epsilon2[_qp]=_grad_couple_u2[_qp](1);
    _epsilon3[_qp]=0.5*_grad_couple_u2[_qp](0)+0.5*_grad_couple_u1[_qp](1);
    
    double lambda,mu,sigmaa,delta_sigma,epsilonbar,epsilona,sigmamean,test,psie_p,psie_n,p_tree,n_tree,p_tree1,n_tree1,tree,tree1,Kn,p_tree_sign,n_tree_sign,f1,f2,f3,pf,fsigma,fsigma1,dfds,f21,Pswt21;
    float Iden[3],epsilon_dev[3],sig_p[3],sig_pd[3],sig[3],sig_n[3],IxI[3][3],Cmat[3][3],Cm[3][3],p_edev_p_e[3][3],rsigma[3];
    int i,j;
    //i=1;
    //_alpha[_qp]=1.0;
    
    lambda=_E*_nu/((1.0+_nu)*(1.0-2.0*_nu)); //plane strain
    mu=_E/(2.0*(1.0+_nu));
    Kn=lambda+2.0*mu/2.0;
    tree=_epsilon1[_qp]+_epsilon2[_qp];
    
//     rsigma[0]=_rsigma1[_qp];
//     rsigma[1]=_rsigma2[_qp];
//     rsigma[2]=_rsigma3[_qp];
//     
//     Cm[0][0]=_cm11;
//     Cm[1][1]=_cm22;
//     Cm[2][2]=_cm33;
//     Cm[0][1]=Cm[1][0]=_cm12;
//     Cm[0][2]=Cm[2][0]=Cm[1][2]=Cm[2][1]=0.0;
    
//     Eigen::Matrix<float,3,3> c_matrix;
//     c_matrix << Cm[0][0],Cm[0][1],Cm[0][2],Cm[1][0],Cm[1][1],Cm[1][2],Cm[2][0],Cm[2][1],Cm[2][2];
//     Eigen::Matrix<float,3,1> eigen_strain, res_stress;
//     res_stress << rsigma[0],rsigma[1],rsigma[2];
//     eigen_strain = c_matrix.inverse()*res_stress;
    

//     tree1=_epsilon1[_qp]+_epsilon2[_qp]+eigen_strain[0]+eigen_strain[1];
    
    //cout << res_stress[0] <<" " << res_stress[1] <<" " << eigen_strain[0] <<" " << eigen_strain[1] <<endl;
    
    Iden[0]=1.0;
    Iden[1]=1.0;
    Iden[2]=0.0;

    epsilon_dev[0]=_epsilon1[_qp]-1.0/2.0*tree;
    epsilon_dev[1]=_epsilon2[_qp]-1.0/2.0*tree;
    epsilon_dev[2]=_epsilon3[_qp];
    
//calculate elastic energy    
    if(tree>=0.0)
    {
        p_tree=tree;
        n_tree=0.0;
    }
    else
    {
        p_tree=0.0;
        n_tree=tree;
    }
    
    psie_p=0.5*Kn*pow(p_tree,2)+mu*(pow(epsilon_dev[0],2)+pow(epsilon_dev[1],2)+2.0*pow(epsilon_dev[2],2));
    psie_n=0.5*Kn*pow(n_tree,2);
    _hel[_qp]=0.5*lambda*pow(_epsilon1[_qp]+_epsilon2[_qp],2)+mu*(pow(_epsilon1[_qp],2)+pow(_epsilon2[_qp],2)+2.0*pow(_epsilon3[_qp],2));
    
//calculate stress    
    if(tree>=0.0)
    {
        p_tree=tree;
        n_tree=0.0;
        p_tree_sign=1.0;
        n_tree_sign=0.0;
    }
    else
    {
        p_tree=0.0;
        n_tree=tree;
        p_tree_sign=0.0;
        n_tree_sign=1.0;
    }
    
    for(i=1;i<=3;i++)
    {
         sig_p[i-1]=Kn*p_tree*Iden[i-1]+2.0*mu*epsilon_dev[i-1];
         sig_n[i-1]=Kn*n_tree*Iden[i-1];
         sig_pd[i-1]=(pow(1.0-_couple_d[_qp],2)+_a0)*sig_p[i-1];
         sig[i-1]=(pow(1.0-_couple_d[_qp],2)+_a0)*sig_p[i-1]+sig_n[i-1]; 
    }
         //cout << sig[0] <<" " << sig[1] <<endl;
    
//calculate Cij    
    for(i=1;i<=3;i++)
    {
        for(j=1;j<=3;j++)
        {
            IxI[i-1][j-1]=Iden[i-1]*Iden[j-1];
        }
    }
    
    for(i=1;i<=3;i++)
    {
        for(j=1;j<=3;j++)
        {
            if(i==j && i<=2)
            {
                p_edev_p_e[i-1][j-1]=1.0/2.0;
            }
            else if(i==j && i>2)
            {
                p_edev_p_e[i-1][j-1]=1.0;
            }
            else if(i>2 || j>2)
            {
                p_edev_p_e[i-1][j-1]=0.0;
            }
            else
            {
                p_edev_p_e[i-1][j-1]=-1.0/2.0;
            }
        }
    }
    
    for(i=1;i<=3;i++)
    {
        for(j=1;j<=3;j++)
        {
            Cmat[i-1][j-1]=(pow(1.0-_couple_d[_qp],2)+_a0)*(Kn*p_tree_sign*IxI[i-1][j-1]+2.0*mu*p_edev_p_e[i-1][j-1])+Kn*n_tree_sign*IxI[i-1][j-1];
        }
    }
    
    _sigma1[_qp]=sig[0];  
    _sigma2[_qp]=sig[1];
    _sigma3[_qp]=sig[2];
    
    _sigmap1[_qp]=sig_p[0];
    _sigmap2[_qp]=sig_p[1];
    _sigmap3[_qp]=sig_p[2];
    
    _c11[_qp]=Cmat[0][0];
    _c12[_qp]=Cmat[0][1];
    _c13[_qp]=Cmat[0][2];
    _c21[_qp]=Cmat[1][0];
    _c22[_qp]=Cmat[1][1];
    _c23[_qp]=Cmat[1][2];
    _c31[_qp]=Cmat[2][0];
    _c32[_qp]=Cmat[2][1];
    _c33[_qp]=Cmat[2][2];

    if (psie_p >= _hist_old[_qp])
    {
      _hist[_qp] = psie_p;
      _sigmapfc1[_qp]=_sigmap1[_qp];
      _sigmapfc2[_qp]=_sigmap2[_qp];
      _sigmapfc3[_qp]=_sigmap3[_qp];
    }
    else
    {
      _hist[_qp] = _hist_old[_qp];
      _sigmapfc1[_qp]=0.0;
      _sigmapfc2[_qp]=0.0;
      _sigmapfc3[_qp]=0.0;
    }
    
//     if (_hel[_qp] >= _hist_old[_qp])
//     {
//       _hist[_qp] = _hel[_qp];
//     }
//     else
//     {
//       _hist[_qp] = _hist_old[_qp];
//     }
    
    _sigmabar[_qp]=pow(pow(_sigma1[_qp],2)+pow(_sigma2[_qp],2)-_sigma1[_qp]*_sigma2[_qp]+3*pow(_sigma3[_qp],2),0.5);  
    
    
    //calculate sigma_a (yes) or sigma_max??
    _sigmareva[_qp]=_sigmabar[_qp];
    
    //fsigma1=pow(_sigmabar[_qp],2)/_sigmareva[_qp]-_sigmareva[_qp]-2.0*_E/(1+_n)*pow(_sigmareva[_qp]/_K,1.0/_n);  //Glinka's rule + Ramberg-Osgood model

//   while(abs(fsigma1) >0.00001)
//    {
//        fsigma1=pow(_sigmabar[_qp],2)/_sigmareva[_qp]-_sigmareva[_qp]-2.0*_E/(1+_n)*pow(_sigmareva[_qp]/_K,1.0/_n);
    
//        dfds=-1.0*pow(_sigmabar[_qp],2)/pow(_sigmareva[_qp],2) - 1.0 - 2.0*_E/(1.0+_n)/(_n*_K)*pow(_sigmareva[_qp]/_K,1.0/_n-1.0);
        
//        _sigmareva[_qp]=_sigmareva[_qp]-fsigma1/dfds;
        
//    }
    
    fsigma1=pow(_sigmabar[_qp],2)/_sigmareva[_qp]-(_sigmareva[_qp] + _E*pow(_sigmareva[_qp]/_K,1/_n));  //Neuber's rule + Ramberg-Osgood model
    
    while(abs(fsigma1) >0.0001)
    {
        fsigma1=pow(_sigmabar[_qp],2)/_sigmareva[_qp] - (_sigmareva[_qp] + _E*pow(_sigmareva[_qp]/_K,1/_n));  
    
        dfds=-1*pow(_sigmabar[_qp],2)/pow(_sigmareva[_qp],2) - 1 - _E/(_K*_n)*pow(_sigmareva[_qp]/_K,1/_n -1);
        
        _sigmareva[_qp]=_sigmareva[_qp]-fsigma1/dfds;
        
    }
    
    if( _sigmareva[_qp] >= _sigmamax_old[_qp])
    {
        _sigmamax[_qp]=_sigmareva[_qp];
    }
    else
    { 
        _sigmamax[_qp]=_sigmamax_old[_qp];
    }

    
    //calculate delta_sigma=2.0*sigma_a=2.0*_sigmareva[_qp]*(1-_R)/2;
    delta_sigma=2.0*_sigmamax[_qp]; //for R=-1
    sigmamean=(1+_R)/(1-_R)*_sigmamax[_qp];
    _sigma_max[_qp]=2.0*_sigmamax[_qp]/(1-_R);
    
    //calculate 1/2*delta_epsilon=epsilona for R=-1 
    epsilona=_sigmamax[_qp]/_E+pow(_sigmamax[_qp]/_K,1.0/_n);
    
    //calculate the damage for various models
    _Pswt1[_qp]=1.2*(_sigmamax[_qp]/_E+ pow(_sigmamax[_qp]/_K,1/_n));
    //Pswt2=_sigmamax[_qp]*_E*epsilona*2.0/(1-_R);
    _Pswt2[_qp]=_sigma_max[_qp]*(_sigmamax[_qp]/_E+pow(_sigmamax[_qp]/_K,1.0/_n));
    _Pswt3[_qp]=_sigma_max[_qp];
     Pswt21=_sigma_max[_qp]*_E*epsilona*2.0;
    
    //_Ni[_qp]=1e6;
    if (_sigmamax[_qp]>_sigmau)
    {
        _Ni[_qp]=1;
    }
    else
    {
        _Ni[_qp]=1e6;
    }
    
    //calculte Ni for one cycle
    f1=_E*_Pswt1[_qp]-_sigmaf*(1-sigmamean/_sigmaf)*pow(2.0*_Ni[_qp],_b)-_E*_epsilonf*pow(2.0*_Ni[_qp],_c); //model1 strain method//修正为考虑平均应力的莫罗公式。
    //f1=_E*Pswt1-_sigmaf*pow(2.0*_Ni[_qp],_b)-_E*_epsilonf*pow(2.0*_Ni[_qp],_c);
    //f20=pow(_sigmaf,2)+_sigmaf*_epsilonf*_E*pow(2.0*_Ni[_qp],_c-_b)-Pswt2*pow(2.0*_Ni[_qp],-2.0*_b); //model2 stressstrain method
    f2=_Pswt2[_qp]-pow(_sigmaf,2)/_E*pow(2.0*_Ni[_qp],2.0*_b)-_sigmaf*_epsilonf*pow(2.0*_Ni[_qp],_b+_c);
    f21=pow(_sigmaf,2)+_sigmaf*_epsilonf*_E*pow(2.0*_Ni[_qp],_c-_b)-Pswt21*pow(2.0*_Ni[_qp],-2.0*_b); //other split
    //f3=pow(_sigmaf,2)+_sigmaf*_epsilonf*_E*pow(2.0*_Ni[_qp],_c-_b)-_Pswt2[_qp]*pow(2.0*_Ni[_qp],-2.0*_b);
    f3=_Pswt3[_qp]-_sigmaf*pow(2.0*_Ni[_qp],_b); //model3 stress method
    
   if (_model==1) //model1
   {
    
    while(abs(f1) >0.00001)
    {
        f1=_E*_Pswt1[_qp]-_sigmaf*(1-sigmamean/_sigmaf)*pow(2.0*_Ni[_qp],_b)-_E*_epsilonf*pow(2.0*_Ni[_qp],_c);
        //f1=_E*Pswt1-_sigmaf*pow(2.0*_Ni[_qp],_b)-_E*_epsilonf*pow(2.0*_Ni[_qp],_c);
    
        pf=-_sigmaf*(1-sigmamean/_sigmaf)*_b*2.0*pow(2.0*_Ni[_qp],_b-1)-_E*_epsilonf*_c*2.0*pow(2.0*_Ni[_qp],_c-1);
        
        _Ni[_qp]=_Ni[_qp]-f1/pf;
//         if(_Ni[_qp] < 0.0 || _Ni[_qp] > 1e8)
//         {
//             _Ni[_qp]=1e8;
//             break;
//         }
    }
   }
   
   else if (_model==2)
   {
    while(fabs(f2) >0.0001)
    {
        f2=_Pswt2[_qp]-pow(_sigmaf,2)/_E*pow(2.0*_Ni[_qp],2.0*_b)-_sigmaf*_epsilonf*pow(2.0*_Ni[_qp],_b+_c);
        //f20=pow(_sigmaf,2)+_sigmaf*_epsilonf*_E*pow(2.0*_Ni[_qp],_c-_b)-Pswt2*pow(2.0*_Ni[_qp],-2.0*_b);
        
        pf=-4.0*_b*pow(_sigmaf,2)/_E*pow(2.0*_Ni[_qp],2.0*_b-1)-2.0*(_b+_c)*_sigmaf*_epsilonf*pow(2.0*_Ni[_qp],_b+_c-1);
        //pf0=_sigmaf*_epsilonf*_E*(_c-_b)*pow(2.0*_Ni[_qp],_c-_b-1)*2.0+4.0*Pswt2*_b*pow(2.0*_Ni[_qp],-2.0*_b-1);
        
        _Ni[_qp]=_Ni[_qp]-f2/pf;
    }
   }

   else if (_model==3)
   {
         
    while(abs(f3) >0.0001)
    {
        
        f3=_Pswt3[_qp]-_sigmaf*pow(2.0*_Ni[_qp],_b);
        //f3=pow(_sigmaf,2)+_sigmaf*_epsilonf*_E*pow(2.0*_Ni[_qp],_c-_b)-_Pswt2[_qp]*pow(2.0*_Ni[_qp],-2.0*_b);
    
        pf=-_sigmaf*2.0*_b*pow(2.0*_Ni[_qp],_b-1);
        //pf=_sigmaf*_epsilonf*_E*(_c-_b)*pow(2.0*_Ni[_qp],_c-_b-1)*2.0+4.0*_Pswt2[_qp]*_b*pow(2.0*_Ni[_qp],-2.0*_b-1);
        
        _Ni[_qp]=_Ni[_qp]-f3/pf;
        
        if(_Ni[_qp] < 0.0 || _Ni[_qp] > 1e8)
         {
             _Ni[_qp]=1e8;
             break;
         }
//     cout << _Ni[_qp] <<" " << _sigmamax[_qp] <<" " << Pswt2 <<" " << test <<endl;
    }
   }
   
   else
   {
    while(fabs(f21) >0.0001)
    {
        
        f21=pow(_sigmaf,2)+_sigmaf*_epsilonf*_E*pow(2.0*_Ni[_qp],_c-_b)-Pswt21*pow(2.0*_Ni[_qp],-2.0*_b);
    
        pf=_sigmaf*_epsilonf*_E*(_c-_b)*pow(2.0*_Ni[_qp],_c-_b-1)*2.0+4.0*Pswt21*_b*pow(2.0*_Ni[_qp],-2.0*_b-1);
        
        _Ni[_qp]=_Ni[_qp]-f21/pf;
        
//          if(_Ni[_qp] < 0.0 || _Ni[_qp] > 1e8)
//          {
//              _Ni[_qp]=1e8;
//              break;
//          }
    }
   }
   
   
   //count error
    if(_Ni[_qp]<1e8)
    {
    if (_print_convergence_message)
      mooseWarning(
          "seiler2_Material:",
          _current_elem->id(),
          " and qp ",
          _Ni[_qp],',', _sigmareva[_qp], 
          "\n");
     }
   
   //cout << _Ni[_qp] <<endl;
   
//     if(_Ni[_qp] < 0.0 || _Ni[_qp] > 1e7)
//     {
//             _Ni[_qp]=1e7;
//           //  break;
//     }

    _D[_qp]=_D_old[_qp]+1.0/_Ni[_qp]*_cycle[_qp];

    
    if (_D[_qp]>1.0)
    {
        _D[_qp]=1.0;
    }
    
    if (_D[_qp]>_D0)
    {
        _alpha[_qp]=_alpha1;
    }
    else
    {
        _alpha[_qp]=(1.0-_alpha0)*pow(1.0-_D[_qp],_xi)+_alpha0;        
    }
    
    //_alpha[_qp]=(1.0-_alpha0)*pow(1.0-_D[_qp],_xi)+_alpha0;
    
     
}
