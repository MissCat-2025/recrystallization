// ADBulkFreeEnergyDensity3.C
#include "ADBulkFreeEnergyDensity3.h"

registerMooseObject("RecrystallizationApp", ADBulkFreeEnergyDensity3);

InputParameters
ADBulkFreeEnergyDensity3::validParams()
{
  InputParameters params = Material::validParams();
  
  params.addRequiredCoupledVar("c_g", "气体浓度变量");
  params.addRequiredCoupledVar("grain_etas", "晶粒序参数列表，例如 'eta_1 eta_2 eta_3'");
  params.addRequiredCoupledVar("bubble_etas", "气泡序参数列表，例如 'eta_4 eta_5'");
  
  params.addRequiredParam<Real>("A", "自由能系数A");
  params.addRequiredParam<Real>("C_p", "多晶结构相互作用系数C_p");
  params.addRequiredParam<Real>("C_q", "气泡与多晶相互作用系数C_q");
  params.addRequiredParam<Real>("S_I", "气体自由能系数S_I");
  params.addRequiredParam<Real>("c_m_e", "颗粒结构中的Xe平衡浓度");
  params.addRequiredParam<Real>("c_b_e", "气泡中的Xe平衡浓度");
  
  params.addParam<MaterialPropertyName>("f_name", "f_bulk", "体自由能密度材料属性名称");
  
  params.addClassDescription("计算多晶与气泡相的体自由能密度材料，并计算对变量的导数");
  
  return params;
}

ADBulkFreeEnergyDensity3::ADBulkFreeEnergyDensity3(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _c_g(coupledValue("c_g")),
    _c_g_name(getVar("c_g", 0)->name()),
    _num_grains(coupledComponents("grain_etas")),
    _grain_vars(_num_grains),
    _grain_var_names(_num_grains),
    _num_bubbles(coupledComponents("bubble_etas")),
    _bubble_vars(_num_bubbles),
    _bubble_var_names(_num_bubbles),
    _A(getParam<Real>("A")),
    _C_p(getParam<Real>("C_p")),
    _C_q(getParam<Real>("C_q")),
    _S_I(getParam<Real>("S_I")),
    _c_m_e(getParam<Real>("c_m_e")),
    _c_b_e(getParam<Real>("c_b_e")),
    _f_bulk(declareProperty<Real>(getParam<MaterialPropertyName>("f_name")))
{
  // 存储晶粒序参数和变量名
  for (unsigned int i = 0; i < _num_grains; ++i)
  {
    _grain_vars[i] = &coupledValue("grain_etas", i);
    _grain_var_names[i] = getVar("grain_etas", i)->name();
  }
  
  // 存储气泡序参数和变量名
  for (unsigned int i = 0; i < _num_bubbles; ++i)
  {
    _bubble_vars[i] = &coupledValue("bubble_etas", i);
    _bubble_var_names[i] = getVar("bubble_etas", i)->name();
  }
}

void
ADBulkFreeEnergyDensity3::computeQpProperties()
{
  // 计算多晶结构相序参数平方和
  Real sum_grain_eta_sq = 0.0;
  for (unsigned int i = 0; i < _num_grains; ++i)
    sum_grain_eta_sq += (*_grain_vars[i])[_qp] * (*_grain_vars[i])[_qp];
  
  // 计算所有序参数平方和
  Real sum_all_eta_sq = sum_grain_eta_sq;
  for (unsigned int i = 0; i < _num_bubbles; ++i)
    sum_all_eta_sq += (*_bubble_vars[i])[_qp] * (*_bubble_vars[i])[_qp];
  
  // 计算θ_m和θ_b
  Real theta_m = sum_grain_eta_sq / sum_all_eta_sq;
  Real theta_b = 1.0 - theta_m;
  
  // 计算f_g^m和f_g^b
  Real c_g_diff_m = _c_g[_qp] - _c_m_e;
  Real c_g_diff_b = _c_g[_qp] - _c_b_e;
  Real f_g_m = _S_I * c_g_diff_m * c_g_diff_m;
  Real f_g_b = _S_I * c_g_diff_b * c_g_diff_b;
  
  // 计算各项总和
  Real sum_eta2_eta4 = 0.0;
  for (unsigned int i = 0; i < _num_grains; ++i)
    sum_eta2_eta4 += -0.5 * std::pow((*_grain_vars[i])[_qp], 2) + 0.25 * std::pow((*_grain_vars[i])[_qp], 4);
  
  for (unsigned int i = 0; i < _num_bubbles; ++i)
    sum_eta2_eta4 += -0.5 * std::pow((*_bubble_vars[i])[_qp], 2) + 0.25 * std::pow((*_bubble_vars[i])[_qp], 4);
  
  // 计算晶粒间相互作用
  Real sum_grain_interaction = 0.0;
  for (unsigned int i = 0; i < _num_grains; ++i)
    for (unsigned int j = 0; j < _num_grains; ++j)
      if (i != j)
        sum_grain_interaction += std::pow((*_grain_vars[i])[_qp], 2) * std::pow((*_grain_vars[j])[_qp], 2);
  
  // 计算气泡与晶粒间相互作用
  Real sum_bubble_grain_interaction = 0.0;
  for (unsigned int i = 0; i < _num_bubbles; ++i)
    for (unsigned int j = 0; j < _num_grains; ++j)
      sum_bubble_grain_interaction += std::pow((*_bubble_vars[i])[_qp], 2) * std::pow((*_grain_vars[j])[_qp], 2);
  
  // 计算体自由能密度
  _f_bulk[_qp] = _A * (sum_eta2_eta4 + _C_p * sum_grain_interaction + 
                      _C_q * sum_bubble_grain_interaction + 0.25 + 
                      f_g_m * theta_m + f_g_b * theta_b);
  
}