// ... existing code ...
#include "ADGasMobilityMaterial.h"

registerMooseObject("RecrystallizationApp", ADGasMobilityMaterial);

InputParameters
ADGasMobilityMaterial::validParams()
{
  InputParameters params = ADMaterial::validParams();
  params.addRequiredParam<Real>("M_b", "气体在基体中的迁移率");
  params.addRequiredParam<Real>("M_g", "气体在晶界上的迁移率");
  params.addRequiredCoupledVar("etas", "所有晶粒序参数变量名列表");
  params.addClassDescription("根据M(η)=M_b*h(η)+M_g*[1-h(η)]自动适配任意个etas的气体迁移率材料");
  params.addParam<std::string>("property_name", "M_eta", "输出材料属性名");
  return params;
}

ADGasMobilityMaterial::ADGasMobilityMaterial(const InputParameters & parameters)
  : ADMaterial(parameters),
    _M_b(getParam<Real>("M_b")),
    _M_g(getParam<Real>("M_g")),
    _num_etas(coupledComponents("etas")),
    _eta_vars(_num_etas),
    _M_eta(declareADProperty<Real>(getParam<std::string>("property_name")))
{
  for (unsigned int i = 0; i < _num_etas; ++i)
    _eta_vars[i] = &adCoupledValue("etas", i);
}

void
ADGasMobilityMaterial::computeQpProperties()
{
  // 计算h(eta)的平均值
  ADReal h_sum = 0.0;
  for (unsigned int i = 0; i < _num_etas; ++i)
  {
    const ADReal & eta = (*_eta_vars[i])[_qp];
    h_sum += eta * eta * eta * (6.0 * eta * eta - 15.0 * eta + 10.0);
  }
  ADReal h_eta = h_sum / _num_etas;

  // 计算M(eta)
  _M_eta[_qp] = _M_b * h_eta + _M_g * (1.0 - h_eta);
}