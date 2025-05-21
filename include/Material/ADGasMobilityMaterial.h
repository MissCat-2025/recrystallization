// ... existing code ...
#pragma once

#include "ADMaterial.h"

class ADGasMobilityMaterial : public ADMaterial
{
public:
  static InputParameters validParams();

  ADGasMobilityMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  // 输入参数
  const Real _M_b;
  const Real _M_g;
  const unsigned int _num_etas;
  std::vector<const ADVariableValue *> _eta_vars;

  // 输出材料属性
  ADMaterialProperty<Real> & _M_eta;
};