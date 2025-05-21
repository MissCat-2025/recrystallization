// ADDislocEnergyDensity.h
#pragma once

#include "DerivativeMaterialInterface.h"
#include "Material.h"

/**
 * 材料类，用于计算位错引起的储存弹性能量密度
 * f_stored(η_i) = 1/2 η_i^2 G b_v^2 ρ_i, (i=1,2,...,p)
 * 其中位错密度计算为：ρ_d(t) = ρ_d(0)⋅exp[λḟ⋅t]
 */
class ADDislocEnergyDensity : public DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();
  ADDislocEnergyDensity(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// 计算位错密度
  Real computeDislocationDensity() const;

  /// 输入的剪切模量
  const Real _shear_modulus;
  
  /// 输入的伯格斯矢量
  const Real _burgers_vector;
  
  /// 初始位错密度
  const Real _rho_d_0;
  
  /// 位错累积速率
  const Real _lambda;
  
  /// 裂变速率
  const Real _fission_rate;
  
  /// 相场变量η的数量
  const unsigned int _num_etas;
  
  /// 相场变量值
  std::vector<const VariableValue *> _eta_values;
  
  /// 材料属性：位错密度
  MaterialProperty<Real> & _dislocation_density;
  
  /// 材料属性：储存弹性能量密度
  MaterialProperty<Real> & _f_stored;
  
  /// 储存弹性能量密度对各相场变量的导数
  std::vector<MaterialProperty<Real> *> _df_stored_deta;
  
  /// 储存弹性能量密度对相场变量的二阶导数
  std::vector<std::vector<MaterialProperty<Real> *>> _d2f_stored_deta;
  
  /// 时间
  const Real & _t;
};