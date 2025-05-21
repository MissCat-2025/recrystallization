// ADGasSourceSink.h
#pragma once

#include "ADKernel.h"

/**
 * 核函数实现气体浓度方程中的源项: G ̇ - R ̇
 */
class ADGasSourceSink : public ADKernel
{
public:
  static InputParameters validParams();

  ADGasSourceSink(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  /// 裂变率
  const Real _fission_rate;
  
  /// 原子体积
  const Real _atomic_volume;
  
  /// 气体分辨率系数
  const Real _resolution_rate;
  
  /// 气泡相序参数
  unsigned int _num_bubbles;
  std::vector<const ADVariableValue *> _bubble_vars;
  
  /// 随机数种子
  unsigned int _seed;
  
  /// 是否在每个时间步中更新随机数
  const bool _update_seed_per_timestep;
  
  /// 存储随机数值
  std::vector<Real> _random_values;
};