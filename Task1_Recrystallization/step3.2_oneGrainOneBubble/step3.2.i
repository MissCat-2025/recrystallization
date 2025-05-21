################################################################################
## 模型step4：简化的耦合位错密度演化、气泡形成与晶粒生长
## mpirun -n 9 /home/yp/projects/recrystallization/recrystallization-opt -i step3.2.i
################################################################################

[Mesh]
    type        = GeneratedMesh
    dim         = 2
    nx          = 40
    ny          = 40
    xmin        = 0
    xmax        = 10
    ymin        = 0
    ymax        = 10
  []
  
  [GlobalParams]
    outputs     = exodus
    derivative_order = 2
  []
  
  #------------------------------------------------------------------------------#
  # 1. 辅助变量：位错密度, 气泡成核概率, 局部能量
  #------------------------------------------------------------------------------#
  [AuxVariables]
    [./dislocation_density]  # 位错密度 ρ_d(t)
      order      = CONSTANT
      family     = MONOMIAL
      [./InitialCondition]
        type     = ConstantIC
        value    = 5.0e10    # 降低初始位错密度
      [../]
    [../]
    [./stored_energy]  # 位错储存弹性能
      order      = CONSTANT
      family     = MONOMIAL
    [../]
    [./local_energy]  # 局部自由能
      order      = CONSTANT
      family     = MONOMIAL
    [../]
    [./bnds]   # 晶界指示
      order     = CONSTANT
      family    = MONOMIAL
    [../]
    [./bubble_sum_squared]   # 气泡相场平方和
      order     = CONSTANT
      family    = MONOMIAL
    [../]
    []
  
  [AuxKernels]
    # 计算位错密度演化 dρ_d(t)/dt = [λ·f_dot - g(T)]·ρ_d(t)
    [./dislocation_evolution]
      type      = ParsedAux
      variable  = dislocation_density
      coupled_variables = 'dislocation_density'
      expression = 'dislocation_density * exp(4.0e3 * 5.0e20 * 0.05)'  # 直接计算，无需使用函数
      execute_on = 'timestep_begin'
    [../]
    
    # 计算储存弹性能 f_stored(η_i) = 1/2·η_i^2·G·b_v^2·ρ_d
    [./stored_elastic_energy]
      type      = ParsedAux
      variable  = stored_energy
      coupled_variables = 'dislocation_density eta1'
      expression = '0.5 * 36.0e9 * (3.42e-10)^2 * dislocation_density * eta1^2'
    [../]
    
    # 计算局部自由能
    [./local_free_energy]
      type      = ParsedAux
      variable  = local_energy
      coupled_variables = 'eta1 bubble_eta1'
      expression = 'eta1^2 + bubble_eta1^2'  # 简化的自由能表达式
    [../]
    
    # 计算晶界
    [./bnds_auxiliary]
      type = BndsCalcAux
      variable = bnds
      execute_on = 'timestep_end'
      v = 'eta1'
    [../]
    
    [./bubble_sum_squared_aux]
      type = ParsedAux
      variable = bubble_sum_squared
      coupled_variables = 'bubble_eta1'
      expression = 'bubble_eta1^2'
      execute_on = 'timestep_begin timestep_end'
    [../]
    []
  
  #------------------------------------------------------------------------------#
  # 2. 主变量：晶粒相场、气泡相场、气体浓度
  #------------------------------------------------------------------------------#
  [Variables]
    # 晶粒相场变量
    [./eta1]
      order      = FIRST
      family     = LAGRANGE
      [./InitialCondition]
        type     = RandomIC
        seed     = 1
      [../]
    [../]
    
    # 气泡相场变量
    [./bubble_eta1]
      order      = FIRST
      family     = LAGRANGE
      [./InitialCondition]
        type     = ConstantIC
        value    = 0.0    # 初始无气泡
      [../]
    [../]
    
    # 气体浓度变量 c_g
    [./c_g]
      order      = FIRST
      family     = LAGRANGE
      [./InitialCondition]
        type     = ConstantIC
        value    = 0.0005  # 初始气体浓度很低
      [../]
    [../]
    []
  
  #------------------------------------------------------------------------------#
  # 3. Kernels：Allen-Cahn + Cahn-Hilliard
  #------------------------------------------------------------------------------#
  [Kernels]
    #-------------------#
    # 晶粒相场变量方程 #
    #-------------------#
    [./deta1dt]
      type      = TimeDerivative
      variable  = eta1
    [../]
    
    [./ACBulk1]
      type      = AllenCahn
      variable  = eta1
      f_name    = F_total
      mob_name  = L_eta
    [../]
    
    [./ACInterface1]
      type      = ACInterface
      variable  = eta1
      mob_name  = L_eta
      kappa_name = kappa_eta
    [../]
    
    [./elastic_energy1]
      type      = AllenCahn
      variable  = eta1
      f_name    = F_stored
      mob_name  = L_eta
    [../]
  
    #-------------------#
    # 气泡相场变量方程 #
    #-------------------#
    [./dbubble_etadt1]
      type      = TimeDerivative
      variable  = bubble_eta1
    [../]
    
    [./ACBulk_bubble1]
      type      = AllenCahn
      variable  = bubble_eta1
      f_name    = F_total
      mob_name  = L_bubble
    [../]
    
    [./ACInterface_bubble1]
      type      = ACInterface
      variable  = bubble_eta1
      mob_name  = L_bubble
      kappa_name = kappa_bubble
    [../]
    
    [./bubble_nucleation1]
      type      = LangevinNoise
      variable  = bubble_eta1
      amplitude = 1.0  # 降低核化强度
    [../]
    
    [./bubble1_grain_interaction]
      type      = AllenCahn
      variable  = bubble_eta1
      f_name    = F_grain_bubble_interaction
      mob_name  = L_bubble
      coupled_variables = 'eta1'
    [../]
    
    #-------------------#
    # 气体浓度变量方程 #
    #-------------------#
    [./dc_gdt]
      type      = TimeDerivative
      variable  = c_g
    [../]
    
    [./CHBulk]
      type      = CahnHilliard
      variable  = c_g
      f_name    = F_gas
      mob_name  = M_c_g
    [../]
    
    [./CHInterface]
      type      = CHInterface
      variable  = c_g
      mob_name  = M_c_g
      kappa_name = kappa_c
    [../]
    
    [./gas_generation]
      type      = BodyForce
      variable  = c_g
      value     = 0.005  # 降低气体产生率
    [../]
    
    [./gas_resolution]
      type      = MatReaction
      variable  = c_g
      mob_name  = resolution_rate
    [../]
    
    [./bubble1_concentration_coupling]
      type      = CoupledTimeDerivative
      variable  = bubble_eta1
      v         = c_g
    [../]
    []
  
  [BCs]
    [./Periodic]
      [./All]
        auto_direction = 'x y'
      [../]
    [../]
    []
  
  #------------------------------------------------------------------------------#
  # 4. Materials：常数、自由能、动力学参数
  #------------------------------------------------------------------------------#
  [Materials]
    # 物理常数 - 调整参数提高稳定性
    [./consts]
      type         = GenericConstantMaterial
      prop_names   = 'A Cq L_eta L_bubble M_bulk M_gb
                      kappa_eta kappa_bubble kappa_c
                      G b_v lambda fission_rate
                      c_m_e c_b_e S_I
                      resolution_rate'
      prop_values  = '1.0e7 1.0 5.0e-14 5.0e-14 2.5e-25 2.5e-23
                      3.75e-8 3.75e-8 2.74e-7
                      36.0e9 3.42e-10 4.0e3 5.0e20
                      1.0e-7 1.0 1.0
                      -0.0001'
    [../]
    
    # 气体浓度迁移率 M(η) = M_bulk * h(η) + M_gb * [1-h(η)]
    [./gas_mobility]
      type         = DerivativeParsedMaterial
      property_name = M_c_g
      coupled_variables = 'eta1'
      expression   = 'M_bulk * (eta1^3 * (6*eta1^2 - 15*eta1 + 10)) + 
                      M_gb * (1 - (eta1^3 * (6*eta1^2 - 15*eta1 + 10)))'
      material_property_names = 'M_bulk M_gb'
      constant_on  = 'ELEMENT'
      outputs      = exodus
    [../]
    
    # 晶粒切换函数 h_i
    [./switching1]
      type         = SwitchingFunctionMaterial
      function_name= h1
      eta          = eta1
      h_order      = SIMPLE
    [../]
    
    # 气泡切换函数
    [./switching_bubble1]
      type         = SwitchingFunctionMaterial
      function_name= h_bubble1
      eta          = bubble_eta1
      h_order      = SIMPLE
    [../]
    
    # 体积分数计算
    [./theta_m]
      type         = DerivativeParsedMaterial
      property_name = theta_m
      coupled_variables = 'eta1 bubble_eta1'
      expression   = 'eta1^2 / (eta1^2 + bubble_eta1^2 + 1e-10)'
      outputs      = exodus
    [../]
    
    [./theta_b]
      type         = DerivativeParsedMaterial
      property_name = theta_b
      coupled_variables = 'eta1 bubble_eta1'
      expression   = 'bubble_eta1^2 / (eta1^2 + bubble_eta1^2 + 1e-10)'
      outputs      = exodus
    [../]
    
    # 气体自由能 f_g^m = S_I * (c_g - c_m_e)^2, f_g^b = S_I * (c_g - c_b_e)^2
    [./gas_free_energy_m]
      type         = DerivativeParsedMaterial
      property_name = f_g_m
      coupled_variables = 'c_g'
      expression   = 'S_I * (c_g - c_m_e)^2'
      material_property_names = 'S_I c_m_e'
    [../]
    
    [./gas_free_energy_b]
      type         = DerivativeParsedMaterial
      property_name = f_g_b
      coupled_variables = 'c_g'
      expression   = 'S_I * (c_g - c_b_e)^2'
      material_property_names = 'S_I c_b_e'
    [../]
    
    # 位错储存能 f_stored(η_i) = 1/2·η_i^2·G·b_v^2·ρ_d
    [./stored_elastic_energy]
      type         = DerivativeParsedMaterial
      property_name = F_stored
      coupled_variables = 'dislocation_density eta1'
      expression   = '0.5 * G * b_v^2 * dislocation_density * eta1^2'
      material_property_names = 'G b_v'
      outputs      = exodus
    [../]
    
    # 相自由能 F_i：简单双井势
    [./phase_free_energy_1]
      type             = DerivativeParsedMaterial
      property_name    = F1
      coupled_variables= 'eta1'
      expression       = 'A*(-0.5*eta1^2 + 0.25*eta1^4)'
      material_property_names = 'A'
    [../]
    
    # 气泡自由能
    [./bubble_free_energy1]
      type             = DerivativeParsedMaterial
      property_name    = F_bubble1
      coupled_variables= 'bubble_eta1'
      expression       = 'A*(-0.5*bubble_eta1^2 + 0.25*bubble_eta1^4)'
      material_property_names = 'A'
    [../]
    
    # 晶粒-气泡相互作用 Cq*sum(i,j) η_i^2 * η_j^2
    [./grain_bubble_interaction]
      type             = DerivativeParsedMaterial
      property_name    = F_grain_bubble_interaction
      coupled_variables= 'eta1 bubble_eta1'
      expression       = 'Cq * bubble_eta1^2 * eta1^2'
      material_property_names = 'Cq'
    [../]
    
    # 气体化学自由能 f_g^m * θ_m + f_g^b * θ_b
    [./gas_chemical_energy]
      type             = DerivativeParsedMaterial
      property_name    = F_gas
      coupled_variables= 'c_g'
      expression       = 'f_g_m * theta_m + f_g_b * theta_b'
      material_property_names = 'f_g_m f_g_b theta_m theta_b'
      outputs          = exodus
    [../]
    
    # 总自由能 F_total = F_i + F_interaction + F_gas + 0.25
    [./total_free_energy]
      type             = DerivativeParsedMaterial
      property_name    = F_total
      coupled_variables= 'eta1 bubble_eta1 c_g'
      expression       = 'h1*F1 + h_bubble1*F_bubble1 + F_grain_bubble_interaction + F_gas + 0.25'
      material_property_names = 'h1 h_bubble1 F1 F_bubble1 F_grain_bubble_interaction F_gas'
      outputs          = exodus
    [../]
    []
  
  [Postprocessors]
    [./total_energy]
      type      = ElementIntegralVariablePostprocessor
      variable  = local_energy
    [../]
    [./bubble_volume_fraction]
      type      = ElementAverageValue
      variable  = bubble_eta1
    [../]
    [./average_gas_concentration]
      type      = ElementAverageValue
      variable  = c_g
    [../]
    [./average_dislocation_density]
      type      = ElementAverageValue
      variable  = dislocation_density
    [../]
    []
  
  [Preconditioning]
    [./SMP]
      type      = SMP
      full      = true
      petsc_options_iname = '-pc_type -sub_pc_type -pc_asm_overlap -ksp_gmres_restart'
      petsc_options_value = 'asm      lu           2               200'
    [../]
    []
  
  [Executioner]
    type            = Transient
    scheme          = bdf2
    solve_type      = PJFNK
    start_time      = 0.0
    end_time        = 100.0
    [./TimeStepper]
      type = SolutionTimeAdaptiveDT
      dt = 1.0e-5  # 使用非常小的初始时间步长
    [../]
    nl_rel_tol      = 1e-7
    nl_abs_tol      = 1e-9
    nl_max_its      = 25
    l_max_its       = 100
    l_tol           = 1e-6
    line_search     = 'bt'  # 使用回溯线搜索
    [./Adaptivity]
      initial_adaptivity = 2
      refine_fraction = 0.3
      coarsen_fraction = 0.05
    [../]
    automatic_scaling = true
    compute_scaling_once = false
    petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
    petsc_options_value = 'hypre boomeramg 200'
  []
  
  [Outputs]
    [./exodus]
      type = Exodus
      time_step_interval = 1
      show = 'eta1 bubble_eta1 c_g dislocation_density stored_energy bnds bubble_sum_squared'
    [../]
    print_linear_residuals = false
    checkpoint = true
  []