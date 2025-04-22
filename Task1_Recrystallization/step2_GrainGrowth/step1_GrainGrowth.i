################################################################################
## mpirun -n 12 /home/yp/projects/recrystallization/recrystallization-opt -i step1_GrainGrowth.i
# Polycrystalline Grain Growth Simulation
#   10 相多相场模型，初始各相随机分布
################################################################################

[Mesh]
    type        = GeneratedMesh
    dim         = 2
    nx          = 60
    ny          = 60
    xmin        = 0
    xmax        = 12
    ymin        = 0
    ymax        = 12
    elem_type   = Quad4
  []
  
  [GlobalParams]
    outputs     = exodus
    penalty     = 500    # SwitchingFunctionPenalty 强度
  []
  
  #------------------------------------------------------------------------------#
  # 1. 辅助变量：局部自由能及梯度交叉项
  #------------------------------------------------------------------------------#
  [AuxVariables]
    [./local_energy]
      order      = Constant
      family     = Monomial
    [../]
    [./cross_energy]
      order      = Constant
      family     = Monomial
    [../]
  []
  
  [AuxKernels]
    [./local_free_energy]
      type               = TotalFreeEnergy
      variable           = local_energy
      interfacial_vars   = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7 eta8 eta9 eta10'
      kappa_names        = 'kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta'
      additional_free_energy = cross_energy
    [../]
    [./cross_terms]
      type               = CrossTermGradientFreeEnergy
      variable           = cross_energy
      interfacial_vars   = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7 eta8 eta9 eta10'
      # 10×10 个项，全部用同一个 kappa_eta
      kappa_names        = 'kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta
                            kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta
                            kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta
                            kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta
                            kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta
                            kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta
                            kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta
                            kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta
                            kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta
                            kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta'
    [../]
  []
  
  #------------------------------------------------------------------------------#
  # 2. 主变量：10 个相参数 η1…η10，初始随机
  #------------------------------------------------------------------------------#
  [Variables]
    [./eta1]
      order      = First
      family     = Lagrange
      [./InitialCondition]
        type         = RandomIC
        seed         =  1
      [../]
    [../]
    [./eta2]
      order      = First
      family     = Lagrange
      [./InitialCondition]
        type         = RandomIC
        seed         =  2
      [../]
    [../]
    [./eta3]
      order      = First
      family     = Lagrange
      [./InitialCondition]
        type         = RandomIC
        seed         =  3
      [../]
    [../]
    [./eta4]
      order      = First
      family     = Lagrange
      [./InitialCondition]
        type         = RandomIC
        seed         =  4
      [../]
    [../]
    [./eta5]
      order      = First
      family     = Lagrange
      [./InitialCondition]
        type         = RandomIC
        seed         =  5
      [../]
    [../]
    [./eta6]
      order      = First
      family     = Lagrange
      [./InitialCondition]
        type         = RandomIC
        seed         =  6
    [../]
    [../]
    [./eta7]
      order      = First
      family     = Lagrange
      [./InitialCondition]
        type         = RandomIC
        seed         =  7
      [../]
    [../]
    [./eta8]
      order      = First
      family     = Lagrange
      [./InitialCondition]
        type         = RandomIC
        seed         =  8
      [../]
    [../]
    [./eta9]
      order      = First
      family     = Lagrange
      [./InitialCondition]
        type         = RandomIC
        seed         =  9
      [../]
    [../]
    [./eta10]
      order      = First
      family     = Lagrange
      [./InitialCondition]
        type         = RandomIC
        seed         = 10
      [../]
    [../]
  []
  
  #------------------------------------------------------------------------------#
  # 3. Kernels：纯 Allen‐Cahn + ACMultiInterface + SwitchingFunctionPenalty
  #------------------------------------------------------------------------------#
  [Kernels]
    # eta1
    [./deta1dt]
      type      = TimeDerivative
      variable  = eta1
    [../]
    [./ACBulk1]
      type               = AllenCahn
      variable           = eta1
      coupled_variables  = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7 eta8 eta9 eta10'
      mob_name           = L_eta
      f_name             = F
    [../]
    [./ACInterface1]
      type               = ACMultiInterface
      variable           = eta1
      etas               = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7 eta8 eta9 eta10'
      mob_name           = L_eta
      kappa_names        = 'kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta'
    [../]
    [./penalty1]
      type      = SwitchingFunctionPenalty
      variable  = eta1
      etas      = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7 eta8 eta9 eta10'
      h_names   = 'h1 h2 h3 h4 h5 h6 h7 h8 h9 h10'
    [../]
  
    # eta2
    [./deta2dt]
      type      = TimeDerivative
      variable  = eta2
    [../]
    [./ACBulk2]
      type               = AllenCahn
      variable           = eta2
      coupled_variables  = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7 eta8 eta9 eta10'
      mob_name           = L_eta
      f_name             = F
    [../]
    [./ACInterface2]
      type               = ACMultiInterface
      variable           = eta2
      etas               = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7 eta8 eta9 eta10'
      mob_name           = L_eta
      kappa_names        = 'kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta'
    [../]
    [./penalty2]
      type      = SwitchingFunctionPenalty
      variable  = eta2
      etas      = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7 eta8 eta9 eta10'
      h_names   = 'h1 h2 h3 h4 h5 h6 h7 h8 h9 h10'
    [../]
  
    # eta3
    [./deta3dt]
      type      = TimeDerivative
      variable  = eta3
    [../]
    [./ACBulk3]
      type               = AllenCahn
      variable           = eta3
      coupled_variables  = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7 eta8 eta9 eta10'
      mob_name           = L_eta
      f_name             = F
    [../]
    [./ACInterface3]
      type               = ACMultiInterface
      variable           = eta3
      etas               = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7 eta8 eta9 eta10'
      mob_name           = L_eta
      kappa_names        = 'kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta'
    [../]
    [./penalty3]
      type      = SwitchingFunctionPenalty
      variable  = eta3
      etas      = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7 eta8 eta9 eta10'
      h_names   = 'h1 h2 h3 h4 h5 h6 h7 h8 h9 h10'
    [../]
  
    # eta4
    [./deta4dt]
      type      = TimeDerivative
      variable  = eta4
    [../]
    [./ACBulk4]
      type               = AllenCahn
      variable           = eta4
      coupled_variables  = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7 eta8 eta9 eta10'
      mob_name           = L_eta
      f_name             = F
    [../]
    [./ACInterface4]
      type               = ACMultiInterface
      variable           = eta4
      etas               = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7 eta8 eta9 eta10'
      mob_name           = L_eta
      kappa_names        = 'kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta'
    [../]
    [./penalty4]
      type      = SwitchingFunctionPenalty
      variable  = eta4
      etas      = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7 eta8 eta9 eta10'
      h_names   = 'h1 h2 h3 h4 h5 h6 h7 h8 h9 h10'
    [../]
  
    # eta5
    [./deta5dt]
      type      = TimeDerivative
      variable  = eta5
    [../]
    [./ACBulk5]
      type               = AllenCahn
      variable           = eta5
      coupled_variables  = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7 eta8 eta9 eta10'
      mob_name           = L_eta
      f_name             = F
    [../]
    [./ACInterface5]
      type               = ACMultiInterface
      variable           = eta5
      etas               = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7 eta8 eta9 eta10'
      mob_name           = L_eta
      kappa_names        = 'kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta'
    [../]
    [./penalty5]
      type      = SwitchingFunctionPenalty
      variable  = eta5
      etas      = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7 eta8 eta9 eta10'
      h_names   = 'h1 h2 h3 h4 h5 h6 h7 h8 h9 h10'
    [../]
  
    # eta6
    [./deta6dt]
      type      = TimeDerivative
      variable  = eta6
    [../]
    [./ACBulk6]
      type               = AllenCahn
      variable           = eta6
      coupled_variables  = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7 eta8 eta9 eta10'
      mob_name           = L_eta
      f_name             = F
    [../]
    [./ACInterface6]
      type               = ACMultiInterface
      variable           = eta6
      etas               = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7 eta8 eta9 eta10'
      mob_name           = L_eta
      kappa_names        = 'kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta'
    [../]
    [./penalty6]
      type      = SwitchingFunctionPenalty
      variable  = eta6
      etas      = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7 eta8 eta9 eta10'
      h_names   = 'h1 h2 h3 h4 h5 h6 h7 h8 h9 h10'
    [../]
  
    # eta7
    [./deta7dt]
      type      = TimeDerivative
      variable  = eta7
    [../]
    [./ACBulk7]
      type               = AllenCahn
      variable           = eta7
      coupled_variables  = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7 eta8 eta9 eta10'
      mob_name           = L_eta
      f_name             = F
    [../]
    [./ACInterface7]
      type               = ACMultiInterface
      variable           = eta7
      etas               = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7 eta8 eta9 eta10'
      mob_name           = L_eta
      kappa_names        = 'kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta'
    [../]
    [./penalty7]
      type      = SwitchingFunctionPenalty
      variable  = eta7
      etas      = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7 eta8 eta9 eta10'
      h_names   = 'h1 h2 h3 h4 h5 h6 h7 h8 h9 h10'
    [../]
  
    # eta8
    [./deta8dt]
      type      = TimeDerivative
      variable  = eta8
    [../]
    [./ACBulk8]
      type               = AllenCahn
      variable           = eta8
      coupled_variables  = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7 eta8 eta9 eta10'
      mob_name           = L_eta
      f_name             = F
    [../]
    [./ACInterface8]
      type               = ACMultiInterface
      variable           = eta8
      etas               = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7 eta8 eta9 eta10'
      mob_name           = L_eta
      kappa_names        = 'kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta'
    [../]
    [./penalty8]
      type      = SwitchingFunctionPenalty
      variable  = eta8
      etas      = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7 eta8 eta9 eta10'
      h_names   = 'h1 h2 h3 h4 h5 h6 h7 h8 h9 h10'
    [../]
  
    # eta9
    [./deta9dt]
      type      = TimeDerivative
      variable  = eta9
    [../]
    [./ACBulk9]
      type               = AllenCahn
      variable           = eta9
      coupled_variables  = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7 eta8 eta9 eta10'
      mob_name           = L_eta
      f_name             = F
    [../]
    [./ACInterface9]
      type               = ACMultiInterface
      variable           = eta9
      etas               = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7 eta8 eta9 eta10'
      mob_name           = L_eta
      kappa_names        = 'kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta'
    [../]
    [./penalty9]
      type      = SwitchingFunctionPenalty
      variable  = eta9
      etas      = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7 eta8 eta9 eta10'
      h_names   = 'h1 h2 h3 h4 h5 h6 h7 h8 h9 h10'
    [../]
  
    # eta10
    [./deta10dt]
      type      = TimeDerivative
      variable  = eta10
    [../]
    [./ACBulk10]
      type               = AllenCahn
      variable           = eta10
      coupled_variables  = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7 eta8 eta9 eta10'
      mob_name           = L_eta
      f_name             = F
    [../]
    [./ACInterface10]
      type               = ACMultiInterface
      variable           = eta10
      etas               = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7 eta8 eta9 eta10'
      mob_name           = L_eta
      kappa_names        = 'kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta kappa_eta'
    [../]
    [./penalty10]
      type      = SwitchingFunctionPenalty
      variable  = eta10
      etas      = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7 eta8 eta9 eta10'
      h_names   = 'h1 h2 h3 h4 h5 h6 h7 h8 h9 h10'
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
  # 4. Materials：常数、切换函数、多相自由能
  #------------------------------------------------------------------------------#
  [Materials]
    [./consts]
      type         = GenericConstantMaterial
      prop_names   = 'A L_eta kappa_eta'
      prop_values  = '1.0 0.7 2.0'
    [../]
  
    # 阻挡函数 ∑ h_i = 1
    [./barrier]
      type         = MultiBarrierFunctionMaterial
      etas         = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7 eta8 eta9 eta10'
    [../]
  
    # 切换函数 h_i
    [./switching1]
      type         = SwitchingFunctionMaterial
      function_name= h1
      eta          = eta1
      h_order      = SIMPLE
    [../]
    [./switching2]
      type         = SwitchingFunctionMaterial
      function_name= h2
      eta          = eta2
      h_order      = SIMPLE
    [../]
    [./switching3]
      type         = SwitchingFunctionMaterial
      function_name= h3
      eta          = eta3
      h_order      = SIMPLE
    [../]
    [./switching4]
      type         = SwitchingFunctionMaterial
      function_name= h4
      eta          = eta4
      h_order      = SIMPLE
    [../]
    [./switching5]
      type         = SwitchingFunctionMaterial
      function_name= h5
      eta          = eta5
      h_order      = SIMPLE
    [../]
    [./switching6]
      type         = SwitchingFunctionMaterial
      function_name= h6
      eta          = eta6
      h_order      = SIMPLE
    [../]
    [./switching7]
      type         = SwitchingFunctionMaterial
      function_name= h7
      eta          = eta7
      h_order      = SIMPLE
    [../]
    [./switching8]
      type         = SwitchingFunctionMaterial
      function_name= h8
      eta          = eta8
      h_order      = SIMPLE
    [../]
    [./switching9]
      type         = SwitchingFunctionMaterial
      function_name= h9
      eta          = eta9
      h_order      = SIMPLE
    [../]
    [./switching10]
      type         = SwitchingFunctionMaterial
      function_name= h10
      eta          = eta10
      h_order      = SIMPLE
    [../]
  
    # 各相自由能 F_i：简单双井势
    [./phase_free_energy_1]
      type             = DerivativeParsedMaterial
      property_name    = F1
      expression       = 'A*(-0.5*eta1^2 + 0.25*eta1^4)'
      coupled_variables= 'eta1'
      material_property_names = 'A'
    [../]
    [./phase_free_energy_2]
      type             = DerivativeParsedMaterial
      property_name    = F2
      expression       = 'A*(-0.5*eta2^2 + 0.25*eta2^4)'
      coupled_variables= 'eta2'
      material_property_names = 'A'
    [../]
    [./phase_free_energy_3]
      type             = DerivativeParsedMaterial
      property_name    = F3
      expression       = 'A*(-0.5*eta3^2 + 0.25*eta3^4)'
      coupled_variables= 'eta3'
      material_property_names = 'A'
    [../]
    [./phase_free_energy_4]
      type             = DerivativeParsedMaterial
      property_name    = F4
      expression       = 'A*(-0.5*eta4^2 + 0.25*eta4^4)'
      coupled_variables= 'eta4'
      material_property_names = 'A'
    [../]
    [./phase_free_energy_5]
      type             = DerivativeParsedMaterial
      property_name    = F5
      expression       = 'A*(-0.5*eta5^2 + 0.25*eta5^4)'
      coupled_variables= 'eta5'
      material_property_names = 'A'
    [../]
    [./phase_free_energy_6]
      type             = DerivativeParsedMaterial
      property_name    = F6
      expression       = 'A*(-0.5*eta6^2 + 0.25*eta6^4)'
      coupled_variables= 'eta6'
      material_property_names = 'A'
    [../]
    [./phase_free_energy_7]
      type             = DerivativeParsedMaterial
      property_name    = F7
      expression       = 'A*(-0.5*eta7^2 + 0.25*eta7^4)'
      coupled_variables= 'eta7'
      material_property_names = 'A'
    [../]
    [./phase_free_energy_8]
      type             = DerivativeParsedMaterial
      property_name    = F8
      expression       = 'A*(-0.5*eta8^2 + 0.25*eta8^4)'
      coupled_variables= 'eta8'
      material_property_names = 'A'
    [../]
    [./phase_free_energy_9]
      type             = DerivativeParsedMaterial
      property_name    = F9
      expression       = 'A*(-0.5*eta9^2 + 0.25*eta9^4)'
      coupled_variables= 'eta9'
      material_property_names = 'A'
    [../]
    [./phase_free_energy_10]
      type             = DerivativeParsedMaterial
      property_name    = F10
      expression       = 'A*(-0.5*eta10^2 + 0.25*eta10^4)'
      coupled_variables= 'eta10'
      material_property_names = 'A'
    [../]
  
    # 多相自由能组合 F = Σ h_i F_i
    [./free_energy]
      type             = DerivativeMultiPhaseMaterial
      property_name    = F
      fi_names         = 'F1 F2 F3 F4 F5 F6 F7 F8 F9 F10'
      hi_names         = 'h1 h2 h3 h4 h5 h6 h7 h8 h9 h10'
      etas             = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7 eta8 eta9 eta10'
      W                = 2.0
    [../]
  []
  
  [Postprocessors]
    [./total_energy]
      type      = ElementIntegralVariablePostprocessor
      variable  = local_energy
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
    end_time        = 50.0
    [./TimeStepper]
      type = SolutionTimeAdaptiveDT
      dt = 0.00001
    [../]
    nl_rel_tol      = 1e-5
    nl_abs_tol      = 1e-7
    nl_max_its      = 15
    l_max_its       = 30
    l_tol           = 1e-4
    line_search     = 'basic'
    automatic_scaling = true
  []
  
  [Outputs]
    execute_on     = timestep_end
    exodus         = true
    perf_graph     = true
    print_linear_residuals = false
  []