################################################################################
# Recrystallization + Gas Bubble + Dislocation Evolution Phase‐Field Model
################################################################################

[Mesh]
    type = GeneratedMesh
    dim = 2
    nx = 40
    ny = 40
    xmin = -12
    xmax = 12
    ymin = -12
    ymax = 12
    elem_type = QUAD4
  []
  
  [GlobalParams]
    outputs = exodus
    penalty = 1e3   # 切换函数罚项强度
  []
  
  #==============================================================================#
  # 1. 辅助变量：位错密度、局部与交叉自由能
  #==============================================================================#
  [AuxVariables]
    [./rho_d]
      order = CONSTANT
      family = MONOMIAL
    [../]
    [./local_energy]
      order = CONSTANT
      family = MONOMIAL
    [../]
    [./cross_energy]
      order = CONSTANT
      family = MONOMIAL
    [../]
  []
  
  #==============================================================================#
  # 2. AuxKernels：计算局部/交叉自由能，以及气体产生/消解源
  #==============================================================================#
  [AuxKernels]
    [./local_free_energy]
      type = TotalFreeEnergy
      variable = local_energy
      interfacial_vars = 'eta1 eta2 eta3 eta4'
      kappa_names = 'kappa_eta'
      additional_free_energy = cross_energy
    [../]
    [./cross_terms]
      type = CrossTermGradientFreeEnergy
      variable = cross_energy
      interfacial_vars = 'eta1 eta2 eta3 eta4'
      kappa_names = 'kappa11 kappa12 kappa13 kappa14
                     kappa21 kappa22 kappa23 kappa24
                     kappa31 kappa32 kappa33 kappa34
                     kappa41 kappa42 kappa43 kappa44'
    [../]
    # Cahn–Hilliard 的非保守源项：气体产生 \dot G 与消解 \dot R
    [./production]
      type = AnalyticSource
      variable = c_g
      function = 'varpi*rand()'
    [../]
    [./resolution]
      type = AnalyticSource
      variable = c_g
      function = '-Lambda*eta4^2*c_g'
    [../]
  []
  
  #==============================================================================#
  # 3. 主变量：气体浓度 c_g、晶粒相 eta1–3、气泡相 eta4
  #==============================================================================#
  [Variables]
    [./c_g]
      order = FIRST
      family = LAGRANGE
      [./InitialCondition]
        type = SmoothCircleIC
        x1 = 0.0; y1 = 0.0; radius = 5.0
        invalue = c_m_e; outvalue = c_b_e; int_width = 10.0
      [../]
    [../]
  
    # 三个晶粒相
    [./eta1]
      order = FIRST; family = LAGRANGE
      [./InitialCondition]
        type = FunctionIC
        function = 'r:=sqrt(x^2+y^2);if(r<=4,1,0)'
      [../]
    [../]
    [./eta2]
      order = FIRST; family = LAGRANGE
      [./InitialCondition]
        type = FunctionIC
        function = 'r:=sqrt(x^2+y^2);if(r>4&r<=7,1,0)'
      [../]
    [../]
    [./eta3]
      order = FIRST; family = LAGRANGE
      [./InitialCondition]
        type = FunctionIC
        function = 'r:=sqrt(x^2+y^2);if(r>7,1,0)'
      [../]
    [../]
  
    # 气泡相参数
    [./eta4]
      order = FIRST; family = LAGRANGE
      [./InitialCondition]
        type = RandomIC
        seed = 12345; initial_min = 0.0; initial_max = 0.01
      [../]
    [../]
  []
  
  #==============================================================================#
  # 4. Kernels：Cahn–Hilliard + TimeDerivative；Allen–Cahn + ACMultiInterface + 惩罚
  #==============================================================================#
  [Kernels]
    # ——— 气体浓度场 c_g ———
    [./c_res]
      type = CahnHilliard
      variable = c_g
      f_name = F
      coupled_variables = 'eta1 eta2 eta3 eta4'
    [../]
    [./time_c]
      type = TimeDerivative
      variable = c_g
    [../]
  
    # ——— 晶粒相 (i=1~3) ———
    [./deta1dt]   type = TimeDerivative; variable = eta1 [../]
    [./ACBulk1]
      type = AllenCahn; variable = eta1
      coupled_variables = 'eta2 eta3 eta4 c_g'
      mob_name = L_eta; f_name = F
    [../]
    [./ACInterface1]
      type = ACMultiInterface; variable = eta1
      etas = 'eta1 eta2 eta3 eta4'
      mob_name = L_eta
      kappa_names = 'kappa11 kappa12 kappa13 kappa14'
    [../]
    [./penalty1]
      type = SwitchingFunctionPenalty; variable = eta1
      etas = 'eta1 eta2 eta3 eta4'
      h_names = 'h1 h2 h3 h4'
    [../]
  
    [./deta2dt]   type = TimeDerivative; variable = eta2 [../]
    [./ACBulk2]
      type = AllenCahn; variable = eta2
      coupled_variables = 'eta1 eta3 eta4 c_g'
      mob_name = L_eta; f_name = F
    [../]
    [./ACInterface2]
      type = ACMultiInterface; variable = eta2
      etas = 'eta1 eta2 eta3 eta4'
      mob_name = L_eta
      kappa_names = 'kappa21 kappa22 kappa23 kappa24'
    [../]
    [./penalty2]
      type = SwitchingFunctionPenalty; variable = eta2
      etas = 'eta1 eta2 eta3 eta4'
      h_names = 'h1 h2 h3 h4'
    [../]
  
    [./deta3dt]   type = TimeDerivative; variable = eta3 [../]
    [./ACBulk3]
      type = AllenCahn; variable = eta3
      coupled_variables = 'eta1 eta2 eta4 c_g'
      mob_name = L_eta; f_name = F
    [../]
    [./ACInterface3]
      type = ACMultiInterface; variable = eta3
      etas = 'eta1 eta2 eta3 eta4'
      mob_name = L_eta
      kappa_names = 'kappa31 kappa32 kappa33 kappa34'
    [../]
    [./penalty3]
      type = SwitchingFunctionPenalty; variable = eta3
      etas = 'eta1 eta2 eta3 eta4'
      h_names = 'h1 h2 h3 h4'
    [../]
  
    # ——— 气泡相 eta4 ———
    [./deta4dt]   type = TimeDerivative; variable = eta4 [../]
    [./ACBulk4]
      type = AllenCahn; variable = eta4
      coupled_variables = 'eta1 eta2 eta3 c_g'
      mob_name = L_eta; f_name = F
    [../]
    [./ACInterface4]
      type = ACMultiInterface; variable = eta4
      etas = 'eta1 eta2 eta3 eta4'
      mob_name = L_eta
      kappa_names = 'kappa41 kappa42 kappa43 kappa44'
    [../]
    [./penalty4]
      type = SwitchingFunctionPenalty; variable = eta4
      etas = 'eta1 eta2 eta3 eta4'
      h_names = 'h1 h2 h3 h4'
    [../]
  []
  
  [BCs]
    [./Periodic]
      [./All]
        auto_direction = 'x y'
      [../]
    [../]
  []
  
  #==============================================================================#
  # 5. Materials：常数属性、位错密度解析、切换函数、势垒、相自由能与多相组合
  #==============================================================================#
  [Materials]
    # — 常数属性 — M_b, M_g, kappa_eta, kappa_ij, L_eta, varpi, Lambda, A, Cp, Cq,
    #             Sl, G, bv, c_m_e, c_b_e, rho_d0, lambda_dotf
    [./consts]
      type = GenericConstantMaterial
      prop_names  = 'M_b M_g kappa_eta kappa11 kappa12 kappa13 kappa14
                     kappa21 kappa22 kappa23 kappa24
                     kappa31 kappa32 kappa33 kappa34
                     kappa41 kappa42 kappa43 kappa44
                     L_eta varpi Lambda A Cp Cq Sl G bv
                     c_m_e c_b_e rho_d0 lambda_dotf'
      prop_values = '0.1 1.0 1.0  ... (请替换为实际值) ... 0.01 0.1'
    [../]
  
    # — 位错密度解析表达（忽略 g(T)）—
    [./rho_d_mat]
      type = ParsedMaterial
      property_name = rho_d
      expression = 'rho_d0*exp(lambda_dotf*t)'
      coupled_variables = 't'
    [../]
  
    # — 切换函数 h_i —
    [./switching1] type = SwitchingFunctionMaterial; function_name = h1; eta = eta1; h_order = HARD [../]
    [./switching2] type = SwitchingFunctionMaterial; function_name = h2; eta = eta2; h_order = HARD [../]
    [./switching3] type = SwitchingFunctionMaterial; function_name = h3; eta = eta3; h_order = HARD [../]
    [./switching4] type = SwitchingFunctionMaterial; function_name = h4; eta = eta4; h_order = HARD [../]
  
    # — 势垒函数 ∑h_i=1 —
    [./barrier]
      type = MultiBarrierFunctionMaterial
      etas = 'eta1 eta2 eta3 eta4'
    [../]
  
    # — 各相自由能 F_i —
    [./phase_free_energy_1]
      type = DerivativeParsedMaterial
      property_name = F1
      expression = 'A*(-0.5*eta1^2+0.25*eta1^4)
                    +A*Cp*(eta1^2*eta2^2+eta1^2*eta3^2+eta2^2*eta3^2)
                    +0.5*G*bv^2*rho_d*eta1^2
                    +Sl*(c_g-c_m_e)^2'
      coupled_variables = 'eta1 eta2 eta3 c_g rho_d'
    [../]
    [./phase_free_energy_2]
      type = DerivativeParsedMaterial
      property_name = F2
      expression = 'A*(-0.5*eta2^2+0.25*eta2^4)
                    +A*Cp*(eta1^2*eta2^2+eta2^2*eta3^2+eta1^2*eta3^2)
                    +0.5*G*bv^2*rho_d*eta2^2
                    +Sl*(c_g-c_m_e)^2'
      coupled_variables = 'eta1 eta2 eta3 c_g rho_d'
    [../]
    [./phase_free_energy_3]
      type = DerivativeParsedMaterial
      property_name = F3
      expression = 'A*(-0.5*eta3^2+0.25*eta3^4)
                    +A*Cp*(eta1^2*eta3^2+eta2^2*eta3^2+eta1^2*eta2^2)
                    +0.5*G*bv^2*rho_d*eta3^2
                    +Sl*(c_g-c_m_e)^2'
      coupled_variables = 'eta1 eta2 eta3 c_g rho_d'
    [../]
    [./phase_free_energy_4]
      type = DerivativeParsedMaterial
      property_name = F4
      expression = 'A*(-0.5*eta4^2+0.25*eta4^4)
                    +Cq*(eta4^2*(eta1^2+eta2^2+eta3^2))
                    +0.25
                    +Sl*(c_g-c_b_e)^2'
      coupled_variables = 'eta1 eta2 eta3 eta4 c_g'
    [../]
  
    # — 多相自由能组合 F = ∑ h_i F_i —
    [./free_energy]
      type = DerivativeMultiPhaseMaterial
      property_name = F
      fi_names = 'F1 F2 F3 F4'
      hi_names = 'h1 h2 h3 h4'
      etas     = 'eta1 eta2 eta3 eta4'
      coupled_variables = 'c_g'
      W = 1
    [../]
  []
  
  [Postprocessors]
    [./total_free_energy]
      type = ElementIntegralVariablePostprocessor; variable = local_energy
    [../]
    [./total_solute]
      type = ElementIntegralVariablePostprocessor; variable = c_g
    [../]
    [./total_dislocation]
      type = ElementIntegralVariablePostprocessor; variable = rho_d
    [../]
  []
  
  [Preconditioning]
    [./SMP]
      type = SMP; full = true
    [../]
  []
  
  [Executioner]
    type = Transient
    scheme = bdf2
    solve_type = NEWTON
    nl_max_its = 50; nl_rel_tol = 1e-6
    l_max_its = 15;   l_tol = 1e-6
    start_time = 0.0; end_time = 150.0
    [./TimeStepper]
      type = SolutionTimeAdaptiveDT; dt = 0.1
    [../]
  []
  
  [Outputs]
    execute_on = timestep_end
    exodus = true
    [./csv]
      type = CSV; delimiter = ' '
    [../]
  []