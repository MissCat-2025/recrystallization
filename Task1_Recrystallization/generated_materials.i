# f_bulk_auto 材料项
[./f_bulk_auto]
  type = ADDerivativeParsedMaterial
  property_name = f_bulk_auto
  coupled_variables = 'c_g eta1 eta2 bubble_eta1 bubble_eta2'
  constant_names = 'A C_p C_q S_I c_m_e c_b_e'
  constant_expressions = '${A} ${C_p} ${C_q} ${S_I} ${c_m_e} ${c_b_e}'
  expression = 'A * ((-0.5*eta1^2 + 0.25*eta1^4 + -0.5*eta2^2 + 0.25*eta2^4 + -0.5*bubble_eta1^2 + 0.25*bubble_eta1^4 + -0.5*bubble_eta2^2 + 0.25*bubble_eta2^4 + 0.25) + (C_p * (eta1^2*eta2^2 + eta2^2*eta1^2)) + (C_q * (bubble_eta1^2*eta1^2 + bubble_eta1^2*eta2^2 + bubble_eta2^2*eta1^2 + bubble_eta2^2*eta2^2)) + (S_I * (c_g - c_m_e)^2) * ((eta1^2 + eta2^2) / (eta1^2 + eta2^2 + bubble_eta1^2 + bubble_eta2^2)) + (S_I * (c_g - c_b_e)^2) * (1 - ((eta1^2 + eta2^2) / (eta1^2 + eta2^2 + bubble_eta1^2 + bubble_eta2^2))))'
[../]

# f_stored_auto 材料项
[./f_stored_auto]
  type = ADDerivativeParsedMaterial
  property_name = f_stored_auto
  coupled_variables = 'eta1 eta2'
  material_property_names = 'rho'
  constant_names = 'G b_v2'
  constant_expressions = '${shear_modulus} ${burgers_vector}'
  expression = '0.5*G*b_v2*rho*(eta1^2 + eta2^2)'
[../]
# M_eta 材料项
[./M_eta]
  type = ADDerivativeParsedMaterial
  property_name = M_eta
  coupled_variables = 'eta1 eta2'
  constant_names = 'M_b M_g'
  constant_expressions = '${M_b} ${M_g}'
  expression = 'M_b*(eta1^3*(6*eta1^2-15*eta1+10)) + M_g*(1-(eta1^3*(6*eta1^2-15*eta1+10))) + M_b*(eta2^3*(6*eta2^2-15*eta2+10)) + M_g*(1-(eta2^3*(6*eta2^2-15*eta2+10)))'
  outputs = exodus
[../]
