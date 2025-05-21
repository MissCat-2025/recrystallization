import os

def parse_etas(etas_str):
    etas = etas_str.strip().split()
    # 假设气泡序参以bubble_开头，其他为普通序参，cg为浓度
    normal_etas = [e for e in etas if not e.startswith('bubble_') and e != 'c_g']
    bubble_etas = [e for e in etas if e.startswith('bubble_')]
    cg = [e for e in etas if e == 'c_g']
    return normal_etas, bubble_etas, cg

def gen_f1(normal_etas, bubble_etas):
    terms = [f'-0.5*{e}^2 + 0.25*{e}^4' for e in normal_etas + bubble_etas]
    f1 = ' + '.join(terms) + ' + 0.25'
    return f1

def gen_f2(normal_etas):
    terms = []
    for i, ei in enumerate(normal_etas):
        for j, ej in enumerate(normal_etas):
            if i != j:
                terms.append(f'{ei}^2*{ej}^2')
    if terms:
        f2 = 'C_p * (' + ' + '.join(terms) + ')'
    else:
        f2 = '0'
    return f2

def gen_f3(normal_etas, bubble_etas):
    terms = []
    for ei in bubble_etas:
        for ej in normal_etas:
            terms.append(f'{ei}^2*{ej}^2')
    if terms:
        f3 = 'C_q * (' + ' + '.join(terms) + ')'
    else:
        f3 = '0'
    return f3

def gen_thetam(normal_etas, bubble_etas):
    # (sum normal_etas^2) / (sum all_etas^2)
    num = ' + '.join([f'{e}^2' for e in normal_etas]) if normal_etas else '0'
    denom = ' + '.join([f'{e}^2' for e in normal_etas + bubble_etas]) if (normal_etas + bubble_etas) else '1'
    return f'(({num}) / ({denom}))'

def gen_thetab(normal_etas, bubble_etas):
    # 1 - thetam
    thetam = gen_thetam(normal_etas, bubble_etas)
    return f'(1 - {thetam})'

def gen_fg_m():
    # S_I * (c_g - c_m_e)^2
    return 'S_I * (c_g - c_m_e)^2'

def gen_fg_b():
    # S_I * (c_g - c_b_e)^2
    return 'S_I * (c_g - c_b_e)^2'

def gen_f_bulk(f1, f2, f3, thetam, thetab):
    # A * (f1 + f2 + f3 + fg_m*thetam + fg_b*thetab)
    fg_m = 'S_I * (c_g - c_m_e)^2'
    fg_b = 'S_I * (c_g - c_b_e)^2'
    return f"A * (({f1}) + ({f2}) + ({f3}) + ({fg_m}) * {thetam} + ({fg_b}) * {thetab})"

def gen_f_stored(normal_etas):
    """
    合并所有eta_i的f_stored项，生成一个材料块
    f_stored_total = 0.5 * G * b_v2 * sum(eta_i^2 * rho)
    这里rho作为material_property_names输入
    """
    if not normal_etas:
        return '', []
    expr = ' + '.join([f"{eta}^2" for eta in normal_etas])
    total_expr = f"0.5*G*b_v2*rho*({expr})"
    return total_expr, normal_etas

def gen_M_eta(normal_etas):
    """
    合并所有eta_i的M(eta_i)项，生成一个材料块
    M_eta = sum(M_b * h(eta_i) + M_g * (1 - h(eta_i)))
    h(eta) = eta^3 * (6*eta^2 - 15*eta + 10)
    """
    if not normal_etas:
        return '', []
    terms = []
    for eta in normal_etas:
        h_expr = f"{eta}^3*(6*{eta}^2-15*{eta}+10)"
        terms.append(f"M_b*({h_expr}) + M_g*(1-({h_expr}))")
    expr = ' + '.join(terms)
    return expr, normal_etas

def generate_materials(etas_str):
    normal_etas, bubble_etas, cg = parse_etas(etas_str)
    f1 = gen_f1(normal_etas, bubble_etas)
    f2 = gen_f2(normal_etas)
    f3 = gen_f3(normal_etas, bubble_etas)
    thetam = gen_thetam(normal_etas, bubble_etas)
    thetab = gen_thetab(normal_etas, bubble_etas)
    f_bulk = gen_f_bulk(f1, f2, f3, thetam, thetab)

    # 输出文件放在脚本同目录
    out_path = os.path.join(os.path.dirname(__file__), 'generated_materials.i')
    with open(out_path, 'w', encoding='utf-8') as f:
        # 1. 先写f_bulk_auto
        f.write('# f_bulk_auto 材料项\n')
        f.write('[./f_bulk_auto]\n')
        f.write('  type = ADDerivativeParsedMaterial\n')
        f.write('  property_name = f_bulk_auto\n')
        args = ' '.join(['c_g'] + normal_etas + bubble_etas)
        f.write(f"  coupled_variables = '{args}'\n")
        f.write("  constant_names = 'A C_p C_q S_I c_m_e c_b_e'\n")
        f.write("  constant_expressions = '${A} ${C_p} ${C_q} ${S_I} ${c_m_e} ${c_b_e}'\n")
        f.write(f"  expression = '{f_bulk}'\n")
        f.write('[../]\n\n')

        # 2. 再写f_stored_auto
        f_stored_expr, etas = gen_f_stored(normal_etas)
        if f_stored_expr:
            f.write('# f_stored_auto 材料项\n')
            f.write('[./f_stored_auto]\n')
            f.write('  type = ADDerivativeParsedMaterial\n')
            f.write('  property_name = f_stored_auto\n')
            f.write(f"  coupled_variables = '{' '.join(etas)}'\n")
            f.write(f"  material_property_names = 'rho'\n")
            f.write("  constant_names = 'G b_v2'\n")
            f.write("  constant_expressions = '${shear_modulus} ${burgers_vector}'\n")
            f.write(f"  expression = '{f_stored_expr}'\n")
            f.write('[../]\n')

        # 3. 再写M_eta材料块（合并所有eta）
        M_expr, etas = gen_M_eta(normal_etas)
        if M_expr:
            f.write('# M_eta 材料项\n')
            f.write('[./M_eta]\n')
            f.write('  type = ADDerivativeParsedMaterial\n')
            f.write('  property_name = M_eta\n')
            f.write(f"  coupled_variables = '{' '.join(etas)}'\n")
            f.write("  constant_names = 'M_b M_g'\n")
            f.write("  constant_expressions = '${M_b} ${M_g}'\n")
            f.write(f"  expression = '{M_expr}'\n")
            f.write('  outputs = exodus\n')
            f.write('[../]\n')
    print(f'已生成 {out_path}，可直接复制到MOOSE输入文件。')

if __name__ == '__main__':
    # 直接运行脚本时，自动生成一个示例材料文件
    generate_materials('eta1 eta2 bubble_eta1 bubble_eta2 c_g') 