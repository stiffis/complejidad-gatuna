#!/usr/bin/env python3
import math
import sys
from itertools import combinations, permutations

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import golden, linprog

PHI = (1 + math.sqrt(5)) / 2
SEED = 20261105
fallos = []


def chequear(nombre, obtenido, esperado, tol=1e-6):
    ok = abs(obtenido - esperado) <= tol
    print(f"{'ok ' if ok else 'MAL'} {nombre}: {obtenido:.6g} (esperado {esperado:.6g}, tol {tol:g})")
    if not ok:
        fallos.append(nombre)


def seccion(titulo):
    print(f"\n== {titulo} ==")


def codigo():
    seccion("Código del maullido")
    p = np.array([0.5, 0.25, 0.125, 0.125])
    largos = np.array([1, 2, 3, 3])
    chequear("entropía en bits", float(-(p * np.log2(p)).sum()), 1.75)
    chequear("longitud esperada", float((p * largos).sum()), 1.75)
    chequear("desigualdad de Kraft", float((2.0 ** -largos).sum()), 1.0)
    k, n, c = 4, 10, 2
    comprimibles = (k ** (n - c) - 1) // (k - 1)
    print(f"    antojos con menos de {n - c} actos: {comprimibles} de {k ** n} = {comprimibles / k ** n:.5f}")
    chequear("cota de incompresibilidad", comprimibles / k ** n, k ** (-c) / (k - 1), 1e-5)


def caida_integrada(curva, periodo, puntos=4000):
    def campo(t, y):
        psi_p, lam = curva(t)
        i1, i2 = 1 - lam, lam
        return [i2 / (i1 + i2) * psi_p, -i1 / (i1 + i2) * psi_p]

    malla = np.linspace(0, periodo, puntos)
    return solve_ivp(campo, (0, periodo), [0, 0], t_eval=malla, rtol=1e-10, atol=1e-12)


def caida():
    seccion("Caída: holonomía igual a área")
    dl, dp = 0.6, 4.0
    lam0 = 0.5
    w = 2 * math.pi

    def elipse(t):
        psi_p = -(dp / 2) * w * math.sin(w * t)
        lam = lam0 + (dl / 2) * math.sin(w * t)
        return psi_p, lam

    sol = caida_integrada(elipse, 1.0)
    area = math.pi * (dp / 2) * (dl / 2)
    print(f"    área de la elipse: {area:.6f}")
    chequear("rotación del cuerpo 1", abs(sol.y[0][-1]), area, 1e-6)
    chequear("rotación del cuerpo 2", abs(sol.y[1][-1]), area, 1e-6)

    i1, i2, om = sp.symbols("I1 I2 omega", positive=True)
    lam_e = i2 / (i1 + i2)
    momento = sp.simplify(i1 * lam_e * om + i2 * (lam_e - 1) * om)
    chequear("momento angular de la conexión", float(momento), 0.0)

    k = 4
    chequear("salto de lambda con k=4", (k - 1) / (k + 1), 0.6)
    dlam = (k - 1) / (k + 1)
    print(f"    torsión mínima rectángulo k=4: {math.degrees(math.pi / dlam):.1f} grados")
    print(f"    torsión mínima elipse k=4: {math.degrees(4 / dlam):.1f} grados")
    chequear("cociente de energía cuadrado/círculo", (4 * math.sqrt(math.pi)) ** 2 / (2 * math.pi) ** 2, 4 / math.pi)
    print(f"    exceso de energía del cuadrado: {100 * (4 / math.pi - 1):.2f} %")
    print(f"    exceso de perímetro del cuadrado: {100 * (2 / math.sqrt(math.pi) - 1):.2f} %")

    eps = (dlam / 2) ** 2
    g = 9.81
    chequear("epsilon para k=4", eps, 0.09)
    h = 0.30
    t_caida = math.sqrt(2 * h / g)
    v = 2 * math.pi / (math.sqrt(eps) * t_caida)
    print(f"    T(h=30 cm) = {t_caida:.4f} s, v necesaria = {v:.1f} rad/s")
    chequear("h_min recupera 30 cm", 2 * math.pi ** 2 * g / (eps * v ** 2), h, 1e-9)
    psi_max = math.pi / 2
    n = math.ceil(4 / (dlam * psi_max))
    print(f"    ciclos con torsión de 90 grados: n = {n}, frecuencia = {n / t_caida:.1f} Hz")


def persuasion():
    seccion("Persuasión")
    mu0, mus = 0.3, 0.5
    q = mu0 * (1 - mus) / ((1 - mu0) * mus)
    chequear("q*", q, 3 / 7)
    chequear("posterior tras maullar", mu0 / (mu0 + (1 - mu0) * q), mus)
    chequear("probabilidad de comer", mu0 + (1 - mu0) * q, mu0 / mus)

    mejor = 0.0
    for qh in np.linspace(0, 1, 401):
        for qs in np.linspace(0, 1, 401):
            pm = mu0 * qh + (1 - mu0) * qs
            valor = 0.0
            if pm > 0 and mu0 * qh / pm >= mus - 1e-12:
                valor += pm
            pc = 1 - pm
            if pc > 0 and mu0 * (1 - qh) / pc >= mus - 1e-12:
                valor += pc
            mejor = max(mejor, valor)
    chequear("óptimo por búsqueda exhaustiva", mejor, mu0 / mus, 2e-3)

    rng = np.random.default_rng(SEED)
    n, p, rep = 6, mu0 / mus, 2_000_000
    sorteos = rng.random((rep, n)) < p
    chequear("cenas sin pizarra", float(sorteos.sum(axis=1).mean()), n * p, 3e-3)
    chequear("cenas con pizarra", float(sorteos.any(axis=1).mean()), 1 - (1 - p) ** n, 3e-3)
    print(f"    sin pizarra {n * p:.3f}, con pizarra {1 - (1 - p) ** n:.4f}")


def puerta():
    seccion("Puerta: opción real")
    b = sp.symbols("beta", positive=True)
    raiz = sp.solve(sp.Eq(b * (b - 1) - 1, 0), b)[0]
    chequear("beta dorada", float(raiz), PHI, 1e-12)
    chequear("umbral dorado", PHI / (PHI - 1), PHI ** 2, 1e-12)
    chequear("nunca sale con x0=I", 1 - PHI ** -2, 1 / PHI, 1e-12)
    print(f"    P(sale)={PHI ** -2:.4f}  P(nunca)={1 / PHI:.4f}")

    v_i = (PHI ** 2 - 1) * PHI ** (-2 * PHI)
    chequear("V(I)/I", v_i, PHI ** (1 - 2 * PHI), 1e-12)
    print(f"    E[e^-rt]={PHI ** (-2 * PHI):.4f}  V(I)/I={v_i:.4f}  costo humano/(c/r)={1 - PHI ** (-2 * PHI):.4f}")

    x, a, r, s2, mu, c = sp.symbols("x A r sigma2 mu c", positive=True)
    beta = sp.symbols("beta", positive=True)
    w = a * x ** beta - c / r
    edo = sp.Rational(1, 2) * s2 * x ** 2 * sp.diff(w, x, 2) + mu * x * sp.diff(w, x) - r * w - c
    edo_car = sp.simplify(edo.subs(a, 1) / x ** beta)
    chequear("EDO con raíz característica", float(edo_car.subs(mu, 0).subs(s2, 0.04).subs(r, 0.02).subs(beta, PHI).subs(c, 0.3)), 0.0, 1e-10)

    rr, sig2, mm = 0.02, 0.04, 0.0
    inv_n, cc = 1.0, 0.005
    bb = 0.5 - mm / sig2 + math.sqrt((mm / sig2 - 0.5) ** 2 + 2 * rr / sig2)
    chequear("beta general", bb, PHI, 1e-12)
    xc = bb / (bb - 1) * inv_n
    xp = bb / (bb - 1) * (inv_n - cc / rr)
    a_w = xp ** (1 - bb) / bb
    chequear("emparejamiento de valor social", a_w * xp ** bb - cc / rr, xp - inv_n, 1e-12)
    chequear("brecha gato menos planificador", xc - xp, bb / (bb - 1) * cc / rr, 1e-12)
    print(f"    x_gato={xc:.4f}  x_social={xp:.4f}  x_cierre={(cc / (rr * a_w)) ** (1 / bb):.4f}")

    grilla = np.linspace(0.05, 6, 4000)
    a_v = (xc - inv_n) / xc ** bb
    v = np.where(grilla < xc, a_v * grilla ** bb, grilla - inv_n)
    chequear("dominancia de V sobre el pago", float((v - np.maximum(grilla - inv_n, 0)).min()), 0.0, 1e-9)
    chequear("pegado suave en x*", a_v * bb * xc ** (bb - 1), 1.0, 1e-12)

    def valor_umbral(u):
        return -(u - inv_n) * (inv_n / u) ** bb

    chequear("umbral óptimo numérico", float(golden(valor_umbral, brack=(1.5, 2.6, 3.5))), xc, 1e-5)

    def valor_social(u):
        return -((inv_n / u) ** bb * (u - inv_n + cc / rr) - cc / rr)

    chequear("umbral social numérico", float(golden(valor_social, brack=(1.2, 1.9, 3.0))), xp, 1e-5)

    rng = np.random.default_rng(SEED)
    dt, pasos, rep = 0.5, 6000, 20000
    nu = mm - sig2 / 2
    meta = math.log(xc / inv_n)
    posicion = np.zeros(rep)
    alcanzo = np.zeros(rep, dtype=bool)
    for _ in range(pasos):
        activos = ~alcanzo
        cuantos = int(activos.sum())
        if cuantos == 0:
            break
        nueva = posicion[activos] + nu * dt + math.sqrt(sig2 * dt) * rng.standard_normal(cuantos)
        puente = np.exp(-2 * np.maximum(meta - posicion[activos], 0) * np.maximum(meta - nueva, 0) / (sig2 * dt))
        cruza = (nueva >= meta) | (rng.random(cuantos) < puente)
        indices = np.flatnonzero(activos)
        alcanzo[indices[cruza]] = True
        posicion[indices] = nueva
    chequear("Monte Carlo de P(sale)", float(alcanzo.mean()), PHI ** -2, 0.01)
    print(f"    fracción que sale: {alcanzo.mean():.4f} en {rep} gatos")


def sol():
    seccion("Parche de sol")
    v, w, r, c = 1e-4, 0.6, 30.0, 2000.0
    t_opt = (3 * c * w ** 2 / (2 * r * v ** 2)) ** (1 / 3)

    def tasa(t):
        return r * (1 - (v * t) ** 2 / (3 * w ** 2)) - c / t

    chequear("T* numérico contra cerrado", float(golden(lambda t: -tasa(t), brack=(500.0, 1500.0, 4000.0))), t_opt, 1e-2)
    chequear("teorema del valor marginal", r * (1 - (v * t_opt / w) ** 2), tasa(t_opt), 1e-9)
    print(f"    T*={t_opt:.1f} s = {t_opt / 60:.1f} min, desplazamiento={v * t_opt * 100:.1f} cm, g(T*)={tasa(t_opt):.3f} W")

    d = np.linspace(-v * t_opt, 0, 2001)
    medias = [np.mean(1 - (2 * (d0 + v * t_opt * np.linspace(0, 1, 2001)) / w) ** 2) for d0 in d]
    chequear("inicio simétrico óptimo", float(d[int(np.argmax(medias))]), -v * t_opt / 2, 2e-5)


def costo_noche(j, t, b):
    return t if t < j else j - 1 + b


def noche():
    seccion("Noche: ski rental")
    for b in (2, 3, 5, 10, 20, 50):
        c = 1 / (1 - (1 - 1 / b) ** b)
        tope = 6 * b
        peor = min(max(costo_noche(j, t, b) / min(t, b) for t in range(1, tope)) for j in range(1, tope))
        chequear(f"determinista B={b}", peor, 2 - 1 / b, 1e-12)
        adaptativo = min((j - 1 + b) / min(j, b) for j in range(1, 10 * b))
        chequear(f"adaptativo B={b}", adaptativo, 2 - 1 / b, 1e-12)

        p = [(c / b) * (1 - 1 / b) ** (b - j) for j in range(1, b + 1)]
        chequear(f"suma del dado B={b}", sum(p), 1.0, 1e-12)
        razones = [sum(pj * costo_noche(j, t, b) for j, pj in enumerate(p, 1)) / min(t, b) for t in range(1, 4 * b)]
        chequear(f"igualación primal B={b}", max(razones), c, 1e-10)
        chequear(f"igualación primal mínima B={b}", min(razones), c, 1e-10)

        cola = [(1 - 1 / b) ** (j - 1) for j in range(1, b + 1)]
        q = {t: cola[t - 1] / b for t in range(1, b)}
        m = max(2 * b, math.ceil(b / cola[b - 1]))
        q[m] = cola[b - 1]
        chequear(f"masa de la distribución dual B={b}", sum(q.values()), 1.0, 1e-12)
        opt_medio = sum(qt * min(t, b) for t, qt in q.items())
        minimo = min(sum(qt * costo_noche(j, t, b) for t, qt in q.items()) for j in range(1, 3 * m))
        chequear(f"cota dual B={b}", minimo / opt_medio, c, 1e-10)

        horizonte = 40 * b
        ts = list(range(1, horizonte))
        a = [[costo_noche(j, t, b) for j in range(1, horizonte)] + [-min(t, b)] for t in ts]
        costo = [0.0] * (horizonte - 1) + [1.0]
        res = linprog(costo, A_ub=a, b_ub=[0.0] * len(ts), A_eq=[[1.0] * (horizonte - 1) + [0.0]], b_eq=[1.0], bounds=[(0, None)] * horizonte)
        chequear(f"programa lineal B={b}", float(res.fun), c, 1e-6)
    print(f"    B=10: determinista {2 - 1 / 10:.4f}, con dado {1 / (1 - 0.9 ** 10):.4f}, límite {math.e / (math.e - 1):.4f}")
    print(f"    B=10: p_1={(1 / (1 - 0.9 ** 10)) / 10 * 0.9 ** 9:.4f}, p_10={(1 / (1 - 0.9 ** 10)) / 10:.4f}")


def mtf(seq, n):
    lista = list(range(n))
    total = 0
    for x in seq:
        i = lista.index(x)
        total += i + 1
        lista.insert(0, lista.pop(i))
    return total


def inversiones(a, b):
    pos = {v: i for i, v in enumerate(b)}
    return sum(1 for i in range(len(a)) for j in range(i + 1, len(a)) if pos[a[i]] > pos[a[j]])


def opt_lista(seq, n, estados, dist, adelante):
    d = np.full(len(estados), np.inf)
    d[estados.index(tuple(range(n)))] = 0
    for x in seq:
        d = (d[:, None] + dist).min(axis=0)
        d = d + np.array([e.index(x) + 1 for e in estados])
        nuevo = d.copy()
        for s, destinos in enumerate(adelante[x]):
            for u in destinos:
                nuevo[u] = min(nuevo[u], d[s])
        d = nuevo
    return float(d.min())


def lista():
    seccion("¿Dónde está el gato?: list update")
    rng = np.random.default_rng(SEED)
    n = 4
    estados = list(permutations(range(n)))
    indice = {e: i for i, e in enumerate(estados)}
    dist = np.array([[inversiones(a, b) for b in estados] for a in estados], dtype=float)
    adelante = {}
    for x in range(n):
        filas = []
        for e in estados:
            i = e.index(x)
            filas.append([indice[e[:k] + (x,) + e[k:i] + e[i + 1:]] for k in range(i)])
        adelante[x] = filas
    peor = -1.0
    pruebas = 0
    for _ in range(300):
        m = int(rng.integers(5, 40))
        seq = [int(v) for v in rng.integers(0, n, m)]
        opt = opt_lista(seq, n, estados, dist, adelante)
        c = mtf(seq, n)
        pruebas += 1
        if c > 2 * opt - m + 1e-9:
            chequear("cota de Sleator y Tarjan", c, 2 * opt - m)
        peor = max(peor, c / opt)
    chequear("violaciones de 2 OPT - m en pruebas aleatorias", 0.0, 0.0)
    print(f"    {pruebas} secuencias contra OPT exacto; razón máxima observada {peor:.4f}")

    for n2 in (2, 3, 6, 10):
        m = 6000 * n2
        lst = list(range(n2))
        total = 0
        conteo = [0] * n2
        for _ in range(m):
            x = lst[-1]
            conteo[x] += 1
            total += n2
            lst.insert(0, lst.pop())
        orden = sorted(conteo, reverse=True)
        estatico = sum((i + 1) * c for i, c in enumerate(orden))
        chequear(f"gato en último lugar n={n2}", total / estatico, 2 - 2 / (n2 + 1), 1e-9)
        chequear(f"promedio estático n={n2}", estatico / m, (n2 + 1) / 2, 1e-9)
    print(f"    n=6: {2 - 2 / 7:.4f}")


def cubre(rutas, gatos):
    conjunto = set(gatos)
    return sum(1 for r in rutas if conjunto & set(r))


def voraz(rutas, k, vertices):
    elegidos = []
    for _ in range(k):
        elegidos.append(max(vertices, key=lambda v: cubre(rutas, elegidos + [v])))
    return elegidos


def arbol_busqueda(rutas, k):
    if not rutas:
        return True
    if k == 0:
        return False
    return any(arbol_busqueda([r for r in rutas if v not in r], k - 1) for v in rutas[0])


def colonia():
    seccion("Colonia: interdicción")
    rng = np.random.default_rng(SEED)
    peor = 2.0
    for _ in range(1500):
        n = int(rng.integers(4, 10))
        m = int(rng.integers(3, 25))
        d = int(rng.integers(2, 5))
        vertices = list(range(n))
        rutas = [[int(v) for v in rng.choice(n, size=int(rng.integers(1, min(d, n) + 1)), replace=False)] for _ in range(m)]
        k = int(rng.integers(1, 5))
        opt = max(cubre(rutas, c) for c in combinations(vertices, k))
        g = cubre(rutas, voraz(rutas, k, vertices))
        peor = min(peor, g / opt)
        if g < (1 - (1 - 1 / k) ** k) * opt - 1e-9:
            chequear("cota voraz", g, (1 - (1 - 1 / k) ** k) * opt)
        factible = any(cubre(rutas, c) == m for c in combinations(vertices, k))
        if arbol_busqueda(rutas, k) != factible:
            chequear("árbol de búsqueda contra fuerza bruta", 1.0, 0.0)
        conteo = [sum(1 for r in rutas if v in r) for v in vertices]
        if max(conteo) != max(cubre(rutas, [v]) for v in vertices):
            chequear("un gato por conteo", 1.0, 0.0)
    chequear("violaciones de voraz y FPT", 0.0, 0.0)
    print(f"    peor razón voraz/óptimo observada: {peor:.4f}")

    for _ in range(300):
        n = int(rng.integers(3, 8))
        aristas = [(u, v) for u, v in combinations(range(n), 2) if rng.random() < 0.5]
        if not aristas:
            continue
        rutas = [list(e) for e in aristas]
        for k in range(n + 1):
            cobertura = any(all(u in c or v in c for u, v in aristas) for c in combinations(range(n), k))
            if cobertura != arbol_busqueda(rutas, k):
                chequear("reducción desde cobertura de vértices", 1.0, 0.0)
    chequear("violaciones de la reducción", 0.0, 0.0)
    for b in (2, 10, 50):
        chequear(f"producto de constantes B={b}", (1 - (1 - 1 / b) ** b) * (1 / (1 - (1 - 1 / b) ** b)), 1.0, 1e-12)


def menores():
    seccion("Resultados menores")
    h, mu_f, d = 0.75, 0.3, 0.05
    print(f"    ganancia del vaso: {h / (mu_f * d):.1f}")
    chequear("ganancia del vaso", h / (mu_f * d), 50.0)

    u, hh, rr, g = 0.78, 0.03, 0.005, 9.81
    fr = u / math.sqrt(g * rr)
    fr_star = (rr / hh) * fr ** (2 / 3)
    print(f"    Fr={fr:.3f}  Fr*={fr_star:.3f}")
    chequear("Fr* de Reis et al.", fr_star, 0.4, 0.02)
    f4 = 4.6 * 4 ** -0.181
    print(f"    f(4 kg)={f4:.3f} s^-1")
    chequear("frecuencia a 4 kg dentro de 3.5±0.4", f4, 3.5, 0.4)
    chequear("exponente -1/6 dentro de -0.181±0.024", -1 / 6, -0.181, 0.024)

    masa = 4.0
    vol = masa / 1000
    r_e = (3 * vol / (4 * math.pi)) ** (1 / 3)
    a_e = 4 * math.pi * r_e ** 2
    largo = 0.45
    r_c = math.sqrt(vol / (math.pi * largo))
    a_c = 2 * math.pi * r_c * largo + 2 * math.pi * r_c ** 2
    chequear("isoperimétrica", a_e, (36 * math.pi) ** (1 / 3) * vol ** (2 / 3), 1e-12)
    m_w = 70 * masa ** 0.75 * 4184 / 86400
    tb, tl = 38.5, 30.0
    conductancia = m_w / (a_e * (tb - tl))
    t_estirado = tb - m_w / (conductancia * a_c)
    deficit = conductancia * a_e * (tb - 20) - m_w
    print(f"    r_e={r_e * 100:.2f} cm  A_e={a_e:.4f} m2  r_c={r_c * 100:.2f} cm  A_c={a_c:.4f} m2  A_e/A_c={a_e / a_c:.3f}")
    print(f"    M={m_w:.2f} W  h={conductancia:.2f} W/m2K  T_estirado={t_estirado:.1f} C  déficit a 20 C={deficit:.1f} W")
    chequear("sin déficit en T_inferior", conductancia * a_e * (tb - tl) - m_w, 0.0, 1e-9)


def main():
    codigo()
    caida()
    persuasion()
    puerta()
    sol()
    noche()
    lista()
    colonia()
    menores()
    print()
    if fallos:
        print(f"FALLARON {len(fallos)}: {', '.join(fallos)}")
        return 1
    print("Todos los chequeos pasaron.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
