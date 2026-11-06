# Sobre la complejidad gatuna

> *Se suele decir que el gato doméstico es perezoso. Muestro que esta afirmación es un error de modelado.*

Un paper de 15 páginas, escrito como si fuera en serio, que demuestra que el gato no evita el trabajo: lo externaliza. Es una parodia, pero los teoremas tienen demostración. Lo caricaturesco son los modelos.

- [`paper.pdf`](paper.pdf), en español
- [`paper_en.pdf`](paper_en.pdf), *On Feline Complexity*, en inglés

## Qué se demuestra

| Conducta | Herramienta | Resultado |
|---|---|---|
| Te mira fijo | Complejidad de Kolmogorov | Casi ningún antojo es comprimible |
| Cae de pie | Holonomía, gauge abeliano | La rotación neta es el área encerrada en el espacio de formas |
| Maúlla ya saciado | Persuasión bayesiana | Come con probabilidad $\mu_0/\mu^\star$ sin mentir |
| Pide la puerta y no sale | Opción real, frontera libre | Se queda en el umbral para siempre con probabilidad $1/\varphi$ |
| Persigue el sol | Teorema del valor marginal | $T^\star=(3cW^2/2Rv^2)^{1/3}$ |
| Maúlla a las 3 a.m. | *Ski rental* | Razón $2-1/B$, y el azar no ayuda si te está mirando |
| Siempre está en el último lugar donde buscas | *List update* | Cota inferior $2-2/(n+1)$ |
| Una colonia te bloquea el paso | Cobertura de vértices | NP-completo; voraz $1-1/e$, sin mejora posible salvo que P = NP |
| Se sienta en tu teclado | Desigualdad isoperimétrica | Le faltan unos 11 W, más o menos lo que disipa una laptop |

## Revisión por patas

`revision_por_patas.py` hace de revisión por pares: comprueba cada número del paper con cálculo simbólico, simulación de Monte Carlo, programación lineal y fuerza bruta contra el óptimo exacto.

```sh
pip install -r requirements.txt
python3 revision_por_patas.py
```

Para recompilar los PDF hace falta una distribución de TeX con `pgfplots` y `babel-spanish`:

```sh
make          # revisa y compila ambos papers
make check    # solo la revisión
make clean    # borra los auxiliares de LaTeX
```

## Honestidad

Los parámetros marcados con ‡ están calibrados a ojo. Las afirmaciones empíricas citadas vienen de artículos reales, como el de *Science* sobre cómo beben los gatos. Ningún gato fue dejado caer.

## Créditos

A Miko, mi gato, que inspiró el trabajo y no ha leído el manuscrito. La idea surgió de un reel de [blackboardguys](https://www.instagram.com/blackboardguys/) que aplicaba teoría de juegos a las *situationships*.

---

## English

*On Feline Complexity* is a 15-page parody paper whose theorems are actually proved: cats are not lazy, they outsource. The results cover holonomy in mid-air righting, Bayesian persuasion in meowing, real options at the door, ski rental at 3 a.m., list update when searching for the cat, and NP-completeness for colonies. Read [`paper_en.pdf`](paper_en.pdf). `revision_por_patas.py` (Spanish for *review by paws*) checks every number in the paper.

## Licencia

[GPL-2.0](LICENSE).
