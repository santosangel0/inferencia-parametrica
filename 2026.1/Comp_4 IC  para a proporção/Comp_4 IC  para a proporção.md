
```#| label: libs
#| message: false
#| warning: false
library(gt)
library(tibble)
```

## 1. Fundamentação Teórica

Seja $X_1, X_2, \dots, X_n \stackrel{\text{i.i.d.}}{\sim} \text{Bernoulli}(p)$, com $p \in (0,1)$. Queremos construir intervalos de confiança para a probabilidade de sucesso $p$.

O estimador natural é a **proporção amostral**

$$
\hat{p} = \frac{1}{n}\sum_{i=1}^n X_i,
$$

que é não viesado ($\mathbb{E}[\hat{p}] = p$), tem variância $\operatorname{Var}(\hat{p}) = p(1-p)/n$ e é o **EMV** (estimador de máxima verossimilhança) para $p$.

### Distribuição assintótica de $\hat{p}$

Pelo **Teorema Central do Limite** (TCL), para $n$ suficientemente grande:

$$
\frac{\hat{p} - p}{\sqrt{\frac{p(1-p)}{n}}} \xrightarrow{d} N(0,1).
$$

Esse resultado é a base de todos os ICs assintóticos para proporção. O problema é que a variância $p(1-p)/n$ depende do próprio parâmetro $p$, que é desconhecido. Cada método resolve essa dependência de forma diferente.

### O IC de Wald

O método de **Wald** substitui $p$ por $\hat{p}$ no erro padrão:

$$
\text{ep}_{\text{wald}} = \sqrt{\frac{\hat{p}(1-\hat{p})}{n}}.
$$

O intervalo é dado por

$$
IC_{\text{Wald}} = \left[ \hat{p} - z_{\alpha/2} \cdot \text{ep}_{\text{wald}},\; \hat{p} + z_{\alpha/2} \cdot \text{ep}_{\text{wald}} \right],
$$

onde $z_{\alpha/2} = \Phi^{-1}(1 - \alpha/2)$.

**Problema crítico:** Quando $p$ é próximo de $0$ ou $1$ e $n$ é pequeno, a probabilidade de $\hat{p}=0$ (ou $\hat{p}=1$) é alta. Nesse caso:

- $\hat{p}=0 \Rightarrow \text{ep}_{\text{wald}} = 0 \Rightarrow \text{intervalo} = [0, 0]$
- Esse intervalo degenerado nunca contém $p>0$

A taxa de cobertura real pode cair para valores **muito abaixo** do nominal. Esse fenômeno é conhecido como "colapso de Wald".

### Estabilização da variância — Arcoseno

A transformação **arcoseno** (ou arco-seno) usa o **método delta** para estabilizar a variância.

Pelo método delta, se $\sqrt{n}(\hat{p} - p) \xrightarrow{d} N(0, p(1-p))$, então para uma função $g$ diferenciável:

$$
\sqrt{n}\bigl(g(\hat{p}) - g(p)\bigr) \xrightarrow{d} N\!\left(0,\; [g'(p)]^2 p(1-p)\right).
$$

Escolhendo $g(p) = \arcsin(\sqrt{p})$, temos:

$$
g'(p) = \frac{1}{2\sqrt{p(1-p)}} \quad\Rightarrow\quad [g'(p)]^2 p(1-p) = \frac{1}{4}.
$$

Portanto:

$$
\sqrt{n}\bigl(\arcsin(\sqrt{\hat{p}}) - \arcsin(\sqrt{p})\bigr) \xrightarrow{d} N\!\left(0,\; \frac{1}{4}\right),
$$

ou equivalentemente:

$$
\operatorname{Var}\bigl[\arcsin(\sqrt{\hat{p}})\bigr] \approx \frac{1}{4n}.
$$

A variância **não depende mais de $p$** — daí o nome "estabilização da variância".

**Intervalo:** Construímos o IC na escala transformada e depois aplicamos a transformação inversa:

1. Na escala transformada:
   $$
   \arcsin(\sqrt{\hat{p}}) \pm \frac{z_{\alpha/2}}{2\sqrt{n}}.
   $$

2. Transformação reversa ($\sin^2$):
   $$
   IC_{\text{arc}} = \left[ \sin^2\!\bigl(\arcsin(\sqrt{\hat{p}}) - \tfrac{z_{\alpha/2}}{2\sqrt{n}}\bigr),\;
                     \sin^2\!\bigl(\arcsin(\sqrt{\hat{p}}) + \tfrac{z_{\alpha/2}}{2\sqrt{n}}\bigr) \right].
   $$

**Vantagem:** Mesmo quando $\hat{p}=0$, temos $\arcsin(0) = 0$ e o limite inferior se torna:

$$
\sin^2\!\left(-\frac{z_{\alpha/2}}{2\sqrt{n}}\right) > 0,
$$

evitando o intervalo degenerado $[0,0]$.

### Método Conservador (Variância Máxima)

O método **conservador** substitui $p(1-p)$ por seu valor máximo possível:

$$
p(1-p) \leq \frac{1}{4} \quad\text{(atingido em $p=0.5$)}.
$$

Assim:

$$
\text{ep}_{\text{cons}} = \frac{1}{2\sqrt{n}},
$$

e o intervalo é

$$
IC_{\text{cons}} = \left[ \hat{p} - \frac{z_{\alpha/2}}{2\sqrt{n}},\; \hat{p} + \frac{z_{\alpha/2}}{2\sqrt{n}} \right].
$$

**Interpretação:** Ao usar a maior variância possível, o método garante **supercobertura** — a taxa real de cobertura nunca cai abaixo do nominal. O custo é uma **largura excessiva**, especialmente quando $p$ está longe de $0.5$.

## 2. Simulação

```#| label: simulacao
#| message: false
#| warning: false

# ============================================================
# Parâmetros da simulação
# ============================================================

set.seed(42)

p_valores <- c(0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.95, 0.99)
n_valores <- c(10, 30, 50, 100, 500)
gama      <- 0.95
M         <- 1e4
z         <- qnorm((1 + gama) / 2)

cenarios <- expand.grid(p = p_valores, n = n_valores)

# ============================================================
# Função de simulação
# ============================================================

simular_todos_metodos <- function(p, n, M, z, gama) {
  p_chapeu <- rbinom(n = M, size = n, prob = p) / n

  # (a) Wald
  ep_wald  <- sqrt((p_chapeu * (1 - p_chapeu)) / n)
  inf_wald <- p_chapeu - z * ep_wald
  sup_wald <- p_chapeu + z * ep_wald
  cob_wald <- mean(inf_wald <= p & p <= sup_wald)
  amp_wald <- mean(sup_wald - inf_wald)
  vies_wald <- cob_wald - gama

  # (b) Arcoseno
  termo_escala <- z / (2 * sqrt(n))
  inf_arc <- sin(asin(sqrt(p_chapeu)) - termo_escala)^2
  sup_arc <- sin(asin(sqrt(p_chapeu)) + termo_escala)^2
  cob_arc <- mean(inf_arc <= p & p <= sup_arc)
  amp_arc <- mean(sup_arc - inf_arc)
  vies_arc <- cob_arc - gama

  # (c) Conservador
  ep_cons  <- 1 / (2 * sqrt(n))
  inf_cons <- p_chapeu - z * ep_cons
  sup_cons <- p_chapeu + z * ep_cons
  cob_cons <- mean(inf_cons <= p & p <= sup_cons)
  amp_cons <- mean(sup_cons - inf_cons)
  vies_cons <- cob_cons - gama

  return(c(
    wald_cobertura = cob_wald, wald_largura = amp_wald, wald_vies = vies_wald,
    arc_cobertura  = cob_arc,  arc_largura  = amp_arc,  arc_vies  = vies_arc,
    cons_cobertura = cob_cons, cons_largura = amp_cons, cons_vies = vies_cons
  ))
}

# ============================================================
# Execução
# ============================================================

resultados_brutos <- sapply(1:nrow(cenarios), function(i) {
  simular_todos_metodos(
    p = cenarios$p[i], n = cenarios$n[i],
    M = M, z = z, gama = gama
  )
})

df_final <- cbind(cenarios, t(resultados_brutos))
```

## 3. Resultados — Tabela

```#| label: tabela-completa
#| message: false
#| warning: false

df_final %>%
  gt() %>%
  tab_header(
    title = md("**Resultados da Simulação: IC para Proporção**"),
    subtitle = md("Cobertura, Largura e Viés para métodos Wald, Arcoseno e Conservador")
  ) %>%
  fmt_number(
    columns = -c(p, n),
    decimals = 4
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(columns = c(p, n))
  ) %>%
  opt_align_table_header(align = "left")
```

## 4. Discussão dos Resultados

### O Colapso de Wald ($p = 0.01$, $n = 10$)

Na primeira linha do `df_final`, a coluna `wald_cobertura` exibe um valor **assustadoramente baixo**: **9,68%** (para um nível nominal de 95%).

**Causa:** Com $p$ muito pequeno e $n$ pequeno, a probabilidade de obter $\hat{p}=0$ (nenhum sucesso observado) é extremamente alta:

$$
P(\hat{p}=0) = (1-p)^n = 0.99^{10} \approx 0.9044 \quad (90,44\%).
$$

Quando $\hat{p}=0$, o erro padrão estimado de Wald é $\text{ep}_{\text{wald}} = 0$, gerando o intervalo degenerado $[0,0]$. Como $p = 0.01 > 0$, esse intervalo não contém o verdadeiro parâmetro. Isso limita a taxa máxima de cobertura possível para Wald a $1 - 0.9044 = 9.56\%$ (muito próxima do resultado empírico de **9,68%**, com a variação devida ao erro amostral de simulação).

Consequências:
- `wald_cobertura` = **9,68%** (viesado negativamente com `wald_vies` = **-85,32%**)

> O IC de Wald é **inutilizável** para proporções próximas a 0 ou 1 com amostras pequenas.

### O Comportamento do Arcoseno

A transformação arcoseno usa a estabilização de variância ($\operatorname{Var}[\arcsin(\sqrt{\hat{p}})] \approx 1/(4n)$) para contornar a dependência direta de $p$ no erro padrão. Os resultados empíricos nos trazem dois fatos importantes:

1. **Fronteira extrema e amostra minúscula ($p=0.01, n=10$):**
   A cobertura empírica do Arcoseno foi de **9,29%** (muito similar à de Wald).
   
   *Por que isso ocorre?* Quando $\hat{p}=0$, a transformação reversa $\sin^2\left(\arcsin(0) \pm \frac{z}{2\sqrt{n}}\right)$ gera os limites inferior e superior idênticos a $\sin^2\left(\frac{z}{2\sqrt{n}}\right) \approx \sin^2(0.3099) \approx 0.0930$. O intervalo resultante é $[0.0930, 0.0930]$. Como $p = 0.01 < 0.0930$, esse intervalo de largura zero também erra o alvo. Portanto, sob $n=10$ e $p=0.01$, o método do arcoseno também é penalizado pela alta probabilidade de $\hat{p}=0$ ($90,44\%$), não conseguindo superar o colapso nesse caso extremo.

2. **Superioridade em amostras moderadas e médias ($n \ge 50$):**
   À medida que o tamanho da amostra cresce, a aproximação normal na escala transformada do Arcoseno se mostra muito superior:
   - Para **$n = 100$ e $p = 0.05$**: Wald obteve **87,76%** de cobertura, enquanto o Arcoseno alcançou **95,17%** (aderindo ao nível nominal).
   - Para **$n = 500$ e $p = 0.01$**: Wald obteve **87,05%** de cobertura, enquanto o Arcoseno obteve **94,87%**.

### O Custo do Conservador

O método **conservador** garante que a cobertura real fique sempre igual ou acima da nominal (gerando supercobertura) ao preço de assumir a variância máxima ($p=0.5$).

- **Cenário central ($p = 0.5$):** Conservador e Wald se comportam de forma idêntica porque $\hat{p}(1-\hat{p})$ atinge o máximo $0.25$ justamente em $p=0.5$. Para $n=500$ e $p=0.5$, ambos têm cobertura de **94,62%** e largura de **0.0877**.
- **Cenário de cauda ($p = 0.01$):**
  - A cobertura trava em **100,00%** (intervalo superprotegido).
  - No entanto, a largura média revela o desperdício de precisão:
    - Com **$n = 10$**: a largura do Conservador é **0.6198** (fixa), o que é **17 vezes maior** que as larguras de Wald (**0.0365**) e do Arcoseno (**0.0342**).
    - Com **$n = 500$**: a largura do Conservador é **0.0877**, que é **5,2 vezes maior** que a largura dos outros métodos (**0.0169**).

### Considerações Práticas e de Consultoria

Ao analisar o uso prático de cada método, destaca-se que:

- **O método de Wald**, apesar de apresentar o pior desempenho empírico nas caudas e em amostras pequenas (subcobertura acentuada), continua sendo o **mais usual e difundido**. Isso se deve exclusivamente à sua **simplicidade teórica e de cálculo**, sendo a fórmula clássica ensinada nos cursos introdutórios de estatística e implementada por padrão na maioria das ferramentas e pacotes computacionais.
- **O método do Arcoseno** oferece um excelente ganho de precisão e correção da assimetria para amostras moderadas nas caudas. Contudo, possui uma **barreira de comunicação na prática**: a transformação de estabilização de variância baseada no método delta ($\arcsin(\sqrt{\hat{p}})$) é de difícil compreensão e justificativa em uma consultoria para tomadores de decisão ou clientes sem formação quantitativa avançada.
- **O método Conservador**, por outro lado, é de simples explicação estatística (baseia-se na proteção contra o pior cenário de variabilidade), mas seu uso prático é severamente limitado pela grande perda de precisão (intervalos largos demais), o que frequentemente o torna pouco informativo para tomadas de decisão que exijam estimativas mais refinadas.
