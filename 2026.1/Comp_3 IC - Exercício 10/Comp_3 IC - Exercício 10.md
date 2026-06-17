

## 1. Fundamentação Teórica

Seja $X_1, X_2, \dots, X_n \stackrel{\text{i.i.d.}}{\sim} N(\mu, \sigma^2)$, com $n = 3$.

O melhor estimador pontual para $\mu$ é a **média amostral**

$$
\bar{X} = \frac{1}{n}\sum_{i=1}^n X_i,
$$

que é não viesado ($\mathbb{E}[\bar{X}] = \mu$), tem variância $\sigma^2/n$ e atinge o limite inferior de Cramér-Rao, sendo portanto o **UMVUE** (Uniformly Minimum Variance Unbiased Estimator).

### Por que esses candidatos a pivô funcionam?

A distribuição normal pertence à **família locação-escala** (_location-scale family_). Isso significa que, se $X \sim N(\mu, \sigma^2)$, então para quaisquer constantes $a \in \mathbb{R}$ e $b > 0$,

$$
a + bX \sim N(a + b\mu,\; b^2\sigma^2).
$$

Em particular, a família é fechada sob **translação** (subtrair $\mu$) e **mudança de escala** (dividir por $\sigma$).

Nosso primeiro candidato a quantidade pivotal é

$$
Q = \bar{X} - \mu.
$$

Três observações são cruciais:

**1.** $\bar{X} - \mu$ **pertence à família escala (locação-escala)**. Como $\bar{X} \sim N(\mu,\, \sigma^2/n)$, a transformação $Q = \bar{X} - \mu$ é uma translação que desloca a média para $0$ sem alterar a forma da distribuição:

$$
Q = \bar{X} - \mu \sim N\!\left(0,\, \frac{\sigma^2}{n}\right).
$$

**2.** **A esperança de $Q$ é $0$**:

$$
\mathbb{E}[\bar{X} - \mu] = \mathbb{E}[\bar{X}] - \mu = \mu - \mu = 0.
$$

**3.** **A subtração de $\mu$ preserva a normalidade.** Este é o resultado-chave: se $X \sim N(\mu, \sigma^2)$, então $X - \mu \sim N(0, \sigma^2)$. A normal é fechada sob translação — subtrair $\mu$ apenas recentra a distribuição em $0$, mantendo a forma exatamente normal.

Isso mostra que $Q = \bar{X} - \mu$ é uma **quantidade pivotal**: sua distribuição não depende de $\mu$ (embora ainda dependa de $\sigma^2$).

Se **$\sigma$ for conhecido**, podemos ir além: dividir $Q$ pelo seu desvio padrão $\sigma/\sqrt{n}$ elimina também a dependência de $\sigma$, gerando o pivô canônico $Z \sim N(0,1)$:

$$
Z = \frac{\bar{X} - \mu}{\sigma/\sqrt{n}} \sim N(0,1).
$$

Se **$\sigma$ for desconhecido**, substituímos $\sigma$ por seu estimador $S$. A distribuição do pivô muda para $t$ de Student (como veremos a seguir), mas o princípio permanece o mesmo: a quantidade $\frac{\bar{X} - \mu}{S/\sqrt{n}}$ é pivô porque sua distribuição não depende de $\mu$ nem de $\sigma$.



### Pivô e distribuições

#### a) $\sigma$ conhecido

$$
Z = \frac{\bar{X} - \mu}{\sigma / \sqrt{n}} \sim N(0, 1).
$$

Intervalo de confiança a $100(1-\alpha)\%$:

$$
IC_{1-\alpha}(\mu) = \left[ \bar{X} - z_{\alpha/2} \frac{\sigma}{\sqrt{n}},\; \bar{X} + z_{\alpha/2} \frac{\sigma}{\sqrt{n}} \right],
$$

onde $z_{\alpha/2} = \Phi^{-1}(1 - \alpha/2)$.

Para $\alpha = 0.05$, $z_{0.025} \approx 1.96$, $\sigma = 5$, $n = 3$:

$$
IC_{95\%}(\mu) = \bar{X} \pm 1.96 \cdot \frac{5}{\sqrt{3}} \approx \bar{X} \pm 5.658.
$$

> A amplitude total do intervalo é $2 \times 5.658 = 11.316$, sempre positiva. A notação $\pm$ indica simetria em torno de $\bar{X}$.

#### b) $\sigma$ desconhecido

Estimamos $\sigma$ pelo **desvio padrão amostral**

$$
S = \sqrt{\frac{1}{n-1} \sum_{i=1}^n (X_i - \bar{X})^2}.
$$

O pivô resultante é

$$
T = \frac{\bar{X} - \mu}{S / \sqrt{n}} \sim t_{n-1},
$$

onde $t_{n-1}$ denota a distribuição **t de Student** com $n-1$ graus de liberdade.

Intervalo de confiança:

$$
IC_{1-\alpha}(\mu) = \left[ \bar{X} - t_{n-1;\, \alpha/2} \frac{S}{\sqrt{n}},\; \bar{X} + t_{n-1;\, \alpha/2} \frac{S}{\sqrt{n}} \right].
$$

Para $n = 3$, $t_{2;\,0.025} \approx 4.303$, significativamente maior que $z_{0.025} = 1.96$.

### Por que uma t de Student?

A substituição de $\sigma$ por $S$ altera a distribuição do pivô. A demonstração formal usa três fatos:

**1.** Padronização com $\sigma$ conhecido gera uma normal padrão:

$$
Z = \frac{\bar{X} - \mu}{\sigma / \sqrt{n}} \sim N(0,1).
$$

**2.** A soma de quadrados dos desvios, escalonada por $\sigma^2$, segue uma qui-quadrado com $n-1$ graus de liberdade (Teorema de Cochran):

$$
\frac{(n-1)S^2}{\sigma^2} = \frac{\sum_{i=1}^n (X_i - \bar{X})^2}{\sigma^2} \sim \chi^2_{n-1}.
$$

**3.** Independência entre $\bar{X}$ e $S^2$ (válida para amostras normais).

A estatística $T$ pode ser reescrita como

$$
T = \frac{\bar{X} - \mu}{S / \sqrt{n}}
    = \frac{ \displaystyle\frac{\bar{X} - \mu}{\sigma / \sqrt{n}} }
           { \displaystyle\sqrt{\frac{(n-1)S^2}{\sigma^2} \bigg/ (n-1)} }
    = \frac{Z}{\sqrt{\chi^2_{n-1} / (n-1)}} \sim t_{n-1}.
$$

Ou seja, $T$ é a razão entre uma **normal padrão** e a **raiz de uma qui-quadrado** escalada por seus graus de liberdade, onde numerador e denominador são independentes — exatamente a definição da distribuição t de Student.

**Consequência prática:** com $n=3$, a t possui **apenas 2 graus de liberdade**, o que resulta em caudas muito mais pesadas que a normal. O quantil $t_{2;0.025} \approx 4.303$ é mais que o dobro de $z_{0.025} = 1.96$, e a amplitude do IC dependerá também da variabilidade amostral de $S$.

### $S^2$ é UMVUE para $\sigma^2$?

Sim. A justificativa segue por **Lehmann-Scheffé**:

A família normal $N(\mu, \sigma^2)$ é da **família exponencial** com estatística suficiente $(\sum X_i, \sum X_i^2)$ — equivalente a $(\bar{X}, S^2)$. Essa estatística é **completa**.

Como $S^2$ é não viesado ($\mathbb{E}[S^2] = \sigma^2$) e função dessa suficiente completa, pelo Teorema de Lehmann-Scheffé, $S^2$ é **UMVUE** para $\sigma^2$.


### Independência $\bar{X} \perp S^2$ (via Basu)

O argumento usa o **Teorema de Basu**:

**1.** $\bar{X}$ é **suficiente completa** para $\mu$ (na família normal com $\sigma^2$ fixo, a forma de família exponencial garante isso).

**2.** $S^2$ é **ancilar** para $\mu$ — sua distribuição independe de $\mu$:

$$
\frac{(n-1)S^2}{\sigma^2} \sim \chi^2_{n-1}.
$$

Note que $\mu$ não aparece. Isso faz sentido: os desvios $(X_i - \bar{X})$ eliminam $\mu$, então qualquer função deles (como $S^2$) tem distribuição livre de $\mu$.

**3.** Pelo **Teorema de Basu**: se $T$ é suficiente completa e $U$ é ancilar, então $T \perp U$.

Portanto $\bar{X} \perp S^2$.

E esse é o terceiro ingrediente que fecha a demonstração do pivô $t$ — junto com a normalidade de $\bar{X}$ e a qui-quadrado de $(n-1)S^2/\sigma^2$, a independência garante que

$$
T = \frac{Z}{\sqrt{V/(n-1)}} \sim t_{n-1},
$$
onde $Z \sim N(0,1)$ e $V \sim \chi^2_{n-1}$.

## 2. Simulação

```#| label: bibliotecas
#| message: false
#| warning: false
library(ggplot2)
library(patchwork)
```

```#| label: simulacao
#| message: false
#| warning: false
# ============================================================
# Parâmetros e simulação
# ============================================================

set.seed(42)

# Parâmetros populacionais
mu    <- 70
sigma <- 5
n     <- 3
N     <- 100
alpha <- 0.05

# Quantis
z_alpha2 <- qnorm(1 - alpha/2)
t_alpha2 <- qt(1 - alpha/2, df = n - 1)

# Data frame para armazenar os resultados
resultados <- data.frame(
  amostra      = integer(N),
  xbar         = numeric(N),
  s            = numeric(N),
  ic_a_inf     = numeric(N),
  ic_a_sup     = numeric(N),
  ic_b_inf     = numeric(N),
  ic_b_sup     = numeric(N),
  captura_a    = logical(N),
  captura_b    = logical(N)
)

for (i in 1:N) {
  x <- rnorm(n, mean = mu, sd = sigma)
  xbar <- mean(x)
  s    <- sd(x)
  erro_a <- z_alpha2 * sigma / sqrt(n)
  erro_b <- t_alpha2 * s / sqrt(n)
  captura_a <- (xbar - erro_a <= mu) & (mu <= xbar + erro_a)
  captura_b <- (xbar - erro_b <= mu) & (mu <= xbar + erro_b)

  resultados[i, ] <- list(i, xbar, s,
                          xbar - erro_a, xbar + erro_a,
                          xbar - erro_b, xbar + erro_b,
                          captura_a, captura_b)
}

# Proporção observada de cobertura
prop_a <- mean(resultados$captura_a)
prop_b <- mean(resultados$captura_b)
```

```#| label: cobertura
cat("=== Cobertura observada ===\n")
cat(sprintf("Sigma conhecido:    %.1f%%\n", prop_a * 100))
cat(sprintf("Sigma desconhecido: %.1f%%\n", prop_b * 100))
cat(sprintf("Nominal (teórica):  %.1f%%\n", (1 - alpha) * 100))
```

## 3. Gráfico dos 100 intervalos

```#| label: grafico-ic
#| fig-cap: "Intervalos de confiança de 95% para mu (100 amostras, n=3)"
#| message: false
#| warning: false

resultados$ordem <- 1:N

# Painel (a): sigma conhecido
p_a <- ggplot(resultados) +
  geom_segment(aes(x = ic_a_inf, xend = ic_a_sup,
                   y = ordem, yend = ordem,
                   color = captura_a), linewidth = 1) +
  geom_vline(xintercept = mu, linetype = "dashed",
             color = "darkred", linewidth = 0.8) +
  scale_color_manual(
    values = c("TRUE" = "steelblue", "FALSE" = "red"),
    labels = c("TRUE" = "Contém mu", "FALSE" = "Não contém mu")
  ) +
  labs(title = expression("IC"[95*"%"]~"para"~mu~"(sigma conhecido)"),
       subtitle = sprintf("Cobertura: %.1f%% (nominal: 95%%)", prop_a * 100),
       x = expression(mu), y = "Amostra", color = "") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

# Painel (b): sigma desconhecido
p_b <- ggplot(resultados) +
  geom_segment(aes(x = ic_b_inf, xend = ic_b_sup,
                   y = ordem, yend = ordem,
                   color = captura_b), linewidth = 1) +
  geom_vline(xintercept = mu, linetype = "dashed",
             color = "darkred", linewidth = 0.8) +
  scale_color_manual(
    values = c("TRUE" = "steelblue", "FALSE" = "red"),
    labels = c("TRUE" = "Contém mu", "FALSE" = "Não contém mu")
  ) +
  labs(title = expression("IC"[95*"%"]~"para"~mu~"(sigma desconhecido - t de Student)"),
       subtitle = sprintf("Cobertura: %.1f%% (nominal: 95%%)", prop_b * 100),
       x = expression(mu), y = "Amostra", color = "") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

# Painel lado a lado
p_a + p_b +
  plot_annotation(title = "Comparação: IC conhecido vs. desconhecido (n = 3)",
                  theme = theme(plot.title = element_text(hjust = 0.5)))
```

## 4. Gráfico das distribuições

```#| label: grafico-dist
#| fig-cap: "Distribuições N(0,1) e t2 com quantis críticos"
x_vals <- seq(-6, 6, length.out = 301)

df_dist <- data.frame(
  x = rep(x_vals, 2),
  densidade = c(dnorm(x_vals), dt(x_vals, df = 2)),
  distribuicao = factor(
    rep(c("N(0,1)", "t[2]"), each = length(x_vals))
  )
)

ggplot(df_dist, aes(x = x, y = densidade, color = distribuicao)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = c(-z_alpha2, z_alpha2),
             color = "steelblue", linetype = "dotted", linewidth = 0.6) +
  geom_vline(xintercept = c(-t_alpha2, t_alpha2),
             color = "tomato", linetype = "dotted", linewidth = 0.6) +
  annotate("text", x = z_alpha2 + 0.4, y = 0.05,
           label = sprintf("z[0.025] == %.2f", z_alpha2),
           parse = TRUE, color = "steelblue", size = 3.5, hjust = 0) +
  annotate("text", x = t_alpha2 + 0.4, y = 0.10,
           label = sprintf("t[2, 0.025] == %.2f", t_alpha2),
           parse = TRUE, color = "tomato", size = 3.5, hjust = 0) +
  scale_color_manual(
    values = c("N(0,1)" = "steelblue", "t[2]" = "tomato"),
    labels = expression(N(0, 1), t[2])
  ) +
  labs(title = "Distribuições envolvidas nos pivôs",
       x = "Estatística", y = "Densidade", color = "") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")
```

## 5. Discussão dos Resultados

| Aspecto | (a) $\sigma$ conhecido | (b) $\sigma$ desconhecido |
|---------|------------------------|---------------------------|
| Distribuição do pivô | $N(0,1)$ | $t_{2}$ |
| Quantil crítico | $z_{0.025} \approx 1.96$ | $t_{2;0.025} \approx 4.30$ |
| Amplitude do IC | Constante: $\pm 5.66$ | Variável (depende de $S$) |
| Cobertura esperada | ~95% | ~95% |
| Observação | $\sigma$ raramente conhecido na prática | Mais realista; IC mais largo |

Com $n=3$, a $t_2$ tem caudas muito pesadas. O quantil 4.30 significa que, para capturar 95% de confiança, o intervalo precisa se estender a mais de 4 erros padrão da média. Como consequência, os ICs em (b) serão **substancialmente mais largos** que em (a), embora ambos mantenham (em média) a cobertura nominal de 95%.

Note que a principal diferença prática entre as duas situações não está na cobertura (ambos os métodos produzem intervalos que contêm $\mu$ em aproximadamente 95% das repetições), mas sim na amplitude: enquanto em (a) a semi-amplitude é fixa em 5.66, em (b) ela varia de amostra para amostra porque $S$ é estimado, refletindo a incerteza adicional de não conhecermos $\sigma$.
