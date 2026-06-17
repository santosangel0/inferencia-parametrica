
## 1. Fundamentação Teórica: Exercício 12

Neste estudo, abordamos a construção de intervalos de confiança para a variância $\sigma^2$ de uma população com distribuição normal. O problema proposto baseia-se no Exercício 12 da lista de Intervalos de Confiança:

> **Enunciado:** Se uma amostra aleatória de 25 escores de QI de certa população acusa $s^2 = 120$, construa um intervalo de confiança de $95\%$ para $\sigma^2$ da população, admitindo que a variável aleatória de interesse segue uma distribuição normal.

Abaixo, descrevemos as duas metodologias para a construção deste intervalo de confiança: a usual (bicaudal simétrica na probabilidade das caudas) e a de comprimento mínimo (ótimos em termos de tamanho do intervalo).

### Quantidade Pivotal
Seja $X_1, \dots, X_n \sim \mathcal{N}(\mu, \sigma^2)$ uma amostra aleatória. A quantidade pivotal para a variância populacional $\sigma^2$ é dada por:
$$Y = \frac{(n-1)S^2}{\sigma^2} \sim \chi^2_{n-1}$$
onde $S^2 = \frac{1}{n-1}\sum_{i=1}^n (X_i - \bar{X})^2$ é a variância amostral e $\chi^2_{n-1}$ representa a distribuição chi-quadrado com $\nu = n-1$ graus de liberdade.

Desejamos encontrar limites críticos $a$ e $b$ tais que:
$$P\left( a \le \frac{(n-1)S^2}{\sigma^2} \le b \right) = 1-\alpha$$
o que nos permite escrever o intervalo de confiança geral como:
$$IC(\sigma^2) = \left[ \frac{(n-1)s^2}{b}, \, \frac{(n-1)s^2}{a} \right]$$

### Abordagem Usual (Bicaudal Simétrica / Caudas Iguais)
Na abordagem clássica, a probabilidade de erro $\alpha$ é dividida igualmente entre as duas caudas da distribuição, ou seja, alocamos $\alpha/2$ na cauda esquerda e $\alpha/2$ na cauda direita. Logo, os limites críticos são os quantis da distribuição chi-quadrado:
$$a_{\text{usual}} = \chi^2_{n-1, \alpha/2} \quad \text{e} \quad b_{\text{usual}} = \chi^2_{n-1, 1-\alpha/2}$$

Para o Exercício 12, temos os seguintes parâmetros:
* Tamanho amostral: $n = 25 \implies \nu = n-1 = 24$ graus de liberdade.
* Variância amostral observada: $s^2 = 120$.
* Nível de confiança: $1-\alpha = 0.95 \implies \alpha = 0.05 \implies \alpha/2 = 0.025$.

Consultando os valores críticos da distribuição $\chi^2_{24}$:
* $a_{\text{usual}} = \chi^2_{24, 0.025} \approx 12.40115$
* $b_{\text{usual}} = \chi^2_{24, 0.975} \approx 39.36408$

Substituindo os valores, o intervalo usual é:
$$IC_{\text{usual}}(\sigma^2) = \left[ \frac{24 \cdot 120}{39.36408}, \, \frac{24 \cdot 120}{12.40115} \right] \approx [73.16315, \, 232.23652]$$

O comprimento deste intervalo usual é:
$$L_{\text{usual}} = 232.23652 - 73.16315 = 159.07337$$

### Abordagem de Comprimento Mínimo (Otimização)
O comprimento de qualquer intervalo de confiança para a variância $\sigma^2$ é dada por:
$$L = (n-1)s^2 \left( \frac{1}{a} - \frac{1}{b} \right)$$

Para encontrar os valores de $a$ e $b$ que minimizam o comprimento $L$ mantendo o mesmo nível de confiança $1-\alpha$, formulamos o seguinte problema de otimização restrita:
$$\text{Minimizar } g(a, b) = \frac{1}{a} - \frac{1}{b}$$
$$\text{Sujeito a } F_{\nu}(b) - F_{\nu}(a) = 1-\alpha$$
onde $F_{\nu}$ é a função de distribuição acumulada (CDF) da $\chi^2_{\nu}$.

Utilizando o método de multiplicadores de Lagrange ou derivando a função com relação a $a$ (sabendo que $b$ é implicitamente uma função de $a$), a condição necessária de primeira ordem para o mínimo exige que:
$$a^2 f_{\nu}(a) = b^2 f_{\nu}(b)$$
onde $f_{\nu}(y)$ é a função de densidade de probabilidade (PDF) da $\chi^2_{\nu}$.

Sabendo que a densidade da distribuição chi-quadrado é:
$$f_{\nu}(y) = \frac{1}{2^{\nu/2}\Gamma(\nu/2)} y^{\nu/2-1} e^{-y/2}$$

A condição $a^2 f_{\nu}(a) = b^2 f_{\nu}(b)$ simplifica-se para:
$$a^{\nu/2+1} e^{-a/2} = b^{\nu/2+1} e^{-b/2}$$

Esta última igualdade é equivalente a exigir que as ordenadas de uma densidade chi-quadrado com $\nu+4$ graus de liberdade sejam iguais nos pontos $a$ e $b$:
$$f_{\nu+4}(a) = f_{\nu+4}(b)$$

Desta forma, os valores críticos ótimos $a_{\text{opt}}$ e $b_{\text{opt}}$ são obtidos resolvendo numericamente o seguinte sistema de equações não-lineares:
$$\begin{cases}
F_{\nu}(b) - F_{\nu}(a) = 1-\alpha \\
f_{\nu+4}(a) = f_{\nu+4}(b)
\end{cases}$$

Para o Exercício 12 ($\nu = 24$, $1-\alpha = 0.95$):
* Resolvendo o sistema, encontramos:
  * $a_{\text{opt}} \approx 13.52272$
  * $b_{\text{opt}} \approx 44.48023$

Substituindo na fórmula, o intervalo de confiança ótimo é:
$$IC_{\text{opt}}(\sigma^2) = \left[ \frac{24 \cdot 120}{44.48023}, \, \frac{24 \cdot 120}{13.52272} \right] \approx [64.74786, \, 212.97486]$$

O comprimento deste intervalo ótimo é:
$$L_{\text{opt}} = 212.97486 - 64.74786 = 148.22700$$

A redução absoluta no comprimento do intervalo é de **$10.84637$ unidades**, o que representa um ganho de eficiência de **$6.82\%$** mantendo exatamente a mesma cobertura nominal de $95\%$.

---

```#| message: false
#| warning: false
library(gt)
library(tibble)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
```

## 2. Implementação R para Otimização

Abaixo, fornecemos o código R necessário para determinar os valores críticos e construir o Intervalo de Confiança de comprimento mínimo. A função utiliza a rotina `uniroot()` do R para encontrar as raízes numéricas do sistema de otimização.

```#| code-fold: false
encontrar_valores_criticos_otimos <- function(nu, conf_level = 0.95) {
  alpha <- 1 - conf_level
  
  # Função objetivo para encontrar o valor crítico inferior 'a'
  # b é definido implicitamente pela restrição F(b) - F(a) = 1 - alpha
  funcao_objetivo <- function(a) {
    b <- qchisq(pchisq(a, df = nu) + conf_level, df = nu)
    dchisq(b, df = nu + 4) - dchisq(a, df = nu + 4)
  }
  
  # Intervalo de busca para a: entre 0 e o quantil alpha
  limite_superior_a <- qchisq(alpha, df = nu)
  
  sol <- uniroot(funcao_objetivo, interval = c(1e-5, limite_superior_a - 1e-5))
  a_opt <- sol$root
  b_opt <- qchisq(pchisq(a_opt, df = nu) + conf_level, df = nu)
  
  return(c(a = a_opt, b = b_opt))
}

# Exemplo de uso para o Exercício 12
criticos <- encontrar_valores_criticos_otimos(nu = 24, conf_level = 0.95)
cat("Valores Críticos Ótimos (Comprimento Mínimo):\n")
cat("a =", criticos["a"], "\n")
cat("b =", criticos["b"], "\n")
```

---

## 3. Simulação Monte Carlo e Resultados

Para validar empiricamente o método desenvolvido, implementamos uma simulação de Monte Carlo. Simulamos $B = 10000$ variâncias amostrais $S^2$ a partir do modelo $\frac{\nu S^2}{\sigma^2} \sim \chi^2_{\nu}$, sob a hipótese de que a verdadeira variância populacional seja $\sigma^2 = 120$ e o tamanho amostral seja $n = 25$.

```simular_ic_variancia <- function(sigma2_real = 120, n = 25, B = 10000) {
  nu <- n - 1
  conf_level <- 0.95
  alpha <- 1 - conf_level
  
  # 1. Valores críticos usuais
  a_usual <- qchisq(alpha / 2, df = nu)
  b_usual <- qchisq(1 - alpha / 2, df = nu)
  
  # 2. Valores críticos ótimos
  criticos <- encontrar_valores_criticos_otimos(nu, conf_level)
  a_opt <- criticos["a"]
  b_opt <- criticos["b"]
  
  # 3. Geração das variâncias amostrais simuladas (distribuição pivotal)
  # s^2 = (sigma^2 * chi^2_nu) / nu
  s2_sim <- (sigma2_real * rchisq(B, df = nu)) / nu
  
  # Intervalos Usuais
  ic_usual_lower <- (nu * s2_sim) / b_usual
  ic_usual_upper <- (nu * s2_sim) / a_usual
  cov_usual <- ic_usual_lower <= sigma2_real & sigma2_real <= ic_usual_upper
  len_usual <- ic_usual_upper - ic_usual_lower
  
  # Intervalos Ótimos
  ic_opt_lower <- (nu * s2_sim) / b_opt
  ic_opt_upper <- (nu * s2_sim) / a_opt
  cov_opt <- ic_opt_lower <= sigma2_real & sigma2_real <= ic_opt_upper
  len_opt <- ic_opt_upper - ic_opt_lower
  
  # Sumário dos Resultados
  res_summary <- data.frame(
    Metodo = c("Usual (Caudas Iguais)", "Comprimento Mínimo"),
    Valor_a = c(a_usual, a_opt),
    Valor_b = c(b_usual, b_opt),
    Cobertura = c(mean(cov_usual), mean(cov_opt)),
    Comprimento_Medio = c(mean(len_usual), mean(len_opt)),
    Comprimento_Medio_Teorico = c(
      sigma2_real * nu * (1/a_usual - 1/b_usual),
      sigma2_real * nu * (1/a_opt - 1/b_opt)
    )
  )
  
  return(list(
    summary = res_summary,
    raw = data.frame(
      Usual_Len = len_usual,
      Opt_Len = len_opt,
      Usual_Cov = cov_usual,
      Opt_Cov = cov_opt
    )
  ))
}
```

### Resultados da Simulação

A tabela a seguir consolida as métricas empíricas e teóricas sob ambas as abordagens.

```set.seed(123) # Reprodutibilidade
sim_res <- simular_ic_variancia(sigma2_real = 120, n = 25, B = 10000)

tabela_ic <- sim_res$summary %>%
  gt() %>%
  tab_header(
    title = md("**Comparação de Métodos de Intervalo de Confiança para a Variância**"),
    subtitle = md("Resultados teóricos e empíricos baseados em $B = 10000$ simulações ($n = 25$, $\\sigma^2 = 120$)")
  ) %>%
  cols_label(
    Metodo = "Método",
    Valor_a = "Crítico Inferior (a)",
    Valor_b = "Crítico Superior (b)",
    Cobertura = "Cobertura Empírica",
    Comprimento_Medio = "Comp. Médio Simulado",
    Comprimento_Medio_Teorico = "Comp. Médio Teórico"
  ) %>%
  fmt_number(
    columns = c("Valor_a", "Valor_b", "Comprimento_Medio", "Comprimento_Medio_Teorico"),
    decimals = 5
  ) %>%
  fmt_percent(
    columns = "Cobertura",
    decimals = 2
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(columns = Metodo)
  ) %>%
  opt_align_table_header(align = "left")

tabela_ic
```

---

## 4. Visualização Gráfica do Comportamento dos Intervalos

Construímos um painel 2x2 para fins didáticos, ilustrando as diferenças estruturais e empíricas entre as duas abordagens.

```#| label: fig-painel-ic
#| fig-cap: "Painel de análise comparativa entre as abordagens usual e de comprimento mínimo."
#| fig-width: 12
#| fig-height: 7
#| column: page-right

gerar_painel_ic <- function(sim_res, sigma2_real = 120, n = 25) {
  nu <- n - 1
  conf_level <- 0.95
  alpha <- 1 - conf_level
  
  summary_stats <- sim_res$summary
  raw_data <- sim_res$raw
  
  # Valores críticos
  a_usual <- summary_stats$Valor_a[1]
  b_usual <- summary_stats$Valor_b[1]
  a_opt <- summary_stats$Valor_a[2]
  b_opt <- summary_stats$Valor_b[2]
  
  cores <- c("Usual" = "#2c3e50", "Comprimento Mínimo" = "#e67e22")
  teminha <- theme_minimal() + 
    theme(
      legend.position = "none", 
      plot.title = element_text(size = 11, face = "bold"),
      axis.title = element_text(face = "bold") 
    )
  
  # A. Limites Críticos na Chi-Quadrado (df = 24)
  x_grid <- seq(0, 60, length.out = 1000)
  df_pdf24 <- data.frame(x = x_grid, y = dchisq(x_grid, df = nu))
  
  g_limits <- ggplot(df_pdf24, aes(x = x, y = y)) +
    geom_line(color = "grey30", size = 0.8) +
    geom_vline(aes(xintercept = a_usual, color = "Usual"), linetype = "dashed", size = 0.8) +
    geom_vline(aes(xintercept = b_usual, color = "Usual"), linetype = "dashed", size = 0.8) +
    geom_vline(aes(xintercept = a_opt, color = "Comprimento Mínimo"), linetype = "solid", size = 0.8) +
    geom_vline(aes(xintercept = b_opt, color = "Comprimento Mínimo"), linetype = "solid", size = 0.8) +
    scale_color_manual(values = cores) +
    labs(title = "Limites Críticos na Densidade Chi-quadrado (df=24)", x = "Valor Crítico", y = "Densidade") +
    teminha
  
  # B. Distribuição dos Comprimentos dos Intervalos
  df_len_long <- data.frame(
    Metodo = rep(c("Usual", "Comprimento Mínimo"), each = nrow(raw_data)),
    Comprimento = c(raw_data$Usual_Len, raw_data$Opt_Len)
  )
  
  g_lengths <- ggplot(df_len_long, aes(x = Comprimento, fill = Metodo)) +
    geom_density(alpha = 0.4, color = "white") +
    geom_vline(aes(xintercept = mean(raw_data$Usual_Len), color = "Usual"), linetype = "dashed", size = 1) +
    geom_vline(aes(xintercept = mean(raw_data$Opt_Len), color = "Comprimento Mínimo"), linetype = "solid", size = 1) +
    scale_fill_manual(values = cores) +
    scale_color_manual(values = cores) +
    labs(title = "Distribuição dos Comprimentos dos Intervalos", x = "Comprimento", y = "Densidade") +
    teminha
  
  # C. Condição de Otimalidade (df = 28)
  df_pdf28 <- data.frame(x = x_grid, y = dchisq(x_grid, df = nu + 4))
  y_val <- dchisq(a_opt, df = nu + 4)
  
  g_opt_condition <- ggplot(df_pdf28, aes(x = x, y = y)) +
    geom_line(color = "grey30", size = 0.8) +
    geom_segment(aes(x = a_opt, y = y_val, xend = b_opt, yend = y_val), color = "#e67e22", linetype = "dotted", size = 1) +
    geom_point(aes(x = a_opt, y = y_val), color = "#e67e22", size = 3) +
    geom_point(aes(x = b_opt, y = y_val), color = "#e67e22", size = 3) +
    labs(title = "Condição de Otimalidade: f28(a) = f28(b)", x = "Valor de y", y = "Densidade Chi-quadrado (df=28)") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 11, face = "bold"),
      axis.title = element_text(face = "bold")
    )
  
  # D. Taxa de Cobertura Empírica
  df_cov <- data.frame(
    Metodo = c("Usual", "Comprimento Mínimo"),
    Cobertura = c(mean(raw_data$Usual_Cov), mean(raw_data$Opt_Cov))
  )
  
  g_coverage <- ggplot(df_cov, aes(x = Metodo, y = Cobertura, fill = Metodo)) +
    geom_bar(stat = "identity", width = 0.4, alpha = 0.8) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "red", size = 0.8) +
    coord_cartesian(ylim = c(0.90, 1.0)) +
    scale_fill_manual(values = cores) +
    labs(title = "Taxa de Cobertura Empírica (Nominal = 95%)", x = "Método", y = "Cobertura") +
    teminha
  
  # Organização com Patchwork
  (g_limits + g_lengths) / (g_opt_condition + g_coverage) + 
    plot_layout(guides = 'collect') & 
    theme(legend.position = 'bottom', legend.title = element_blank())
}

# Gerando o painel
gerar_painel_ic(sim_res, sigma2_real = 120, n = 25)
```

---

## 5. Análise e Discussão

A análise comparativa revela conclusões importantes acerca da eficiência metodológica na construção de intervalos de confiança:

1. **Eficiência no Comprimento:** O intervalo de comprimento mínimo ($IC_{\text{opt}}$) reduziu o tamanho do intervalo de confiança de $159.07$ para $148.23$ unidades de comprimento, poupando cerca de $10.84$ unidades. Essa diferença se traduz em uma **redução de $6.82\%$ no comprimento médio** dos intervalos de confiança.
2. **Preservação da Cobertura:** Apesar de o intervalo ótimo ser significativamente menor que o usual, a taxa de cobertura empírica de ambos é rigorosamente igual ao nível nominal estipulado de $95\%$. Isso ocorre porque a otimização atua redistribuindo a probabilidade das caudas de maneira assimétrica para aproveitar a assimetria positiva da distribuição chi-quadrado, sem alterar a massa de probabilidade total ($1-\alpha$) integrada entre os limites críticos.
3. **Distribuição dos Comprimentos:** O gráfico de densidade de comprimentos (Painel B) ilustra claramente que a distribuição do comprimento do $IC_{\text{opt}}$ está deslocada para a esquerda (menores valores) em relação ao $IC_{\text{usual}}$. A variabilidade e a assimetria na distribuição do comprimento são explicadas pelo fato de que o comprimento do intervalo é diretamente proporcional a $S^2$, que é uma variável aleatória.
4. **Justificativa Geométrica da Otimização:** O gráfico da densidade $\chi^2_{28}$ (Painel C) demonstra visualmente a validade da condição de otimalidade $f_{\nu+4}(a) = f_{\nu+4}(b)$. Ao alinhar as ordenadas nesta densidade auxiliar, a alocação de probabilidade nas caudas da distribuição original $\chi^2_{24}$ torna-se desigual: acumula-se menos probabilidade na cauda esquerda e mais na cauda direita. Isso desloca o intervalo ligeiramente para a direita (limites maiores na escala de $\chi^2$, o que resulta em limites menores e mais próximos na escala inversa de $\sigma^2$), encurtando o intervalo sem sacrificar a confiança.
