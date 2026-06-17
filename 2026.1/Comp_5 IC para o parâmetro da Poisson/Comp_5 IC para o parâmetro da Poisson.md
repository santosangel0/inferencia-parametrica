
```#| label: setup-libraries
#| message: false
#| warning: false
library(gt)
library(tibble)
```

## 1. Enunciado do Exercício

Um hospital universitário deseja estudar a demanda noturna em seu pronto atendimento. Seja $X_i$ o número de pacientes que chegam ao pronto atendimento entre 00h e 1h na $i$-ésima noite. Suponha que $X_1, \dots, X_n$ são independentes e identicamente distribuídos (iid) com distribuição $\text{Poisson}(\lambda)$, onde $\lambda$ representa o número médio de chegadas por hora nesse período.

O hospital coletou dados durante $n = 20$ noites e obteve os seguintes números de pacientes:
$$3, 2, 4, 1, 5, 3, 2, 4, 6, 3, 2, 1, 5, 4, 3, 2, 4, 3, 1, 2$$

## 2. Resolução dos Itens

Abaixo realizamos a importação dos dados no R e a resolução passo a passo.

```#| label: setup-dados
#| message: false
#| warning: false
# Inserindo os dados como vetor conforme os requisitos
dados <- c(3, 2, 4, 1, 5, 3, 2, 4, 6, 3, 2, 1, 5, 4, 3, 2, 4, 3, 1, 2)
n <- length(dados)
T_total <- sum(dados)
```

### a) Estimativa de Máxima Verossimilhança (EMV)

Para uma amostra iid de uma distribuição de Poisson, a função de verossimilhança é dada por:
$$L(\lambda) = \prod_{i=1}^n \frac{e^{-\lambda} \lambda^{x_i}}{x_i!} = \frac{e^{-n\lambda} \lambda^{\sum x_i}}{\prod x_i!}$$

A log-verossimilhança é:
$$\ln L(\lambda) = -n\lambda + \left(\sum_{i=1}^n x_i\right) \ln(\lambda) - \sum_{i=1}^n \ln(x_i!)$$

Derivando em relação a $\lambda$ e igualando a zero, encontramos o Estimador de Máxima Verossimilhança (EMV):
$$\frac{d}{d\lambda} \ln L(\lambda) = -n + \frac{\sum x_i}{\lambda} = 0 \implies \hat{\lambda} = \bar{X} = \frac{T}{n}$$

Calculando o valor numérico em R:

```#| label: emv
#| message: false
#| warning: false
lambda_hat <- mean(dados)
cat("O Estimador de Máxima Verossimilhança (EMV) para lambda é:", lambda_hat)
```

Portanto, a estimativa pontual para a taxa média de chegadas de pacientes é de **`lambda_hat` pacientes por hora**.

### b) Intervalo de Confiança Assintótico de 95% (Wald)

Pela normalidade assintótica do EMV, temos que para um tamanho de amostra $n$ suficientemente grande:
$$\hat{\lambda} \sim N\left(\lambda, \frac{\lambda}{n}\right)$$

Substituindo a variância desconhecida por sua estimativa baseada no EMV, obtemos o erro padrão de $\hat{\lambda}$ como $\sqrt{\hat{\lambda}/n}$. O intervalo de confiança de $95\%$ é construído como:
$$\text{IC}_{\text{Wald}}(\lambda) = \hat{\lambda} \pm z_{1-\alpha/2} \sqrt{\frac{\hat{\lambda}}{n}}$$

Onde $z_{1-\alpha/2} \approx 1,96$ para um nível de significância de $5\%$ (confiança de 95%).

Calculando em R:

```#| label: ic-wald
#| message: false
#| warning: false
alpha <- 0.05
z_val <- qnorm(1 - alpha/2)

se_wald <- sqrt(lambda_hat / n)
ic_wald_inf <- lambda_hat - z_val * se_wald
ic_wald_sup <- lambda_hat + z_val * se_wald

cat(sprintf("IC Wald (95%%): [%.4f, %.4f]\n", ic_wald_inf, ic_wald_sup))
```

### c) Intervalo de Confiança via Transformação Estabilizadora de Variância

A variância da distribuição de Poisson depende diretamente do parâmetro $\lambda$. Para estabilizar a variância, aplicamos a transformação $g(t) = \sqrt{t}$. 

Pelo Método Delta, a distribuição assintótica de $\sqrt{\hat{\lambda}}$ é:
$$\sqrt{\hat{\lambda}} \sim N\left(\sqrt{\lambda}, \frac{1}{4n}\right)$$

A variância assintótica torna-se independente de $\lambda$. Assim, podemos obter um intervalo de confiança assintótico de 95% para $\sqrt{\lambda}$:
$$\text{IC}(\sqrt{\lambda}) = \sqrt{\hat{\lambda}} \pm z_{1-\alpha/2} \frac{1}{2\sqrt{n}}$$

Elevando ao quadrado os limites inferior e superior deste intervalo, obtemos o intervalo final para $\lambda$:
$$\text{IC}_{\text{Estabilizada}}(\lambda) = \left[ \left(\sqrt{\hat{\lambda}} - z_{1-\alpha/2} \frac{1}{2\sqrt{n}}\right)^2 \,\, , \,\, \left(\sqrt{\hat{\lambda}} + z_{1-\alpha/2} \frac{1}{2\sqrt{n}}\right)^2 \right]$$

Calculando em R:

```#| label: ic-estabilizadora
#| message: false
#| warning: false
lim_trans_inf <- sqrt(lambda_hat) - z_val * (1 / (2 * sqrt(n)))
lim_trans_sup <- sqrt(lambda_hat) + z_val * (1 / (2 * sqrt(n)))

ic_est_inf <- lim_trans_inf^2
ic_est_sup <- lim_trans_sup^2

cat(sprintf("IC Estabilizador de Variância (95%%): [%.4f, %.4f]\n", ic_est_inf, ic_est_sup))
```

### d) Intervalo de Confiança Exato de 95%

Utilizando a estatística suficiente $T = \sum_{i=1}^n X_i$, sabemos que $T \sim \text{Poisson}(n\lambda)$. A relação entre a cauda da distribuição de Poisson e os quantis da distribuição Qui-Quadrado nos permite expressar o intervalo de confiança exato de 95% de forma fechada como:
$$\text{IC}_{\text{Exato}}(\lambda) = \left[ \frac{\chi^2_{2T, \alpha/2}}{2n} \,\, , \,\, \frac{\chi^2_{2(T+1), 1-\alpha/2}}{2n} \right]$$

Onde:
* $T = \sum x_i$ é a soma observada de pacientes.
* $\chi^2_{\nu, p}$ representa o quantil $p$ da distribuição Qui-Quadrado com $\nu$ graus de liberdade.

Calculando em R:

```#| label: ic-exato
#| message: false
#| warning: false
ic_exato_inf <- qchisq(alpha/2, df = 2 * T_total) / (2 * n)
ic_exato_sup <- qchisq(1 - alpha/2, df = 2 * (T_total + 1)) / (2 * n)

cat(sprintf("IC Exato (95%%): [%.4f, %.4f]\n", ic_exato_inf, ic_exato_sup))
```

### e) Comparação dos Intervalos Obtidos

A tabela abaixo compara os limites inferior, superior e a amplitude total de cada um dos métodos propostos:

```#| label: comparacao-tabela
#| message: false
#| warning: false

metodos <- c("Assintótico (Wald)", "Estabilizador de Variância", "Exato (Qui-Quadrado)")
limites_inf <- c(ic_wald_inf, ic_est_inf, ic_exato_inf)
limites_sup <- c(ic_wald_sup, ic_est_sup, ic_exato_sup)
amplitudes <- limites_sup - limites_inf

tabela_comp <- data.frame(
  Método = metodos,
  `Limite Inferior` = round(limites_inf, 4),
  `Limite Superior` = round(limites_sup, 4),
  Amplitude = round(amplitudes, 4),
  check.names = FALSE
)

tabela_comp %>%
  gt() %>%
  tab_header(
    title = md("**Comparação dos Intervalos de Confiança**"),
    subtitle = md("Limites e amplitudes para cada método")
  ) %>%
  fmt_number(
    columns = -Método,
    decimals = 4
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(columns = Método)
  ) %>%
  opt_align_table_header(align = "left")
```

### Análise Comparativa

* **Intervalo Assintótico (Wald)**: Apresenta limites ligeiramente deslocados para a esquerda em comparação com o exato. O método assume normalidade da média amostral e erro padrão fixado na estimativa pontual. Por ser uma amostra pequena ($n = 20$), esta aproximação pode sofrer alguma distorção (subestimando a variabilidade real na cauda direita).
* **Intervalo com Transformação Estabilizadora**: Corrige a dependência entre a variância e a média por meio da transformação $\sqrt{X}$. Apresenta limites muito próximos ao do intervalo exato, demonstrando ser uma excelente alternativa para corrigir assimetrias comuns na distribuição de Poisson quando o tamanho amostral é moderado.
* **Intervalo Exato**: Fornece a cobertura exata de 95% garantida teoricamente para qualquer tamanho amostral $n$, sendo ligeiramente mais conservador (maior amplitude) e assimétrico à direita para acomodar a natureza assimétrica da distribuição de Poisson subjacente.
