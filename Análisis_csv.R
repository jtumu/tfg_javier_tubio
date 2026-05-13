#### Carga de datos: cambia las ubicaciones a donde estén tus datos

Ecoli <- read.delim("EColiK12_tabla.tsv")

####  Resultados FANTASIA

EColi_random100 <- read.csv("~/Desktop/TFG/Result_fantasia_Lite/E_Coli_random_100_resultados/EColi_random_100.csv")
EColi_random200 <- read.csv("~/Desktop/TFG/Result_fantasia_Lite/E_Coli_random_200_resultados/EColi_random_200.csv")
EColi_random300 <- read.csv("~/Desktop/TFG/Result_fantasia_Lite/E_Coli_random_300_resultados/EColi_random_300.csv")
EColi_random400 <- read.csv("~/Desktop/TFG/Result_fantasia_Lite/E_Coli_random_400_resultados/EColi_random_400.csv")
EColi_random500 <- read.csv("~/Desktop/TFG/Result_fantasia_Lite/E_Coli_random_500_resultados/EColi_random_500.csv")
EColi_control <- read.csv("~/Desktop/TFG/Result_fantasia_Lite/resultados_fantasiaLite_control/results_control_EColi/results_control_EColi.csv")

## Nota: quitar los # para instalar paquetes

#install.packages("tidyverse")
#install.packages("tidytext")
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# if (!require("GOSemSim", quietly = TRUE)) BiocManager::install("GOSemSim")
# if (!require("org.EcK12.eg.db", quietly = TRUE)) BiocManager::install("org.EcK12.eg.db")
# if (!require("UniProt.ws", quietly = TRUE)) BiocManager::install("UniProt.ws")
library(UniProt.ws)
library(tidyverse)
library(tidytext)
library(GOSemSim)
library(org.EcK12.eg.db)
library(GO.db)

# Limpiar el ID de Fantasia
EColi_random100 <- EColi_random100 %>%
  mutate(Entry = sapply(strsplit(query_accession, "\\|"), `[`, 2))
EColi_random200 <- EColi_random200 %>%
  mutate(Entry = sapply(strsplit(query_accession, "\\|"), `[`, 2))
EColi_random300 <- EColi_random300 %>%
  mutate(Entry = sapply(strsplit(query_accession, "\\|"), `[`, 2))
EColi_random400 <- EColi_random400 %>%
  mutate(Entry = sapply(strsplit(query_accession, "\\|"), `[`, 2))
EColi_random500 <- EColi_random500 %>%
  mutate(Entry = sapply(strsplit(query_accession, "\\|"), `[`, 2))
EColi_control <- EColi_control %>%
  mutate(Entry = sapply(strsplit(query_accession, "\\|"), `[`, 2))

# Dividir en filas los datos de UniProt 
EColi_largo <- Ecoli %>%
  dplyr::select(Entry, `Gene.Ontology..GO.`) %>% 
  separate_rows(`Gene.Ontology..GO.`, sep = "; ") %>%
  mutate(go_id_real = str_extract(`Gene.Ontology..GO.`, "GO:[0-9]+")) %>%
  filter(!is.na(go_id_real))

#### Análisis del Control

##  Control VS UNIPROT
evaluacion_control <- EColi_control %>%
  dplyr::select(model_key, Entry, go_id) %>%
  distinct() %>%
  left_join(EColi_largo, by = "Entry", relationship = "many-to-many") %>%
  mutate(coincide = (go_id == go_id_real))

resumen_control_exacto <- evaluacion_control %>%
  group_by(model_key) %>%
  summarise(
    total_predicciones = n(),
    aciertos_exactos = sum(coincide, na.rm = TRUE),
    precision_control = (aciertos_exactos / total_predicciones) * 100
  )

print(resumen_control_exacto)

#### GO TERMS ####
# Contar frecuencia de predicciones GO en los proteomas random
top_funciones_random_100 <- EColi_random100 %>%
  group_by(model_key, go_id, go_description) %>%
  summarise(veces_predicha = n(), .groups = "drop") %>%
  arrange(model_key, desc(veces_predicha))

top_funciones_random_200 <- EColi_random200 %>%
  group_by(model_key, go_id, go_description) %>%
  summarise(veces_predicha = n(), .groups = "drop") %>%
  arrange(model_key, desc(veces_predicha))

top_funciones_random_300 <- EColi_random300 %>%
  group_by(model_key, go_id, go_description) %>%
  summarise(veces_predicha = n(), .groups = "drop") %>%
  arrange(model_key, desc(veces_predicha))

top_funciones_random_400 <- EColi_random400 %>%
  group_by(model_key, go_id, go_description) %>%
  summarise(veces_predicha = n(), .groups = "drop") %>%
  arrange(model_key, desc(veces_predicha))

top_funciones_random_500 <- EColi_random500 %>%
  group_by(model_key, go_id, go_description) %>%
  summarise(veces_predicha = n(), .groups = "drop") %>%
  arrange(model_key, desc(veces_predicha))

top_funciones_control <- EColi_control %>% 
  group_by(model_key, go_id, go_description) %>%
  summarise(veces_predicha = n(), .groups = "drop") %>%
  arrange(model_key, desc(veces_predicha))

todos_los_datos_random <- bind_rows(
  mutate(top_funciones_random_100, randomizacion = "R100"),
  mutate(top_funciones_random_200, randomizacion = "R200"),
  mutate(top_funciones_random_300, randomizacion = "R300"),
  mutate(top_funciones_random_400, randomizacion = "R400"),
  mutate(top_funciones_random_500, randomizacion = "R500")
)
resumen_tabla_full <- todos_los_datos_random %>%
  group_by(model_key, go_id, go_description) %>%
  summarise(
    apariciones_replicas = n(),               
    media_veces_predicha = mean(veces_predicha), 
    .groups = "drop") %>% 
  arrange(model_key,  desc(media_veces_predicha))

top_ankh3_control <- top_funciones_control %>% filter(model_key == "Ankh3-Large") %>% head(10)
top_prott5_control <- top_funciones_control %>% filter(model_key == "Prot-T5") %>% head(10)
todos_los_tops_control <- bind_rows(top_ankh3_control,top_prott5_control)



## GO TERMS ORIGINAL
EColi_largo_limpio <- EColi_largo %>%
  mutate(
    go_description_real = str_trim(str_remove(`Gene.Ontology..GO.`, "\\[GO:[0-9]+\\]"))
  )

# Veces que aparece cada GO en el proteoma original
frecuencias_reales <- EColi_largo_limpio %>%
  group_by(go_id_real, go_description_real) %>%
  summarise(veces_real = n(), .groups = "drop") %>%
  arrange(desc(veces_real)) 

#  Top 10 real
top10_real <- frecuencias_reales %>% head(10)
print(top10_real)

grafico_real <- ggplot(top10_real, aes(
  x = reorder(go_description_real, veces_real), 
  y = veces_real, 
  fill = veces_real
)) +
  geom_col(color = "black", linewidth = 0.2) + 
  geom_text(aes(label = veces_real), hjust = -0.15, size = 7, color = "black") + 
  coord_flip() + 
  scale_y_continuous(limits = c(0, max(top10_real$veces_real) * 1.15)) + 
  scale_fill_viridis_c(option = "mako", direction = -1) +
  labs(
    title = "Top 10 términos GO en el proteoma real de E. Coli",
    subtitle = "Distribución natural según UniProt",
    x = "Término GO",
    y = "Número de proteínas"
  ) +
  theme_minimal(base_size = 25) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "none"
  )

# Gráfico
print(grafico_real)


#### GRÁFICO GO TERMS ####

datos_grafico <- resumen_tabla_full %>% 
  group_by(model_key) %>%
  slice_max(order_by = media_veces_predicha, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(go_wrap = str_wrap(go_description, width = 35))

grafico_randomizados_ordenado <- ggplot(
  datos_grafico,
  aes(
    x = tidytext::reorder_within(go_wrap, media_veces_predicha, model_key), 
    y = media_veces_predicha, 
    fill = as.numeric(apariciones_replicas))) +
  geom_col(color = "black", linewidth = 0.3) +
  geom_text(aes(label = round(media_veces_predicha, 1)), 
            hjust = -0.2, size = 5) +
  coord_flip() +
  facet_wrap(~ model_key, scales = "free_y") + 
  tidytext::scale_x_reordered() +
  scale_fill_viridis_c(
    option = "plasma", 
    direction = 1,
    name = "Nº de\nréplicas",
    breaks = 1:5,
    limits = c(1, 5)
  ) +
  labs(
    title = paste0("Top 10 términos GO en secuencias aleatorizadas"),
    x = "Término GO asignado",
    y = "Número de proteínas (Media)"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  theme_minimal(base_size = 18) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(grafico_randomizados_ordenado)

datos_grafico_control <- todos_los_tops_control%>% 
  mutate(
    score_orden = veces_predicha
  ) %>%
  group_by(model_key) %>%
  slice_max(order_by = score_orden, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(go_description_wrap = str_wrap(go_description, width = 35))

grafico_control_ordenado <- ggplot(datos_grafico_control, aes(
  x = reorder_within(go_description_wrap, veces_predicha, model_key), 
  y = veces_predicha, 
  fill = veces_predicha )) +
  geom_col(color = "black", linewidth = 0.2) + 
  geom_text(aes(label = veces_predicha), 
            hjust = -0.2, size = 5, color = "black") + 
  coord_flip() + 
  facet_wrap(~ model_key, scales = "free_y") + 
  scale_x_reordered() + 
  scale_fill_viridis_c(option = "plasma", name = "Nº Proteínas") +
  labs(
    title = "Top 10 términos GO en el control E. Coli (proteoma original)",
    x = "Término GO asignado",
    y = "Número de proteínas") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"), 
    plot.subtitle = element_text(hjust = 0.5),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold"),
    legend.position = "right",
    panel.grid.minor = element_blank())
  
print(grafico_control_ordenado)

#### H2: ¿Existe memorización de la composición de aa? ####
todas_las_randomizaciones <- bind_rows(
  mutate(EColi_random100, randomizacion = "R100"),
  mutate(EColi_random200, randomizacion = "R200"),
  mutate(EColi_random300, randomizacion = "R300"),
  mutate(EColi_random400, randomizacion = "R400"),
  mutate(EColi_random500, randomizacion = "R500")
)

predicciones_ia <- todas_las_randomizaciones %>%
  dplyr::select(randomizacion, model_key, Entry, go_id) %>%
  distinct()

# Comparar con términos reales de UniProt
evaluacion_h2 <- predicciones_ia %>%
  left_join(EColi_largo, by = "Entry", relationship = "many-to-many") %>%
  mutate(coincide = (go_id == go_id_real)) %>%
  group_by(randomizacion, model_key, Entry, go_id) %>%
  summarise(mantuvo_funcion_original = any(coincide %in% TRUE), .groups = "drop")

# Porcentajes por cada random y modelo
resultados_h2_por_replica <- evaluacion_h2 %>%
  group_by(randomizacion, model_key) %>%
  summarise(
    total_predicciones = n(),
    funciones_mantenidas = sum(mantuvo_funcion_original),
    porcentaje_mantenido = (funciones_mantenidas / total_predicciones) * 100,
    .groups = "drop"
  )

# Resumen final
resumen_h2_final <- resultados_h2_por_replica %>%
  group_by(model_key) %>%
  summarise(
    media_porcentaje = mean(porcentaje_mantenido),
    desviacion_estandar = sd(porcentaje_mantenido)
  )

print(resumen_h2_final)

comparativa_h2_completa <- resumen_h2_final %>% 
  mutate(Tipo = "Secuencias Random") %>%
  bind_rows(
    resumen_control_exacto %>% 
      dplyr::select(model_key, media_porcentaje = precision_control) %>% 
      mutate(Tipo = "Control")
  )
print(comparativa_h2_completa)


#### ORGANISMOS APORTADORES DE FUNCIÓN ####
up <- UniProt.ws(taxId=9606)
obtener_organismos <- function(df) {
  ids_unicos <- unique(df$uniprot_accession)
  
  mapping <- select(up, 
                    keys = ids_unicos, 
                    columns = c("organism_name"), 
                    keytype = "UniProtKB")
  
  df_anotado <- df %>%
    left_join(mapping, by = c("uniprot_accession" = "Entry"))
  
  return(df_anotado)
}

# Control
Orgs_Control<-obtener_organismos(EColi_control)
conteo_control <- Orgs_Control %>%
  count(Organism) %>% 
  arrange(desc(n))
conteo_control <- Orgs_Control %>%
  count(model_key, Organism) %>% 
  group_by(model_key) %>%
  mutate(porcentaje = (n / sum(n)) * 100) %>%
  arrange(model_key, desc(porcentaje)) %>%
  group_by(model_key) %>%
  slice_max(order_by = porcentaje, n = 10, with_ties = FALSE)

 view(conteo_control)

# Randoms 
Orgs_EColi100<-obtener_organismos(EColi_random100)
Orgs_EColi200<-obtener_organismos(EColi_random200)
Orgs_EColi300<-obtener_organismos(EColi_random300)
Orgs_EColi400<-obtener_organismos(EColi_random400)
Orgs_EColi500<-obtener_organismos(EColi_random500)


conteo_Ecoli_100 <- Orgs_EColi100 %>%
  count(model_key, Organism) %>% 
  group_by(model_key) %>%
  mutate(porcentaje = (n / sum(n)) * 100) %>%
  ungroup()
#view(conteo_Ecoli_100)

conteo_Ecoli_200 <- Orgs_EColi200 %>%
  count(model_key, Organism) %>% 
  group_by(model_key) %>%
  mutate(porcentaje = (n / sum(n)) * 100) %>%
  ungroup()
# view(conteo_Ecoli_200)

conteo_Ecoli_300 <- Orgs_EColi300 %>%
  count(model_key, Organism) %>% 
  group_by(model_key) %>%
  mutate(porcentaje = (n / sum(n)) * 100) %>%
  ungroup()

conteo_Ecoli_400 <- Orgs_EColi400 %>%
  count(model_key, Organism) %>% 
  group_by(model_key) %>%
  mutate(porcentaje = (n / sum(n)) * 100) %>%
  ungroup()

conteo_Ecoli_500 <- Orgs_EColi500 %>%
  count(model_key, Organism) %>% 
  group_by(model_key) %>%
  mutate(porcentaje = (n / sum(n)) * 100) %>%
  ungroup()

todas_las_org_randoms <- bind_rows(
  conteo_Ecoli_100 %>% mutate(Replica = "R1"),
  conteo_Ecoli_200 %>% mutate(Replica = "R2"),
  conteo_Ecoli_300 %>% mutate(Replica = "R3"),
  conteo_Ecoli_400 %>% mutate(Replica = "R4"),
  conteo_Ecoli_500 %>% mutate(Replica = "R5")
) %>% 
  dplyr::rename(conteo = n)


resumen_org_randoms <- todas_las_org_randoms %>%
  group_by(model_key, Organism) %>%
  summarise(
    n_medio = mean(conteo),
    sd_n = sd(conteo, na.rm = TRUE),
    media_porcentaje = mean(porcentaje),
    sd_porcentaje = sd(porcentaje, na.rm = TRUE),
    apariciones_replicas = n(),
    .groups = "drop"
  ) %>%
  mutate(sd_porcentaje = replace_na(sd_porcentaje, 0)) %>%
  
  group_by(model_key) %>%
  slice_max(order_by = media_porcentaje, n = 10, with_ties = FALSE) %>%
  ungroup()

view(resumen_org_randoms)




#### DISTANCIAS GO ####


########################
#Función Molecular (MF)#
########################
go_data_mf <- godata('org.EcK12.eg.db', ont="MF", computeIC=FALSE)

#  Filtrar FANTASIA, quedarnos solo con Función Molecular ("F")
predicciones_mf <- todas_las_randomizaciones %>%
  filter(category == "F") %>% 
  dplyr::select(randomizacion, model_key, Entry, go_id_pred = go_id) %>%
  distinct()

predicciones_mf_control <- EColi_control %>%
  filter(category == "F") %>% 
  dplyr::select( model_key, Entry, go_id_pred = go_id) %>%
  distinct()

#  Filtrar el UniProt para MF
reales_dist_go <- EColi_largo %>%
  dplyr::select(Entry, go_id_real) %>%
  distinct()

# Listas de GOs por proteína para  compararlas 
lista_predicciones <- predicciones_mf %>%
  group_by(randomizacion, model_key, Entry) %>%
  summarise(go_preds = list(go_id_pred), .groups = "drop")

lista_predicciones_mf_control <- predicciones_mf_control %>% 
  group_by(model_key, Entry) %>%
  summarise(go_preds = list(go_id_pred), .groups = "drop")

lista_reales <- reales_dist_go %>%
  group_by(Entry) %>%
  summarise(go_reals = list(go_id_real), .groups = "drop")

datos_similitud_mf <- lista_predicciones %>%
  inner_join(lista_reales, by = "Entry")
datos_similitud_mf_control <- lista_predicciones_mf_control %>% 
  inner_join(lista_reales, by = "Entry")


#  map2_dbl -> función mgoSim. Compara dos listas de términos GO. 
datos_similitud_mf <- datos_similitud_mf %>%
  mutate(
    similitud_wang = map2_dbl(go_preds, go_reals, function(preds, reals) {
      sim <- tryCatch(
        mgoSim(preds, reals, semData = go_data_mf, measure = "Wang", combine = "BMA"),
        error = function(e) NA
      )
      return(sim)
    })
  )

datos_similitud_mf_control <- datos_similitud_mf_control %>% 
  mutate(
    similitud_wang = map2_dbl(go_preds, go_reals, function(preds, reals) {
      sim <- tryCatch(
        mgoSim(preds, reals, semData = go_data_mf, measure = "Wang", combine = "BMA"),
        error = function(e) NA
      )
      return(sim)
    })
  )

# Resumen de resultados (quitar los NAs)
resumen_similitud_mf <- datos_similitud_mf %>%
  filter(!is.na(similitud_wang)) %>%
  group_by(model_key) %>%
  summarise(
    similitud_media = mean(similitud_wang),
    mediana = median(similitud_wang),
    desviacion = sd(similitud_wang),
    n_proteinas = n()
  )
print(resumen_similitud_mf)


resumen_similitud_mf_control<- datos_similitud_mf_control %>% 
  filter(!is.na(similitud_wang)) %>%
  group_by(model_key) %>%
  summarise(
    similitud_media = mean(similitud_wang),
    mediana = median(similitud_wang),
    desviacion = sd(similitud_wang),
    n_proteinas = n()
  )

print(resumen_similitud_mf_control)


#########################
#Biological process (BP)#
#########################

go_data_bp <- godata('org.EcK12.eg.db', ont="BP", computeIC=FALSE)

predicciones_bp <- todas_las_randomizaciones %>%
  filter(category == "P") %>% 
  dplyr::select(randomizacion, model_key, Entry, go_id_pred = go_id) %>%
  distinct()

predicciones_bp_control <- EColi_control %>% 
  filter(category == "P") %>% 
  dplyr::select(model_key, Entry, go_id_pred = go_id) %>% 
  distinct()

lista_predicciones_bp <- predicciones_bp %>%
  group_by(randomizacion, model_key, Entry) %>%
  summarise(go_preds = list(go_id_pred), .groups = "drop")

lista_predicciones_bp_control <- predicciones_bp_control %>% 
  group_by(model_key, Entry) %>% 
  summarise(go_preds = list(go_id_pred), .groups = "drop")

datos_similitud_bp <- lista_predicciones_bp %>%
  inner_join(lista_reales, by = "Entry")

datos_similitud_bp_control <- lista_predicciones_bp_control %>% 
  inner_join(lista_reales, by = "Entry")

datos_similitud_bp <- datos_similitud_bp %>%
  mutate(
    similitud_wang = map2_dbl(go_preds, go_reals, function(preds, reals) {
      sim <- tryCatch(
        mgoSim(preds, reals, semData = go_data_bp, measure = "Wang", combine = "BMA"),
        error = function(e) NA
      )
      return(sim)
    })
  )

datos_similitud_bp_control <- datos_similitud_bp_control %>% 
  mutate(
    similitud_wang = map2_dbl(go_preds, go_reals, function(preds, reals) {
      sim <- tryCatch(
        mgoSim(preds, reals, semData = go_data_bp, measure = "Wang", combine = "BMA"),
        error = function(e) NA
      )
      return(sim)
    })
  )

resumen_similitud_bp <- datos_similitud_bp %>%
  filter(!is.na(similitud_wang)) %>%
  group_by(model_key) %>%
  summarise(
    similitud_media = mean(similitud_wang),
    mediana = median(similitud_wang),
    desviacion = sd(similitud_wang),
    n_proteinas = n())

print(resumen_similitud_bp)

resumen_similitud_bp_control <- datos_similitud_bp_control %>% 
  filter(!is.na(similitud_wang)) %>%
  group_by(model_key) %>%
  summarise(
    similitud_media = mean(similitud_wang),
    mediana = median(similitud_wang),
    desviacion = sd(similitud_wang),
    n_proteinas = n())
print(resumen_similitud_bp_control)

########################
#Celular component (CC)#
########################

go_data_cc <- godata('org.EcK12.eg.db', ont="CC", computeIC=FALSE)

predicciones_cc <- todas_las_randomizaciones %>%
  filter(category == "C") %>% 
  dplyr::select(randomizacion, model_key, Entry, go_id_pred = go_id) %>%
  distinct()

predicciones_cc_control <- EColi_control %>%
  filter(category == "C") %>% 
  dplyr::select(model_key, Entry, go_id_pred = go_id) %>%
  distinct()

lista_predicciones_cc <- predicciones_cc %>%
  group_by(randomizacion, model_key, Entry) %>%
  summarise(go_preds = list(go_id_pred), .groups = "drop")

lista_predicciones_cc_control <- predicciones_cc_control %>%
  group_by( model_key, Entry) %>%
  summarise(go_preds = list(go_id_pred), .groups = "drop")

datos_similitud_cc <- lista_predicciones_cc %>%
  inner_join(lista_reales, by = "Entry")

datos_similitud_cc_control <- lista_predicciones_cc_control %>%
  inner_join(lista_reales, by = "Entry")

datos_similitud_cc <- datos_similitud_cc %>%
  mutate(
    similitud_wang = map2_dbl(go_preds, go_reals, function(preds, reals) {
      sim <- tryCatch(
        mgoSim(preds, reals, semData = go_data_cc, measure = "Wang", combine = "BMA"),
        error = function(e) NA
      )
      return(sim)
    })
  )

datos_similitud_cc_control <- datos_similitud_cc_control %>%
  mutate(
    similitud_wang = map2_dbl(go_preds, go_reals, function(preds, reals) {
      sim <- tryCatch(
        mgoSim(preds, reals, semData = go_data_cc, measure = "Wang", combine = "BMA"),
        error = function(e) NA
      )
      return(sim)
    })
  )

resumen_similitud_cc <- datos_similitud_cc %>%
  filter(!is.na(similitud_wang)) %>%
  group_by(model_key) %>%
  summarise(
    similitud_media = mean(similitud_wang),
    mediana = median(similitud_wang),
    desviacion = sd(similitud_wang),
    n_proteinas = n()
  )

print(resumen_similitud_cc)

resumen_similitud_cc_control <- datos_similitud_cc_control %>%
  filter(!is.na(similitud_wang)) %>%
  group_by(model_key) %>%
  summarise(
    similitud_media = mean(similitud_wang),
    mediana = median(similitud_wang),
    desviacion = sd(similitud_wang),
    n_proteinas = n()
  )

print(resumen_similitud_cc_control)


######## TODAS LAS DISTANCIAS ###########
distancias_similitud <- bind_rows(resumen_similitud_mf,resumen_similitud_bp,resumen_similitud_cc)
distancias_similitud <- mutate(distancias_similitud, GO = c("MF","MF","BP","BP","CC","CC"))
distancias_similitud

distancias_similitud_control <- bind_rows(resumen_similitud_mf_control,resumen_similitud_bp_control, resumen_similitud_cc_control) %>% 
  mutate( GO = c("MF","MF","BP","BP","CC","CC"))
distancias_similitud_control


#### EnrichmentGO ####

total_reales <- sum(frecuencias_reales$veces_real)

datos_enriquecimiento <- todas_las_randomizaciones %>%
  group_by(model_key, go_id, go_description) %>%
  summarise(veces_pred_media = n() / 5, .groups = "drop") %>%
  left_join(frecuencias_reales, by = c("go_id" = "go_id_real")) %>%
  mutate(veces_real = replace_na(veces_real, 0)) 

datos_enriquecimiento_control <- EColi_control %>% 
  group_by(model_key, go_id, go_description) %>%
  summarise(veces_pred_media = n() / 5, .groups = "drop") %>%
  left_join(frecuencias_reales, by = c("go_id" = "go_id_real")) %>%
  mutate(veces_real = replace_na(veces_real, 0)) 


datos_enriquecimiento <- datos_enriquecimiento %>%
  group_by(model_key) %>%
  mutate(total_pred = sum(veces_pred_media)) %>%
  ungroup() %>%
  mutate(
    Protein_Ratio = veces_pred_media / total_pred,
    Fold_Enrichment = (veces_pred_media / total_pred) / ((veces_real + 1) / (total_reales + 1)),
    p_value = phyper(round(veces_pred_media) - 1, veces_real, total_reales - veces_real, round(total_pred), lower.tail = FALSE),
    p_adjust = p.adjust(p_value, method = "fdr")
  ) %>%
  filter(veces_pred_media >= 10)

datos_enriquecimiento_control <- datos_enriquecimiento_control %>% 
  group_by(model_key) %>% 
  mutate(total_pred = sum(veces_pred_media)) %>%
  ungroup() %>%
  mutate(
    Protein_Ratio = veces_pred_media / total_pred,
    Fold_Enrichment = (veces_pred_media / total_pred) / ((veces_real + 1) / (total_reales + 1)),
    p_value = phyper(round(veces_pred_media) - 1, veces_real, total_reales - veces_real, round(total_pred), lower.tail = FALSE),
    p_adjust = p.adjust(p_value, method = "fdr")
  ) %>%
  filter(veces_pred_media >= 10)

#Selección del TOP 10 
top_enriquecidos <- datos_enriquecimiento %>%
  group_by(model_key) %>%
  slice_max(order_by = Protein_Ratio, n = 10, with_ties = FALSE) %>%
  ungroup() %>% 
  mutate(go_description_wrap = str_wrap(go_description, width = 35))

top_enriquecidos_control <- datos_enriquecimiento_control %>%
  group_by(model_key) %>%
  slice_max(order_by = Protein_Ratio, n = 10, with_ties = FALSE) %>%
  ungroup() %>% 
  mutate(go_description_wrap = str_wrap(go_description, width = 35))

grafico_enriquecimiento <- ggplot(top_enriquecidos, 
                                  aes(x = Protein_Ratio, 
                                      y = reorder_within(go_description, Protein_Ratio, model_key), 
                                      size = veces_pred_media, 
                                      color = p_adjust)) +
  geom_point(alpha = 0.8) +
  facet_wrap(~ model_key, scales = "free_y") + 
  scale_y_reordered() + 
  scale_color_gradient(low = "firebrick", high = "dodgerblue", name = "p.adjust") +
  scale_size_continuous(name = "Count") +
  labs(
    title = "Enriquecimiento GO en secuencias aleatorizadas",
    subtitle = "Términos sobre-predichos por la IA frente al proteoma real de E. Coli (UniProt)",
    x = "Protein Ratio", 
    y = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    plot.subtitle = element_text(hjust = 0.5, size = 11, margin = margin(b=15)),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold"),
    panel.grid.major.y = element_line(color = "grey90", linetype = "dashed")
  )

print(grafico_enriquecimiento)

grafico_enriquecimiento_control <- ggplot(top_enriquecidos_control, 
                                  aes(x = Protein_Ratio, 
                                      y = reorder_within(go_description, Protein_Ratio, model_key), 
                                      size = veces_pred_media, 
                                      color = p_adjust)) +
  geom_point(alpha = 0.8) +
  facet_wrap(~ model_key, scales = "free_y") + 
  scale_y_reordered() + 
  scale_color_gradient(low = "firebrick", high = "dodgerblue", name = "p.adjust") +
  scale_size_continuous(name = "Count") +
  labs(
    title = "Enriquecimiento GO control",
    subtitle = "Términos sobre-predichos por la IA frente al proteoma real de E. Coli (UniProt)",
    x = "Protein Ratio", 
    y = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    plot.subtitle = element_text(hjust = 0.5, size = 11, margin = margin(b=15)),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold"),
    panel.grid.major.y = element_line(color = "grey90", linetype = "dashed")
  )

print(grafico_enriquecimiento_control)
#### IC

#########################
#    CALCULO IC MF      #
#########################

go_data_mf_ic <- godata('org.EcK12.eg.db', ont="MF", computeIC=TRUE)
vector_ic_mf <- go_data_mf_ic@IC

# Calcular el IC medio de E. Coli

ic_real_mf <- reales_dist_go %>%
  mutate(ic_value = vector_ic_mf[go_id_real]) %>%
  filter(!is.na(ic_value) & !is.infinite(ic_value))

resumen_ic_real <- ic_real_mf %>%
  summarise(
    origen = "E. coli",
    ic_medio = mean(ic_value),
    mediana_ic = median(ic_value),
    desviacion = sd(ic_value),
    n_terminos = n()
  )

# Calcular el IC medio de las predicciones de la IA

ic_pred_mf <- predicciones_mf %>%
  mutate(ic_value = vector_ic_mf[go_id_pred]) %>%
  filter(!is.na(ic_value) & !is.infinite(ic_value))

resumen_ic_pred <- ic_pred_mf %>%
  as.data.frame() %>% 
  dplyr::group_by(model_key) %>%
  dplyr::summarise(
    ic_medio = mean(ic_value),
    mediana_ic = median(ic_value),
    desviacion = sd(ic_value),
    n_terminos = n()) %>%
  dplyr::rename(origen = model_key) 

ic_pred_mf_control <- predicciones_mf_control %>% 
  mutate(ic_value = vector_ic_mf[go_id_pred]) %>%
  filter(!is.na(ic_value) & !is.infinite(ic_value))

resumen_ic_pred_control <- ic_pred_mf_control %>%
  as.data.frame() %>% 
  dplyr::group_by(model_key) %>%
  dplyr::summarise(
    ic_medio = mean(ic_value),
    mediana_ic = median(ic_value),
    desviacion = sd(ic_value),
    n_terminos = n()) %>%
  dplyr::rename(origen = model_key) 

# Unir tablas para la comparación final
tabla_final_ic_mf <- bind_rows(resumen_ic_real,resumen_ic_pred_control, resumen_ic_pred) %>% 
  mutate( muestra= c("UniProt", "Control", "Control", "Aleatorizado", "Aleatorizado"))
print(tabla_final_ic_mf)

#########################
#    CALCULO IC BP      #
#########################

go_data_bp_ic <- godata('org.EcK12.eg.db', ont="BP", computeIC=TRUE)
vector_ic_bp <- go_data_bp_ic@IC

ic_real_bp <- EColi_largo %>%
  dplyr::select(Entry, go_id_real) %>%
  distinct() %>%
  mutate(ic_value = vector_ic_bp[go_id_real]) %>%
  filter(!is.na(ic_value) & !is.infinite(ic_value))

resumen_ic_real_bp <- ic_real_bp %>%
  summarise(
    origen = "E. coli",
    ic_medio = mean(ic_value),
    mediana_ic = median(ic_value),
    desviacion = sd(ic_value),
    n_terminos = n()
  )

ic_pred_bp <- predicciones_bp %>%
  mutate(ic_value = vector_ic_bp[go_id_pred]) %>%
  filter(!is.na(ic_value) & !is.infinite(ic_value))

ic_pred_bp_control <- predicciones_bp_control %>% 
  mutate(ic_value = vector_ic_bp[go_id_pred]) %>%
  filter(!is.na(ic_value) & !is.infinite(ic_value))

resumen_ic_pred_bp <- ic_pred_bp %>%
  dplyr::group_by(model_key) %>%
  dplyr::summarise(
    ic_medio = mean(ic_value),
    mediana_ic = median(ic_value),
    desviacion = sd(ic_value),
    n_terminos = n()
  ) %>%
  dplyr::rename(origen = model_key)

resumen_ic_pred_bp_control <- ic_pred_bp_control %>%
  dplyr::group_by(model_key) %>%
  dplyr::summarise(
    ic_medio = mean(ic_value),
    mediana_ic = median(ic_value),
    desviacion = sd(ic_value),
    n_terminos = n()
  ) %>%
  dplyr::rename(origen = model_key)

tabla_final_ic_bp <- bind_rows(resumen_ic_real_bp, resumen_ic_pred_bp_control, resumen_ic_pred_bp) %>% 
  mutate( muestra= c("UniProt", "Control", "Control", "Aleatorizado", "Aleatorizado"))
print(tabla_final_ic_bp)


#####################
#    CALCULO IC CC  #
#####################

go_data_cc_ic <- godata('org.EcK12.eg.db', ont="CC", computeIC=TRUE)
vector_ic_cc <- go_data_cc_ic@IC

ic_real_cc <- EColi_largo %>%
  dplyr::select(Entry, go_id_real) %>%
  distinct() %>%
  mutate(ic_value = vector_ic_cc[go_id_real]) %>%
  filter(!is.na(ic_value) & !is.infinite(ic_value)) # Borra los que no sean CC

resumen_ic_real_cc <- ic_real_cc %>%
  summarise(
    origen = "E. coli",
    ic_medio = mean(ic_value),
    mediana_ic = median(ic_value),
    desviacion = sd(ic_value),
    n_terminos = n()
  )

ic_pred_cc <- predicciones_cc %>%
  mutate(ic_value = vector_ic_cc[go_id_pred]) %>%
  filter(!is.na(ic_value) & !is.infinite(ic_value))

ic_pred_cc_control <-predicciones_cc_control %>% 
  mutate(ic_value = vector_ic_cc[go_id_pred]) %>%
  filter(!is.na(ic_value) & !is.infinite(ic_value))

resumen_ic_pred_cc <- ic_pred_cc %>%
  dplyr::group_by(model_key) %>%
  dplyr::summarise(
    ic_medio = mean(ic_value),
    mediana_ic = median(ic_value),
    desviacion = sd(ic_value),
    n_terminos = n()
  ) %>%
  dplyr::rename(origen = model_key)

resumen_ic_pred_cc_control <- ic_pred_cc_control %>%
  dplyr::group_by(model_key) %>%
  dplyr::summarise(
    ic_medio = mean(ic_value),
    mediana_ic = median(ic_value),
    desviacion = sd(ic_value),
    n_terminos = n()
  ) %>%
  dplyr::rename(origen = model_key)

tabla_final_ic_cc <- bind_rows(resumen_ic_real_cc, resumen_ic_pred_cc_control, resumen_ic_pred_cc) %>% 
  mutate(muestra = c("UniProt", "Control", "Control", "Aleatorizado", "Aleatorizado"))

print(tabla_final_ic_cc)

## TABLA IC COMPLETA ##
tabla_ic_completa <- bind_rows(tabla_final_ic_mf, tabla_final_ic_bp, tabla_final_ic_cc)
tabla_ic_completa <- mutate(tabla_ic_completa, GO = c("MF","MF","MF","MF","MF","BP","BP","BP","BP","BP","CC","CC","CC","CC","CC"))
tabla_ic_completa



#### TOP 10 GO terms por ontología ####
get_go_term <- function(ids) {
  Term(ids)
}
grafico_real_ontologia <- function(datos_ontologia, ontologia){
  datos_grafico_real<- datos_ontologia%>%
    count(go_id_real, sort = TRUE) %>%
    head(10) %>%
    mutate(Termino = get_go_term(go_id_real)) %>% 
    dplyr::rename(veces = n) %>% 
    mutate(Termino2 = str_wrap(Termino, width = 35))
  
  ggplot(datos_grafico_real, aes(
    x = reorder(Termino2, veces), 
    y = veces, 
    fill = veces)) +
    geom_col(color = "black", linewidth = 0.2) + 
    geom_text(aes(label = veces), hjust = -0.15, size = 7, color = "black") + 
    coord_flip() + 
    scale_y_continuous(limits = c(0, max(datos_grafico_real$veces) * 1.2)) + 
    scale_fill_viridis_c(option = "mako", direction = -1) +
    labs(
      title = paste0("Top 10 términos GO ",ontologia, " en el proteoma real de E. Coli"),
      subtitle = "Distribución natural según UniProt",
      x = "Término GO",
      y = "Número de proteínas"
    ) +
    theme_minimal(base_size = 20) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "none"
    )
}
grafico_real_ontologia(ic_real_mf, "MF")
grafico_real_ontologia(ic_real_bp, "BP")
grafico_real_ontologia(ic_real_cc, "CC")

grafico_control_ontologia <- function(datos_ontologia_control, ontologia){
  datos_grafico_control <- datos_ontologia_control %>% 
    group_by(model_key, go_id_pred) %>%
    summarise(veces_predicha = n(), .groups = "drop") %>%
    group_by(model_key) %>%
    slice_max(order_by = veces_predicha, n = 10, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(Termino = get_go_term(go_id_pred)) %>% 
    mutate(go_description_wrap = str_wrap(Termino, width = 35))
  
  grafico_control_ordenado <- ggplot(datos_grafico_control, aes(
    x = reorder_within(go_description_wrap, veces_predicha, model_key), 
    y = veces_predicha, 
    fill = veces_predicha )) +
    geom_col(color = "black", linewidth = 0.2) + 
    geom_text(aes(label = veces_predicha), 
              hjust = -0.2, size = 5, color = "black") + 
    coord_flip() + 
    facet_wrap(~ model_key, scales = "free_y") + 
    scale_x_reordered() + 
    scale_fill_viridis_c(option = "plasma", name = "Nº Proteínas") +
    labs(
      title = paste0( "Top 10 términos GO ", ontologia ," en el control E. Coli"),
      x = "Término GO asignado",
      y = "Número de proteínas") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    theme_minimal(base_size = 18) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"), 
      plot.subtitle = element_text(hjust = 0.5),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(face = "bold"),
      legend.position = "right",
      panel.grid.minor = element_blank())
  print(grafico_control_ordenado)
}
grafico_control_ontologia(ic_pred_mf_control, "MF")
grafico_control_ontologia(ic_pred_bp_control, "BP")
grafico_control_ontologia(ic_pred_cc_control, "CC")

grafico_random_ontologia <- function (datos_ontologia_random, ontologia) {
  datos_grafico <- datos_ontologia_random %>%
    group_by(model_key, go_id_pred, randomizacion) %>%
    summarise(n_por_replica = n(), .groups = "drop") %>%
    group_by(model_key, go_id_pred) %>%
    summarise(
      n_medio = mean(n_por_replica),
      apariciones_replicas = n(),
      .groups = "drop"
    ) %>%
    group_by(model_key) %>%
    slice_max(order_by = n_medio, n = 10, with_ties = FALSE) %>% 
    ungroup() %>%
    mutate(Termino = get_go_term(go_id_pred),
      Termino_wrap = str_wrap(Termino, width = 35))
    
  
  grafico_randomizados_ordenado <- ggplot(
    datos_grafico, aes(
      x = tidytext::reorder_within(Termino_wrap, n_medio, model_key),
      y = n_medio, 
      fill = as.numeric(apariciones_replicas) )) +
    geom_col(color = "black", linewidth = 0.2) +
    geom_text(aes(label = round(n_medio, 1)), 
              hjust = -0.4, size = 4.5, color = "black") +
    coord_flip() +
    facet_wrap(~ model_key, scales = "free_y") + 
    tidytext::scale_x_reordered() +
    scale_fill_viridis_c(option = "plasma", direction = 1, 
                         name = "Nº de\nréplicas",
                         breaks = 1:5,
                         limits = c(1, 5))+
    labs(
      title = paste0("Top 10 términos GO ", ontologia, " en secuencias aleatorizadas E. Coli"),
      subtitle = "Media de 5 réplicas por modelo",
      x = "Término GO asignado",
      y = "Media de proteínas") +

    scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
    
    theme_minimal(base_size = 18) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(face = "bold"),
      legend.position = "right",
      axis.text.y = element_text(lineheight = 0.8, size = 12),
      panel.grid.minor = element_blank()
    )
  
  print(grafico_randomizados_ordenado)
}
grafico_random_ontologia(ic_pred_mf, "MF")
grafico_random_ontologia(ic_pred_bp, "BP")
grafico_random_ontologia(ic_pred_cc, "CC")
