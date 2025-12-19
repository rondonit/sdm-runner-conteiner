# Pacotes
library(raster)
library(rgbif)   # <- usar occ_data
library(dplyr)
library(tidyr)
library(sp)
library(sf)
library(usdm)
library(maps)
library(terra)
library(sdm)

# --------- util: detectar preditores a partir do CSV ---------
infer_predictors_from_csv <- function(csv_path,
                                      response_col = "sp",
                                      ignore_cols = c("decimalLongitude","decimalLatitude","species",
                                                      "year","day","month","lon","lat","x","y","id")) {
  df <- read.csv(csv_path, header = TRUE)
  # Remover colunas ignoradas (se existirem)
  ignore_cols <- intersect(ignore_cols, names(df))
  base_cols <- unique(c(ignore_cols, response_col))
  # Ficar só com colunas numéricas (boas p/ modelagem) que não estão na base_cols
  num_cols <- names(df)[sapply(df, is.numeric)]
  preds <- setdiff(num_cols, base_cols)
  if (length(preds) == 0) {
    stop("Nenhuma coluna preditora numérica encontrada no CSV após remover ignore_cols e a resposta.")
  }
  list(predictor_names = preds, response_col = response_col, ignore_cols = ignore_cols)
}

# --------- 1) GBIF: igual, só habilitando rgbif ---------
process_gbif_data <- function(species_name, lon_range = c(-85,-32), lat_range = c(-58, 15), 
                              year_range = c(2000, 2024), output_dir = getwd(), limit = 100000) {
  gbif_data <- occ_data(scientificName = species_name, hasCoordinate = TRUE,
                        decimalLongitude = paste(lon_range, collapse=","),
                        decimalLatitude  = paste(lat_range, collapse=","), limit = limit)
  species_coords <- gbif_data$data[, c("decimalLongitude","decimalLatitude","year","day","month")]
  setwd(output_dir)
  species_filtered <- species_coords %>% 
    filter(year >= year_range[1] & year <= year_range[2]) %>% 
    drop_na(decimalLongitude, decimalLatitude, year, day, month)
  unique_species <- species_filtered[!duplicated(species_filtered[, c("decimalLongitude","decimalLatitude","year","day","month")]), ]
  output_file <- paste0(gsub(" ", "_", species_name), "_GBIF_Data.csv")
  write.csv(unique_species, output_file, row.names = FALSE)
  coordinates(unique_species) <- ~ decimalLongitude + decimalLatitude
  return(unique_species)
}

# --------- 2) Treino: usa variáveis do CSV automaticamente ---------
generate_model <- function(data_path,
                           output_path,
                           model_name,
                           response_col = "sp",
                           ignore_cols = c("decimalLongitude","decimalLatitude","species","year","day","month"),
                           methods = c("rf"),
                           csv_filename = "predictor.csv",
                           n_reps = 10) {   # <- novo argumento com default 10
  setwd(data_path)
  df <- read.csv(csv_filename, sep = ",", header = TRUE)
  
  inferred <- infer_predictors_from_csv(file.path(data_path, csv_filename),
                                        response_col = response_col,
                                        ignore_cols  = ignore_cols)
  predictor_names <- inferred$predictor_names
  response_col    <- inferred$response_col
  
  if (!any(response_col %in% names(df)))
    stop(paste("Coluna de resposta não encontrada:", paste(response_col, collapse=", ")))
  
  keep_cols <- c(response_col, predictor_names)
  df <- df[, keep_cols, drop = FALSE]
  
  form <- as.formula(paste(response_col, "~ ."))
  sdmDataObj <- sdmData(form, train = df)
  
  # aqui usa o n_reps
  modelo <- sdm(form, sdmDataObj,
                methods     = methods,
                replication = c("boot"),
                n           = n_reps)
  
  setwd(output_path)
  write.sdm(modelo, model_name, overwrite = TRUE)
  return(modelo)
}


# --------- 3) Predição: reusa os nomes dos preditores do CSV ---------
predict_sdm <- function(model_path,
                        predictors_path,
                        output_file,
                        csv_path,
                        pattern = "\\.nc$",
                        resample = FALSE,
                        ref_raster_path = NULL,
                        model_name = basename(model_path),
                        response_col = "sp",
                        ignore_cols = c("decimalLongitude","decimalLatitude",
                                        "species","year","day","month")) {
  
  categorical_vars <- c("ocup_solo")
  
  inferred <- infer_predictors_from_csv(csv_path,
                                        response_col = response_col,
                                        ignore_cols  = ignore_cols)
  predictor_names <- inferred$predictor_names
  
  path   <- file.path(model_path, paste(model_name, "sdm", sep = "."))
  modelo <- sdm::read.sdm(path)
  
  nc_files <- list.files(path = predictors_path, pattern = pattern, full.names = TRUE)
  if (length(nc_files) == 0)
    stop(paste("Nenhum arquivo encontrado em", predictors_path, "com pattern", pattern))
  
  r_list <- lapply(nc_files, function(f) {
    r <- raster::raster(f)
    nm <- tools::file_path_sans_ext(basename(f))
    names(r) <- nm
    r
  })
  
  if (resample) {
    if (!is.null(ref_raster_path)) {
      if (!file.exists(ref_raster_path)) stop("ref_raster_path não existe: ", ref_raster_path)
      ref_raster <- raster::raster(ref_raster_path)
    } else {
      # escolhe raster base ignorando variáveis categóricas (ex.: ocup_solo)
      cell_areas <- sapply(seq_along(r_list), function(i) {
        x  <- r_list[[i]]
        nm <- names(x)[1]
        if (nm %in% categorical_vars) return(Inf)  # << não deixa ocup_solo ser referência
        rr <- raster::res(x)
        if (any(is.na(rr))) return(Inf)
        rr[1] * rr[2]
      })
      
      best_idx <- which.min(cell_areas)
      if (!is.finite(cell_areas[best_idx])) {
        stop("Não consegui escolher raster base (todos foram ignorados/invalidos).")
      }
      
      ref_raster <- r_list[[best_idx]]
      message(sprintf("Resample automático: usando '%s' como base.", basename(nc_files[best_idx])))
      message(sprintf("Resolução escolhida: %f x %f", raster::res(ref_raster)[1], raster::res(ref_raster)[2]))
    }
    
    r_list <- lapply(r_list, function(x) {
      nm <- names(x)[1]
      
      # Categóricas (ex.: ocup_solo): evitar bilinear.
      # Quando estiver fazendo downsample (fino -> grosso), primeiro agrega por "modal"
      # para representar a classe majoritária dentro da célula de saída.
      if (nm %in% categorical_vars) {
        rx <- raster::res(x)
        rr <- raster::res(ref_raster)
        
        if (all(!is.na(rx)) && all(!is.na(rr)) && (rr[1] > rx[1] || rr[2] > rx[2])) {
          fx <- max(1L, as.integer(round(rr[1] / rx[1])))
          fy <- max(1L, as.integer(round(rr[2] / rx[2])))
          
          if (fx > 1L || fy > 1L) {
            message(sprintf(
              "Categórico '%s': downsample -> aggregate(modal) fact=(%d,%d) antes do ngb.",
              nm, fx, fy
            ))
            x <- raster::aggregate(
              x,
              fact = c(fx, fy),
              fun = raster::modal,
              na.rm = TRUE,
              expand = TRUE,
              filename = raster::rasterTmpFile(),
              overwrite = TRUE
            )
          }
        }
        
        return(raster::resample(
          x, ref_raster,
          method = "ngb",
          filename = raster::rasterTmpFile(),
          overwrite = TRUE
        ))
      }
      
      # Contínuas
      raster::resample(
        x, ref_raster,
        method = "bilinear",
        filename = raster::rasterTmpFile(),
        overwrite = TRUE
      )
    })
  }
  
  r_stack <- raster::stack(r_list)
  current_names <- names(r_stack)
  
  matched_idx <- integer(0)
  for (v in predictor_names) {
    hit <- which(tolower(current_names) == tolower(v))
    if (length(hit) == 0) hit <- grep(v, current_names, ignore.case = TRUE)
    if (length(hit) == 0) stop("Não encontrei camada para: ", v)
    matched_idx <- c(matched_idx, hit[1])
  }
  
  r_stack <- r_stack[[matched_idx]]
  names(r_stack) <- predictor_names
  
  p <- sdm::predict(modelo, r_stack, overwrite = TRUE)
  raster::writeRaster(p, filename = output_file, overwrite = TRUE, datatype = "FLT4S")
}


# --------- 4) VIF genérico: usa TODAS variáveis numéricas do CSV ---------
generate_vif <- function(csv_path,
                         response_col = "sp",
                         ignore_cols = c("decimalLongitude","decimalLatitude","species","year","day","month"),
                         res = 0.5,
                         vif_threshold = 7) {
  df <- read.csv(csv_path, header = TRUE)
  
  # Coordenadas obrigatórias para rasterizar
  coord_candidates <- c("decimalLongitude","decimalLatitude","lon","lat","x","y")
  lon_col <- coord_candidates[coord_candidates %in% names(df)][1]
  lat_col <- coord_candidates[coord_candidates %in% names(df)][2]
  if (is.na(lon_col) || is.na(lat_col)) {
    stop("Não encontrei colunas de coordenadas (ex.: decimalLongitude/decimalLatitude ou lon/lat) no CSV.")
  }
  
  # Detectar variáveis numéricas preditoras (sem nomes fixos)
  inferred <- infer_predictors_from_csv(csv_path, response_col = response_col, ignore_cols = ignore_cols)
  predictor_names <- inferred$predictor_names
  
  # Converter para Spatial para poder rasterizar
  coordinates(df) <- stats::as.formula(paste0("~", lon_col, "+", lat_col))
  ext <- extent(df)
  r_template <- raster(ext, res = res)
  
  # Rasterizar cada variável preditora e empilhar
  r_list <- list()
  for (v in predictor_names) {
    vals <- rasterize(df, r_template, field = v, fun = mean)
    r_list[[v]] <- vals
  }
  env_stack <- stack(r_list)
  names(env_stack) <- predictor_names
  
  # usdm trabalha com Raster*; opcionalmente converter para SpatRaster e de volta
  # data <- rast(env_stack) # (terra) – não é necessário para usdm
  
  message("VIFSTEP")
  print(vifstep(env_stack, th = vif_threshold))
  message("VIFCOR")
  print(vifcor(env_stack, th = vif_threshold))
}

# --------- 5) Pipeline com vários métodos (ex: BRT), usando CSV ---------
generate_5_methods <- function(base_path, species_name,
                               response_col = "sp",
                               ignore_cols = c("decimalLongitude","decimalLatitude","species","year","day","month"),
                               csv_filename = "predictor.csv") {
  methods <- c("brt")  # adicione outros: "rf","svm","glm","maxent", etc. conforme instalado
  predictors_daily <- file.path(base_path, "daily")
  predictors_annual <- file.path(base_path, "anual")
  
  for (method in methods) {
    model_name <- paste0(species_name, "_", method)
    message(model_name)
    
    # Treinar
    model <- generate_model(
      data_path   = base_path,
      output_path = base_path,
      model_name  = model_name,
      response_col = response_col,
      ignore_cols  = ignore_cols,
      methods      = c(method),
      csv_filename = csv_filename
    )
    
    # Predição diária (ex.: NetCDF)
    predict_sdm(
      model_path      = base_path,
      predictors_path = predictors_daily,
      output_file     = file.path(base_path, paste0(model_name, "_diario.tif")),
      csv_path        = file.path(base_path, csv_filename),
      pattern         = "*nc$",
      resample        = TRUE,
      model_name      = model_name,
      response_col    = response_col,
      ignore_cols     = ignore_cols
    )
    
    # Predição anual (ex.: GeoTIFF)
    predict_sdm(
      model_path      = base_path,
      predictors_path = predictors_annual,
      output_file     = file.path(base_path, paste0(model_name, "_anual.tif")),
      csv_path        = file.path(base_path, csv_filename),
      pattern         = "*tif$",
      resample        = TRUE,
      model_name      = model_name,
      response_col    = response_col,
      ignore_cols     = ignore_cols
    )
  }
}

# ------------------------------------------------------------
# 6) Relatório completo GUI: CSV + HTML (inclui todos bootstraps)
# ------------------------------------------------------------

generate_gui_report <- function(model = NULL,
                                model_path = NULL,
                                model_name = NULL,
                                output_csv = NULL,
                                output_html = NULL,
                                title = "Avaliação dos modelos SDM") {
  # Precisamos de: OU um objeto de modelo, OU caminho+nome
  if (is.null(model) && (is.null(model_path) || is.null(model_name))) {
    stop("Informe 'model' OU ('model_path' e 'model_name').")
  }
  
  # Se não veio output_csv/output_html, cria defaults com base no caminho/nome
  if (is.null(output_csv) && !is.null(model_path) && !is.null(model_name)) {
    output_csv <- file.path(model_path, paste0(model_name, "_eval.csv"))
  }
  if (is.null(output_html) && !is.null(model_path) && !is.null(model_name)) {
    output_html <- file.path(model_path, paste0(model_name, "_eval.html"))
  }
  
  # 1) Obter o modelo (em memória ou lendo o .sdm)
  if (!is.null(model)) {
    modelo <- model
  } else {
    sdm_file <- file.path(model_path, paste0(model_name, ".sdm"))
    if (!file.exists(sdm_file)) {
      stop("Arquivo .sdm não encontrado: ", sdm_file)
    }
    modelo <- read.sdm(sdm_file)
  }
  
  # 2) Puxar todas as avaliações (cada linha = 1 replicação/método/tipo de dado)
  ev <- getEvaluation(modelo)
  
  # 3) Salvar CSV com TODAS as linhas (todos os 10 bootstraps)
  if (!is.null(output_csv)) {
    write.csv(ev, output_csv, row.names = FALSE)
  }
  
  # 4) Gerar HTML (resumo + tabela completa)
  if (!is.null(output_html)) {
    
    round_df <- function(df, digits = 3) {
      num_cols <- sapply(df, is.numeric)
      if (any(num_cols)) {
        df[, num_cols] <- lapply(df[, num_cols, drop = FALSE], round, digits = digits)
      }
      df
    }
    
    df_to_html_table <- function(df) {
      df <- round_df(df)
      header <- paste0(
        "<tr>",
        paste(sprintf("<th>%s</th>", names(df)), collapse = ""),
        "</tr>"
      )
      rows <- apply(df, 1, function(r) {
        paste0(
          "<tr>",
          paste(sprintf("<td>%s</td>", r), collapse = ""),
          "</tr>"
        )
      })
      paste0("<table>", header, paste(rows, collapse = "\n"), "</table>")
    }
    
    # Resumo por método/tipo de dado (médias das métricas)
    meta_cols   <- intersect(c("model", "method", "algorithm", "replication", "dataType"), names(ev))
    metric_cols <- setdiff(names(ev), meta_cols)
    group_cols  <- intersect(c("method", "dataType"), names(ev))
    summ <- NULL
    if (length(group_cols) > 0 && length(metric_cols) > 0) {
      formula_str <- paste(
        paste(metric_cols, collapse = " + "),
        "~",
        paste(group_cols, collapse = " + ")
      )
      sum_formula <- as.formula(formula_str)
      summ <- aggregate(sum_formula, data = ev, FUN = mean, na.rm = TRUE)
    }
    
    html <- c(
      "<!DOCTYPE html>",
      "<html>",
      "<head>",
      "<meta charset='UTF-8'>",
      sprintf("<title>%s</title>", title),
      "<style>
      body { font-family: Arial, sans-serif; margin: 20px; }
      h1, h2 { font-family: Arial, sans-serif; }
      table { border-collapse: collapse; margin-bottom: 20px; }
      th, td { border: 1px solid #ccc; padding: 4px 8px; font-size: 12px; }
      th { background-color: #f5f5f5; }
    </style>",
      "</head>",
      "<body>",
      sprintf("<h1>%s</h1>", title)
    )
    
    if (!is.null(model_name)) {
      html <- c(html, sprintf("<h2>Modelo: %s</h2>", model_name))
    }
    
    if (!is.null(summ)) {
      html <- c(
        html,
        "<h2>Resumo por método / tipo de dado (médias)</h2>",
        df_to_html_table(summ)
      )
    }
    
    # Aqui vai a tabela COMPLETA: todos os bootstraps
    html <- c(
      html,
      "<h2>Todas as avaliações (cada linha = 1 bootstrap / replicação)</h2>",
      df_to_html_table(ev),
      "</body>",
      "</html>"
    )
    
    writeLines(html, output_html)
  }
  
  # Devolve o data.frame completo (caso tu queira usar em R)
  return(ev)
}

# ------------------------------------------------------------
# 7) Checar range de variáveis por área/tempo (CSV) e TIFF
# ------------------------------------------------------------
check_sdm_ranges <- function(csv_path,
                             raster_path = NULL,
                             response_col = "sp",
                             ignore_cols = c(
                               "decimalLongitude","decimalLatitude","species",
                               "year","day","month","lon","lat","x","y","id"
                             ),
                             area_col = NULL,   # ex: "regiao"
                             time_col = NULL    # ex: "year"
) {
  if (!file.exists(csv_path)) {
    stop("CSV não encontrado em: ", csv_path)
  }
  
  # --------- 1) Lê o CSV e descobre preditores ---------
  df <- read.csv(csv_path, header = TRUE)
  
  inferred <- infer_predictors_from_csv(
    csv_path     = csv_path,
    response_col = response_col,
    ignore_cols  = ignore_cols
  )
  predictor_names <- inferred$predictor_names
  
  # Garante que as colunas existem
  area_col <- if (!is.null(area_col) && area_col %in% names(df)) area_col else NULL
  time_col <- if (!is.null(time_col) && time_col %in% names(df)) time_col else NULL
  
  group_vars <- c()
  if (!is.null(area_col)) group_vars <- c(group_vars, area_col)
  if (!is.null(time_col)) group_vars <- c(group_vars, time_col)
  
  if (length(group_vars) == 0) {
    # Se não tiver nenhuma coluna de área/tempo, cria um grupo único
    df$..grupo.. <- "tudo"
    group_vars <- "..grupo.."
  }
  
  # Fica só com colunas de interesse
  sel_cols <- c(predictor_names, group_vars)
  sel_cols <- intersect(sel_cols, names(df))
  df_sel   <- df[, sel_cols, drop = FALSE]
  
  # --------- 2) Tabela longa e estatísticas por grupo ---------
  long_df <- tidyr::pivot_longer(
    data      = df_sel,
    cols      = dplyr::all_of(predictor_names),
    names_to  = "variable",
    values_to = "value"
  )
  
  predictor_ranges <- long_df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars)), .data$variable) %>%
    dplyr::summarise(
      n      = sum(!is.na(.data$value)),
      min    = ifelse(n > 0, min(.data$value, na.rm = TRUE), NA_real_),
      q25    = ifelse(n > 0, stats::quantile(.data$value, 0.25, na.rm = TRUE), NA_real_),
      median = ifelse(n > 0, stats::median(.data$value, na.rm = TRUE), NA_real_),
      q75    = ifelse(n > 0, stats::quantile(.data$value, 0.75, na.rm = TRUE), NA_real_),
      max    = ifelse(n > 0, max(.data$value, na.rm = TRUE), NA_real_),
      mean   = ifelse(n > 0, mean(.data$value, na.rm = TRUE), NA_real_),
      sd     = ifelse(n > 1, stats::sd(.data$value, na.rm = TRUE), NA_real_),
      .groups = "drop"
    )
  
  # --------- 3) Estatísticas do raster de predição (opcional) ---------
  raster_ranges <- NULL
  if (!is.null(raster_path)) {
    if (!file.exists(raster_path)) {
      warning("Raster de predição não encontrado em: ", raster_path)
    } else {
      r_stack <- raster::stack(raster_path)
      n_layers <- raster::nlayers(r_stack)
      layer_names <- names(r_stack)
      
      raster_ranges <- data.frame(
        layer = layer_names,
        min   = NA_real_,
        max   = NA_real_,
        mean  = NA_real_,
        sd    = NA_real_,
        stringsAsFactors = FALSE
      )
      
      for (i in seq_len(n_layers)) {
        lay <- r_stack[[i]]
        # AQUI a correção: usar cellStats em tudo
        raster_ranges$min[i]  <- raster::cellStats(lay, stat = "min",  na.rm = TRUE)
        raster_ranges$max[i]  <- raster::cellStats(lay, stat = "max",  na.rm = TRUE)
        raster_ranges$mean[i] <- raster::cellStats(lay, stat = "mean", na.rm = TRUE)
        raster_ranges$sd[i]   <- raster::cellStats(lay, stat = "sd",   na.rm = TRUE)
      }
    }
  }
  
  # --------- 4) Retorno ---------
  out <- list(
    predictor_ranges = predictor_ranges,
    raster_ranges    = raster_ranges
  )
  return(out)
}



# -------------------- EXEMPLOS DE USO --------------------
# 1) Treinar um modelo RF lendo variáveis do predictor.csv automaticamente:
# model <- generate_model(
#   data_path   = "/caminho/para/pasta",
#   output_path = "/caminho/para/pasta",
#   model_name  = "minha_especie_rf",
#   response_col = "sp",
#   ignore_cols  = c("decimalLongitude","decimalLatitude","species","year","day","month"),
#   methods      = c("rf"),
#   csv_filename = "predictor.csv"
# )

# 2) Predizer usando os nomes de variáveis vindos do CSV (.nc):
# predict_sdm(
#   model_path      = "/caminho/para/pasta",
#   predictors_path = "/caminho/para/pasta/daily",
#   output_file     = "/caminho/para/pasta/minha_especie_rf_diario.tif",
#   csv_path        = "/caminho/para/pasta/predictor.csv",
#   pattern         = "*nc$",
#   resample        = TRUE,
#   model_name      = "minha_especie_rf",
#   response_col    = "sp"
# )

# 3) VIF automático a partir do CSV (sem nomes fixos):
# generate_vif(
#   csv_path      = "/caminho/para/pasta/predictor.csv",
#   response_col  = "sp",
#   ignore_cols   = c("decimalLongitude","decimalLatitude","species","year","day","month"),
#   res           = 0.5,
#   vif_threshold = 7
# )

# 4) Pipeline (ex.: BRT) já usando CSV para descobrir variáveis:
# generate_5_methods(
#   base_path    = "/caminho/para/pasta",
#   species_name = "Thunnus_obesus",
#   response_col = "sp"
# )

##########################
########## OLD ###########
##########################
# library(raster)
# #library(rgbif)
# library(dplyr)
# library(tidyr)
# library(sp)
# library(sf)
# library(usdm)
# library(maps)
# library(tidyr)
# library(terra)
# library(sdm)

# process_gbif_data <- function(species_name, lon_range = c(-85,-32), lat_range = c(-58, 15), 
#                               year_range = c(2000, 2024), output_dir = getwd(),limit=100000) {
#   gbif_data <- occ_data(scientificName = species_name, hasCoordinate = TRUE,
#                         decimalLongitude = paste(lon_range, collapse=","), 
#                         decimalLatitude = paste(lat_range, collapse=","),limit=limit)
#   species_coords <- gbif_data$data[, c("decimalLongitude", "decimalLatitude", "year","day","month")]
#   setwd(output_dir)
#   species_filtered <- species_coords %>% 
#     filter(year >= year_range[1] & year <= year_range[2]) %>% 
#     drop_na(decimalLongitude, decimalLatitude,year,day,month)
#   unique_species <- species_filtered[!duplicated(species_filtered[, c("decimalLongitude", "decimalLatitude","year","day","month")]), ]
#   output_file <- paste0(gsub(" ", "_", species_name), "_GBIF_Data.csv")
#   write.csv(unique_species, output_file, row.names = FALSE)
#   coordinates(unique_species) <- ~decimalLongitude + decimalLatitude
#   return(unique_species)
# }


# generate_model <- function(data_path, output_path, model_name,ignore_cols,methods=c('rf')) {
#   setwd(data_path)
#   df <- read.csv("predictor.csv", sep = ",", header = TRUE)
#   head(df)
#   cols.dont.want <- ignore_cols
#   df <- df[, ! names(df) %in% cols.dont.want, drop = FALSE]
#   sdmData <- sdmData(sp~ ., train = df)
#   modelo <- sdm(sp~ ., sdmData, methods = methods, replication = c('boot'), n = 7)
#   setwd(output_path)
#   write.sdm(modelo, model_name, overwrite = TRUE)
#   return (modelo)
# }

# predict_sdm <- function(model_path, predictors_path, output_file,variables,pattern="*nc$",resample=FALSE,model_name=basename(model_path)) {
#   setwd(model_path)

#   path<-file.path(model_path, paste(model_name,"sdm",sep="."))
#   print(path)
#   modelo <- read.sdm(path)
#   #print(modelo)
#   nc_files <- list.files(path = predictors_path, pattern = pattern, full.names = TRUE)
#   r <- lapply(nc_files, raster)
#   if (resample) {
#     ref_raster <- r[[1]]
#     r_resampled <- lapply(r, function(x) resample(x, ref_raster, method = "bilinear"))
#     bio_presente_cop <- stack(r_resampled)
#   } else {
#     bio_presente_cop <- stack(r)
#   }
#   names(bio_presente_cop) <- variables
#   print(bio_presente_cop)
#   p <- predict(modelo, bio_presente_cop,overwrite = TRUE)
#   writeRaster(p, filename = output_file, overwrite = TRUE, datatype = 'FLT4S')
# }



# generate_vif <- function(csv_path){
#   data <- read.csv(csv_path)
#   head(data)
#   coordinates(data) <- ~decimalLongitude + decimalLatitude
#   ext <- extent(data)
#   raster_template <- raster(ext, res = 0.5)
#   thetao_raster <- rasterize(data, raster_template, field = "thetao", fun = mean)
#   so_raster <- rasterize(data, raster_template, field = "so", fun = mean)
#   chl_raster <- rasterize(data, raster_template, field = "CHL", fun = mean)
#   depth_raster <- rasterize(data, raster_template, field = "depth", fun = mean)
#   env_stack <- stack(thetao_raster, so_raster, chl_raster,depth_raster)
#   names(env_stack) <- c("thetao", "so", "CHL","depth")
#   data<-rast(env_stack)
#   print("VIFSTEP")
#   print(vifstep(data,th=7))
#   print("VIFCORE")
#   print(vifcor(data,th=7))
# }


# generate_5_methods<-function(base_path,species_name){
#   methods <- c(
#     'rf', 
#     'glm', 
#     'gam',
#     'brt',
#     'maxent'
#   )
#   predictors_daily <- file.path(base_path, "predictors")
#   predictors_annual <- file.path(base_path, "anual")
#   for (method in methods) {
#     model_name <- paste0(species_name,"_", method)
#     print(model_name)

#     model <- generate_model(
#       base_path,
#       base_path,
#       model_name,
#       c('decimalLongitude','decimalLatitude','species','year','day','month'),
#       c(method)
#     )


#     predict_sdm(
#       model_path = base_path,
#       predictors_path = predictors_daily,
#       output_file = file.path(base_path, paste0(model_name, "_diario.tif")),

#       variables = c('CHL','depth','so','thetao'),
#       pattern = '*nc$',
#       resample = TRUE,
#       model_name = model_name
#     )

#     predict_sdm(
#       model_path = base_path,
#       predictors_path = predictors_annual,
#       output_file = file.path(base_path, paste0(model_name, "_anual.tif")),
#       variables = c('CHL','depth','so','thetao'),
#       pattern = '*tif$',
#       resample = TRUE,
#       model_name = model_name
#     )
#   }

# }
