# r/predict_sdm.r
library(argparse)
source("/app/utils.r")

# --- Argumentos ---
parser <- ArgumentParser(description="Gera predições SDM a partir de um modelo treinado.")
parser$add_argument("--model_path", type="character", required=TRUE, help="Diretório onde o modelo .sdm está salvo.")
parser$add_argument("--model_name", type="character", required=TRUE, help="Nome do arquivo do modelo (sem extensão).")
parser$add_argument("--predictors_path", type="character", required=TRUE, help="Diretório com os rasters preditores (arquivos .nc).")
parser$add_argument("--csv", type="character", required=TRUE, help="Caminho para o arquivo CSV original de treinamento (usado para obter a lista de preditores).")
parser$add_argument("--output_file", type="character", required=TRUE, help="Caminho completo para o arquivo de saída .tif.")

args <- parser$parse_args()

# --- Lógica ---
message("Iniciando a geração de predição para o modelo: ", args$model_name)
message("Usando preditores de: ", args$predictors_path)
message("O resultado será salvo em: ", args$output_file)

tryCatch({
  predict_sdm(
    model_path      = args$model_path,
    model_name      = args$model_name,
    predictors_path = args$predictors_path,
    csv_path        = args$csv,
    output_file     = args$output_file,
    pattern         = "\\.nc$",
    resample        = TRUE
  )
  message("Predição gerada com sucesso: ", args$output_file)
}, error = function(e) {
  stop("Erro durante a predição: ", e$message)
})
