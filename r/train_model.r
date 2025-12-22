# r/train_model.r
library(argparse)
source("/app/utils.r")

# --- Argumentos ---
parser <- ArgumentParser(description="Treina um modelo SDM para uma espécie.")
parser$add_argument("--csv", type="character", required=TRUE, help="Caminho para o arquivo CSV de treinamento.")
parser$add_argument("--output_path", type="character", required=TRUE, help="Diretório onde o modelo será salvo.")
parser$add_argument("--model_name", type="character", required=TRUE, help="Nome do arquivo do modelo (sem extensão).")

args <- parser$parse_args()

# --- Lógica ---
message("Iniciando o treinamento do modelo: ", args$model_name)

tryCatch({
  generate_model(
    data_path    = dirname(args$csv),
    output_path  = args$output_path,
    model_name   = args$model_name,
    csv_filename = basename(args$csv),
    methods      = c("rf") # Usando Random Forest como exemplo
  )
  message("Modelo treinado e salvo com sucesso em: ", file.path(args$output_path, paste0(args$model_name, ".sdm")))
}, error = function(e) {
  stop("Erro durante o treinamento do modelo: ", e$message)
})
