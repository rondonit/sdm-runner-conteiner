#!/bin/bash
# Configuração global:
set -e
S3_BASE_BUCKET="s3://ste-siapesq"

echo "Baixando CSV de treinamento para 'guariba'."
mkdir -p "/app/data/guariba"
aws s3 cp "$S3_BASE_BUCKET/species-table/guariba_train_data.csv" "/app/data/guariba/guariba_train_data.csv"
echo ""

echo "--- Iniciando TREINAMENTO para: guariba ---"

# --- 1. Preparação dos diretórios locais para o TREINAMENTO ---
mkdir -p "/app/output/guariba_model"

# --- 2. Execução do Treinamento ---
echo "[1/1] Treinando o modelo para guariba..."
Rscript /app/r/train_model.r \
  --csv "/app/data/guariba/guariba_train_data.csv" \
  --output_path "/app/output/guariba_model" \
  --model_name "guariba_rf_model"

echo "--- Treinamento de 'guariba' CONCLUÍDO! Modelo salvo em /app/output/guariba_model ---"
echo ""

# ==============================================================================
echo "--- Iniciando PREDIÇÃO para: guariba (Cenário 2050) ---"

# cria diretorios locais para a predição
mkdir -p "/app/predictors/2050"
mkdir -p "/app/output/guariba_2050"

echo "[1/4] Baixando preditores NetCDF para o cenário 2050..."
# Arquivos específicos do cenário
aws s3 cp "$S3_BASE_BUCKET/variables/future/2050/e.nc" "/app/predictors/2050/e.nc"
aws s3 cp "$S3_BASE_BUCKET/variables/future/2050/tp.nc" "/app/predictors/2050/tp.nc"
aws s3 cp "$S3_BASE_BUCKET/variables/future/2050/t2m.nc" "/app/predictors/2050/t2m.nc"
# Arquivos comuns a todos os cenarios
aws s3 cp "$S3_BASE_BUCKET/variables/ocup_solo.nc" "/app/predictors/2050/ocup_solo.nc"
aws s3 cp "$S3_BASE_BUCKET/variables/sdtb_sul_expandido_SUL.nc" "/app/predictors/2050/sdtb_sul_expandido_SUL.nc"
aws s3 cp "$S3_BASE_BUCKET/variables/altitude.nc" "/app/predictors/2050/altitude.nc"

echo "[2/4] Gerando predição para 2050..."
Rscript /app/r/predict_sdm.r \
  --model_path "/app/output/guariba_model" \
  --model_name "guariba_rf_model" \
  --predictors_path "/app/predictors/2050" \
  --csv "/app/data/guariba/guariba_train_data.csv" \
  --output_file "/app/output/guariba_2050/guariba_2050_prediction.tif"

echo "[3/4] Enviando resultado de 2050 para o S3..."
aws s3 cp "/app/output/guariba_2050/guariba_2050_prediction.tif" "$S3_BASE_BUCKET/predictions/guariba/guariba_2050_prediction.tif"

echo "[4/4] Limpando dados da predição 2050..."
rm -rf "/app/predictors/2050"
rm -rf "/app/output/guariba_2050"
echo "--- Predição 2050 para 'guariba' CONCLUÍDA! ---"
echo ""

# ==============================================================================
echo "--- Limpando arquivos restantes salvos localmente... ---"
rm -rf "/app/output/guariba_model"
rm -rf "/app/data/guariba"
echo ""
echo "Todos os processos foram concluídos com sucesso!"