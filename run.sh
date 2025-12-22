#!/bin/bash
# Configuração global:
set -e
S3_BASE_BUCKET="s3://ste-siapesq"
# Define o diretório do projeto dentro da pasta home do usuário atual
PROJECT_DIR="/home/ubuntu/sdm_project"

# --- Preparação Inicial ---
echo "Baixando CSV de treinamento para 'guariba'."
mkdir -p "$PROJECT_DIR/data/guariba"
aws s3 cp "$S3_BASE_BUCKET/species-table/guariba_train_data.csv" "$PROJECT_DIR/data/guariba/guariba_train_data.csv"
echo ""

# --- Treinamento ---
echo "--- Iniciando TREINAMENTO para: guariba ---"
mkdir -p "$PROJECT_DIR/output/guariba_model"

echo "[1/1] Treinando o modelo para guariba..."
# Assumindo que seus scripts R estão no diretório do projeto
Rscript "$PROJECT_DIR/r/train_model.r" \
  --csv "$PROJECT_DIR/data/guariba/guariba_train_data.csv" \
  --output_path "$PROJECT_DIR/output/guariba_model" \
  --model_name "guariba_rf_model"

echo "--- Treinamento de 'guariba' CONCLUÍDO! Modelo salvo em $PROJECT_DIR/output/guariba_model ---"
echo ""

# --- Predição 2050 ---
echo "--- Iniciando PREDIÇÃO para: guariba (Cenário 2050) ---"
mkdir -p "$PROJECT_DIR/predictors/2050"
mkdir -p "$PROJECT_DIR/output/guariba_2050"

echo "[1/4] Baixando preditores NetCDF para o cenário 2050..."
aws s3 cp "$S3_BASE_BUCKET/variables/future/2050/e.nc" "$PROJECT_DIR/predictors/2050/e.nc"
aws s3 cp "$S3_BASE_BUCKET/variables/future/2050/tp.nc" "$PROJECT_DIR/predictors/2050/tp.nc"
aws s3 cp "$S3_BASE_BUCKET/variables/future/2050/t2m.nc" "$PROJECT_DIR/predictors/2050/t2m.nc"
aws s3 cp "$S3_BASE_BUCKET/variables/ocup_solo.nc" "$PROJECT_DIR/predictors/2050/ocup_solo.nc"
aws s3 cp "$S3_BASE_BUCKET/variables/sdtb_sul_expandido_SUL.nc" "$PROJECT_DIR/predictors/2050/sdtb_sul_expandido_SUL.nc"
aws s3 cp "$S3_BASE_BUCKET/variables/altitude.nc" "$PROJECT_DIR/predictors/2050/altitude.nc"

echo "[2/4] Gerando predição para 2050..."
Rscript "$PROJECT_DIR/r/predict_sdm.r" \
  --model_path "$PROJECT_DIR/output/guariba_model" \
  --model_name "guariba_rf_model" \
  --predictors_path "$PROJECT_DIR/predictors/2050" \
  --csv "$PROJECT_DIR/data/guariba/guariba_train_data.csv" \
  --output_file "$PROJECT_DIR/output/guariba_2050/guariba_2050_prediction.tif"

echo "[3/4] Enviando resultado de 2050 para o S3..."
aws s3 cp "$PROJECT_DIR/output/guariba_2050/guariba_2050_prediction.tif" "$S3_BASE_BUCKET/predictions/guariba/guariba_2050_prediction.tif"

echo "[4/4] Limpando dados da predição 2050..."
rm -rf "$PROJECT_DIR/predictors/2050"
rm -rf "$PROJECT_DIR/output/guariba_2050"
echo "--- Predição 2050 para 'guariba' CONCLUÍDA! ---"
echo ""

# --- Limpeza Final ---
echo "--- Limpando arquivos restantes salvos localmente... ---"
rm -rf "$PROJECT_DIR/output/guariba_model"
rm -rf "$PROJECT_DIR/data/guariba"
echo ""
echo "Todos os processos foram concluídos com sucesso!"