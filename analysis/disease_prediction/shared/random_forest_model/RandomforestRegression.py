#!/usr/bin/python

# -*- coding: utf-8 -*-
"""
Este script entrena un modelo de Random Forest para una tarea de REGRESIÓN.
Utiliza datos de características (por ejemplo, OTUs) para predecir un valor numérico continuo.
Los pasos incluyen:
1. Cargar y preprocesar los datos de entrenamiento y prueba.
2. Opcionalmente, realizar un ajuste de hiperparámetros.
3. Entrenar un modelo RandomForestRegressor.
4. Evaluar el modelo en el conjunto de prueba utilizando métricas de regresión (R^2, MSE, MAE).
5. Guardar las importancias de las características, las predicciones y un gráfico de evaluación.
"""

import sys
import biom
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
import matplotlib.pyplot as plt

# --- Argumentos de la línea de comandos ---
# Ejemplo de uso:
# python tu_script.py train.biom test.biom metadata.tsv embedding.txt target_variable plot.png predictions.csv feature_weights.csv
train_file = sys.argv[1]
test_file = sys.argv[2]
metadata = sys.argv[3]
embedding_file = sys.argv[4]
group = sys.argv[5] # ¡ASEGÚRATE de que esta columna en los metadatos sea numérica/continua!
plot_file = sys.argv[6]
PREDICTED_VALUES_FILE = sys.argv[7] # Cambiado de PROBS a VALUES
feature_weights_path = sys.argv[8]


def computeRegressionStats(m,
                           data,
                           y_true,
                           sample_ids=None,
                           predictions_output_file=None,
                           plot=False,
                           plot_file=False,
                           graph_title=None):
    """
    Calcula métricas de regresión y opcionalmente guarda predicciones y un gráfico.
    """
    # Para regresión, predecimos el valor continuo directamente
    y_pred = m.predict(data)

    if sample_ids is not None and predictions_output_file is not None:
        # Crear un DataFrame para guardar los resultados
        predictions_df = pd.DataFrame({
            'SampleID': sample_ids,
            'TrueValue': y_true,
            'PredictedValue': y_pred
        })
        predictions_df.to_csv(predictions_output_file, index=False)
        print(f"Las predicciones de las muestras se han guardado en: {predictions_output_file}")

    # Calcular métricas de regresión estándar
    mse = mean_squared_error(y_true, y_pred)
    mae = mean_absolute_error(y_true, y_pred)
    r2 = r2_score(y_true, y_pred)
    rmse = np.sqrt(mse)  # El Error Cuadrático Medio Raíz es a menudo más interpretable

    if plot:
        plt.figure(figsize=(8, 8))
        plt.scatter(y_true, y_pred, alpha=0.6, edgecolors='k', s=50)
        # Añadir una línea para predicciones perfectas (y=x)
        lims = [
            np.min([np.min(y_true), np.min(y_pred)]) * 0.95,
            np.max([np.max(y_true), np.max(y_pred)]) * 1.05
        ]
        plt.plot(lims, lims, 'r-', alpha=0.75, zorder=0, linewidth=2)
        plt.xlim(lims)
        plt.ylim(lims)

        plt.title(graph_title or 'Valores Reales vs. Predichos', fontsize=16)
        plt.xlabel('Valores Reales', fontsize=12)
        plt.ylabel('Valores Predichos', fontsize=12)
        plt.grid(True, linestyle='--', alpha=0.6)
        # Añadir métricas al gráfico
        plt.text(lims[0], lims[1], f'$R^2 = {r2:.3f}$\nMAE = {mae:.3f}\nRMSE = {rmse:.3f}',
                 verticalalignment='top', horizontalalignment='left',
                 bbox={'boxstyle': 'round', 'facecolor': 'wheat', 'alpha': 0.5})
        plt.tight_layout()
        plt.savefig(plot_file)
        plt.close()

    return r2, mae, mse, rmse, y_pred


def getFeatureImportance(m, data):
    """
    Extrae y ordena las importancias de las características de un modelo entrenado.
    """
    feat_imp = m.feature_importances_
    feat_imp_labeled = zip(data.columns.values, feat_imp)
    feat_imp_sort = sorted(feat_imp_labeled, key=lambda t: t[1], reverse=True)
    return feat_imp_sort


def predictIBD_regression(X_train,
                           y_train,
                           X_test,
                           y_test,
                           graph_title="",
                           max_depth=12,
                           n_estimators=140,
                           plot=False,
                           plot_file_path=None,
                           feat_imp=False):
    """
    Entrena y evalúa un modelo RandomForestRegressor.
    """
    m = RandomForestRegressor(max_depth=max_depth,
                              random_state=0,
                              n_estimators=n_estimators,
                              n_jobs=24)
    m.fit(X_train, y_train)

    r2, mae, mse, rmse, y_pred = computeRegressionStats(
        m,
        data=X_test,
        y_true=y_test,
        plot=plot,
        plot_file=plot_file_path,
        graph_title=graph_title
    )
    
    feat_imp_sort = []
    if feat_imp:
        feat_imp_sort = getFeatureImportance(m, data=X_train)

    return m, r2, mae, mse, feat_imp_sort


def crossValPrediction_regression(otu_use,
                                  y,
                                  max_depth=10,
                                  n_estimators=65,
                                  plot=False,
                                  folds=5):
    """
    Realiza validación cruzada para el modelo de regresión.
    """
    kf = KFold(n_splits=folds, shuffle=True, random_state=42)
    
    R2_scores, MAE_scores, MSE_scores = [], [], []
    feat_imp_crossVal = []
    models = []
    
    for train_index, val_index in kf.split(otu_use):
        otu_train, otu_val = otu_use.iloc[train_index, :], otu_use.iloc[val_index, :]
        y_train, y_val = np.array(y)[train_index], np.array(y)[val_index]

        model, r2, mae, mse, feat_imp = predictIBD_regression(
            otu_train, y_train,
            otu_val, y_val,
            max_depth=max_depth,
            n_estimators=n_estimators,
            feat_imp=True,
            plot=plot # Plot individual fold performance if needed
        )
        models.append(model)
        R2_scores.append(r2)
        MAE_scores.append(mae)
        MSE_scores.append(mse)
        feat_imp_crossVal.append(feat_imp)
        
    return models, R2_scores, MAE_scores, MSE_scores, feat_imp_crossVal


def trainHyperParameters_regression(X_train, y_train):
    """
    Realiza una búsqueda en cuadrícula para encontrar los mejores hiperparámetros para la regresión.
    """
    depths = [3, 5, 7, 10, 15]
    n_estimators_list = [50, 100, 150, 200]

    df_results = []
    
    for depth in depths:
        for trees in n_estimators_list:
            print(f"Probando: Profundidad={depth}, Estimadores={trees}")
            _, r2_scores, _, _, _ = crossValPrediction_regression(
                otu_use=X_train,
                y=y_train,
                max_depth=depth,
                n_estimators=trees,
                plot=False,
                folds=5)
            
            # Optimizar para el R^2 más alto
            df_results.append({
                'R2_mean': np.mean(r2_scores),
                'R2_std': np.std(r2_scores),
                'depth': depth,
                'n_trees': trees
            })
            
    return pd.DataFrame(df_results)


# --- Carga y Preprocesamiento de Datos ---
print("Cargando tablas biom...")
train = biom.load_table(train_file)
test = biom.load_table(test_file)

train_col = train.ids(axis='observation')
test_col = test.ids(axis='observation')
train_indx = train.ids(axis='sample')
test_indx = test.ids(axis='sample')

print("Normalizando datos...")
train = train.rankdata(axis='sample', inplace=False)
train = train.matrix_data.multiply(1 / train.max(axis='sample'))
test = test.rankdata(axis='sample', inplace=False)
test = test.matrix_data.multiply(1 / test.max(axis='sample'))

otu_train = pd.DataFrame(columns=train_col, index=train_indx, data=train.toarray().T)
otu_test = pd.DataFrame(columns=test_col, index=test_indx, data=test.toarray().T)

print("Cargando metadatos...")
map_keep = pd.read_csv(metadata, sep="\t", index_col=0, low_memory=False)

# Asegurarse de que la columna objetivo es numérica
if not pd.api.types.is_numeric_dtype(map_keep[group]):
    print(f"¡ERROR! La columna objetivo '{group}' no es numérica. La regresión requiere un objetivo continuo.")
    sys.exit(1)
# Eliminar filas con valores NaN en la columna objetivo
map_keep.dropna(subset=[group], inplace=True)


otu_train = otu_train.loc[np.intersect1d(otu_train.index.values, map_keep.index.values), :]
print(f"Forma de OTU de entrenamiento después de cruzar con metadatos: {otu_train.shape}")
print(f"Forma de OTU de prueba inicial: {otu_test.shape}")

otu_test = otu_test.reindex(columns=otu_train.columns, fill_value=0)
otu_test = otu_test.loc[np.intersect1d(otu_test.index.values, map_keep.index.values), otu_train.columns]
print(f"Forma de OTU de prueba después de alinear columnas y cruzar con metadatos: {otu_test.shape}")


if embedding_file != 'None':
    print("Aplicando embeddings...")
    embedding_vector = pd.read_csv(embedding_file, sep=" ", index_col=0, header=None, low_memory=False)
    embedding_vector.index = [str(i) for i in embedding_vector.index]
    fid = np.intersect1d(otu_train.columns.values, embedding_vector.index.values)
    otu_train = otu_train.loc[:, fid].dot(embedding_vector.loc[fid, :])
    otu_test = otu_test.loc[:, fid].dot(embedding_vector.loc[fid, :])
    print("Forma de datos de entrenamiento después de embeddings:", otu_train.shape)

X_train = otu_train
X_test = otu_test
y_train = map_keep.loc[X_train.index, group]
y_test = map_keep.loc[X_test.index, group]
test_sample_ids = X_test.index.values

# --- Ajuste de Hiperparámetros (Opcional) ---
# Descomenta las siguientes líneas para ejecutar la búsqueda de hiperparámetros.
# Es computacionalmente intensivo.
# print("Iniciando ajuste de hiperparámetros...")
# df_params = trainHyperParameters_regression(X_train, y_train)
# df_params = df_params.sort_values(by=['R2_mean'], ascending=False)
# print("Resultados del ajuste de hiperparámetros:")
# print(df_params)
# best_depth = int(df_params['depth'].values[0])
# best_trees = int(df_params['n_trees'].values[0])
# print(f"Mejores parámetros encontrados: Profundidad={best_depth}, Estimadores={best_trees}")

# Usar hiperparámetros fijos o los mejores encontrados
# Si no realizas el ajuste, establece los valores aquí.
best_depth = 10
best_trees = 500


# --- Entrenamiento y Evaluación del Modelo Final ---
print("\nEntrenando el modelo final con los mejores parámetros...")
model = RandomForestRegressor(n_estimators=best_trees,
                              max_depth=best_depth,
                              random_state=0,
                              n_jobs=24)
model.fit(X_train, y_train)

print("Guardando importancias de las características...")
importances = model.feature_importances_
feature_weights = pd.DataFrame({
    "feature_id": otu_train.columns,
    "feature_importance": importances
})
feature_weights.sort_values(by='feature_importance', ascending=False, inplace=True)
feature_weights.to_csv(feature_weights_path, index=False)
print(f"Las importancias de las características se han guardado en: {feature_weights_path}")

print("\nEvaluando el modelo en el conjunto de prueba...")
r2, mae, mse, rmse, _ = computeRegressionStats(
    model,
    data=X_test,
    y_true=y_test,
    sample_ids=test_sample_ids,
    predictions_output_file=PREDICTED_VALUES_FILE,
    plot=True,
    plot_file=plot_file,
    graph_title="Rendimiento del Modelo de Regresión en el Conjunto de Prueba"
)

print("\n--- Resultados Finales de la Evaluación ---")
print(f"R-cuadrado (R^2): {r2:.4f}")
print(f"Error Absoluto Medio (MAE): {mae:.4f}")
print(f"Error Cuadrático Medio (MSE): {mse:.4f}")
print(f"Error Cuadrático Medio Raíz (RMSE): {rmse:.4f}")
print("----------------------------------------")