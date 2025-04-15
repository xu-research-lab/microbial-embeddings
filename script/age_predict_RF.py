import biom
import pandas as pd
import numpy as np
from sklearn.model_selection import KFold

from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score, mean_absolute_error

import seaborn as sns
import matplotlib.pyplot as plt

def age_predict_fivefold(metadata_file, biom_path, plot_parameter, save_path):
    
    metadata = pd.read_csv(metadata_file, sep='\t', index_col="#SampleID")
    repotred_age = np.array([])
    predicted_age = np.array([])
    r2 = np.array([])
    MAE = np.array([])
    
    # five fold train and test
    for i in range(1,6):
        train_table = biom.load_table(biom_path + f"train_table_{i}.biom")
        train_labels = np.array(metadata.loc[train_table.ids(axis="sample"),].age)
        train_table = train_table.matrix_data.toarray().T
        
        test_table = biom.load_table(biom_path + f"test_table_{i}.biom")
        test_labels = np.array(metadata.loc[test_table.ids(axis="sample"),].age)
        test_table = test_table.matrix_data.toarray().T
        
        #
        rf_regressor = RandomForestRegressor(n_estimators=500, n_jobs=20,
                                             random_state=42, criterion="squared_error") # squared_error, absolute_error
        rf_regressor.fit(train_table, train_labels)
        y_pred = rf_regressor.predict(test_table)
        repotred_age = np.append(repotred_age, test_labels)
        predicted_age = np.append(predicted_age, y_pred)
        r2 = np.append(r2, r2_score(test_labels, y_pred))
        MAE = np.append(MAE, mean_absolute_error(test_labels, y_pred))
    
    MAE_mean = np.mean(MAE)
    MAE_sd = np.std(MAE)
    r2_mean = np.mean(r2*100)
    r2_sd = np.std(r2*100)
    
    data = pd.DataFrame({'Reported age': repotred_age, 'Predicted age': predicted_age})
    
    plt.figure(figsize=(6, 10))
    sns.set_style("whitegrid")
    sns.lmplot(x='Reported age', y='Predicted age', data=data, fit_reg=True,
               scatter_kws={'alpha':0.5},
               scatter=True, ci=None, order=2, line_kws={'color': 'red'})
    plt.text(plot_parameter['x'][0], plot_parameter['y'][0], plot_parameter['title'], fontsize=10, color='black')
    plt.text(plot_parameter['x'][1], plot_parameter['y'][1], f"MAE: {MAE_mean:.2f}±{MAE_sd:.2f} yr", fontsize=10, color='black')
    plt.text(plot_parameter['x'][2], plot_parameter['y'][2], f"R squared: {r2_mean:.2f}±{r2_sd:.2f}%", fontsize=10, color='black')
    plt.xlabel('Reported age',fontsize=12)
    plt.ylabel('Predicted age',fontsize=12)
    plt.title('RF', fontsize=12, loc='left')
    plt.savefig(plot_parameter['path'], dpi=300, bbox_inches='tight')


## all data
metadata_file = "/home/dongbiao/word_embedding_microbiome/programe_test/age_predict/age-prediction-master/modeldata/gut/metadata.txt"
biom_path = "/home/dongbiao/word_embedding_microbiome/programe_test/age_predict/age-prediction-master/modeldata/gut/five_fold/"
plot_parameter = {}
plot_parameter['x'] = [60, 57, 55]
plot_parameter['y'] = [70, 68, 66]
plot_parameter['title'] = "gut microbiota"
plot_parameter['path'] = "/home/dongbiao/word_embedding_microbiome/programe_test/age_predict/age-prediction-master/modeldata/gut/all_age.png"
save_path = "/home/dongbiao/word_embedding_microbiome/programe_test/age_predict/age-prediction-master/modeldata/gut/five_fold/res/"
age_predict_fivefold(metadata_file, biom_path, plot_parameter, save_path)

## sub data
### men
metadata_file = "/home/dongbiao/word_embedding_microbiome/programe_test/age_predict/age-prediction-master/modeldata/gut/subdata/man/metadata.txt"
biom_path = "/home/dongbiao/word_embedding_microbiome/programe_test/age_predict/age-prediction-master/modeldata/gut/subdata/man/"
plot_parameter = {}
plot_parameter['x'] = [57, 57, 55]
plot_parameter['y'] = [68, 66, 64]
plot_parameter['title'] = "men gut microbiota"
plot_parameter['path'] = "/home/dongbiao/word_embedding_microbiome/programe_test/age_predict/age-prediction-master/modeldata/gut/subdata/man/res/men_age.png"
save_path = "/home/dongbiao/word_embedding_microbiome/programe_test/age_predict/age-prediction-master/modeldata/gut/subdata/man/res/"
age_predict_fivefold(metadata_file, biom_path, plot_parameter, save_path)

### women
metadata_file = "/home/dongbiao/word_embedding_microbiome/programe_test/age_predict/age-prediction-master/modeldata/gut/subdata/woman/metadata.txt"
biom_path = "/home/dongbiao/word_embedding_microbiome/programe_test/age_predict/age-prediction-master/modeldata/gut/subdata/woman/"
plot_parameter = {}
plot_parameter['x'] = [55, 57, 55]
plot_parameter['y'] = [67, 65, 63]
plot_parameter['title'] = "women gut microbiota"
plot_parameter['path'] = "/home/dongbiao/word_embedding_microbiome/programe_test/age_predict/age-prediction-master/modeldata/gut/subdata/woman/res/women_age.png"
save_path = "/home/dongbiao/word_embedding_microbiome/programe_test/age_predict/age-prediction-master/modeldata/gut/subdata/woman/res/"
age_predict_fivefold(metadata_file, biom_path, plot_parameter, save_path)
