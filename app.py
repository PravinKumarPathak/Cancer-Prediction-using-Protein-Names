import warnings
warnings.filterwarnings('ignore')
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix,ConfusionMatrixDisplay,f1_score,jaccard_score
<<<<<<< HEAD
=======
#from matplotlib import pyplot as plt
>>>>>>> 8ebf39e7ac3bb8437b32b8ff7794d173e7365177

df = pd.read_csv('Human_Cancer_Genes_Data.csv')


<<<<<<< HEAD
=======
#print(df.head())
>>>>>>> 8ebf39e7ac3bb8437b32b8ff7794d173e7365177

print(df.columns)

cancer_df = df.drop(['Entry','Entry Name', 'Protein names','Gene Names','Sequence', 'Mass'], axis=1)

print(cancer_df.head())


print(cancer_df['target'].value_counts())



X = cancer_df.drop(columns=['target'])
Y = cancer_df['target']

X_train, X_test, Y_train, Y_test = train_test_split(X,Y,test_size=0.2, random_state=4)

lr = LogisticRegression(C=0.01, solver='liblinear')
lr.fit(X_train, Y_train)
yhat = lr.predict(X_test)

fscore = round(f1_score(Y_test, yhat, average='weighted'), 4)
jscore = round(jaccard_score(Y_test, yhat, pos_label=0), 4)
print("Avg F1-score:%0.4f" %fscore)
print("Jaccard score:%0.4f" %jscore)



from flask import Flask, render_template, request

app = Flask(__name__)

from Bio.SeqUtils.ProtParam import ProteinAnalysis

def compute_protein_features(sequence):
    try:
        analysis = ProteinAnalysis(sequence)
        return {
            "Hydrophobicity": analysis.gravy(),
            "Isoelectric Point": analysis.isoelectric_point(),
            "Aromaticity": analysis.aromaticity(),
        }
    except Exception as e:
        return {"Hydrophobicity": 0, "Isoelectric Point": 0, "Aromaticity": 0}

def amino_Sequence(sequence):
    amino_seq = {}
    for aa in "ACDEFGHIKLMNPQRSTVWY" :
        if len(sequence) > 0 :
            amino_seq[f"AA_{aa}"] = sequence.count(aa) / len(sequence)
    return amino_seq


@app.route('/', methods=['GET','POST'])
def hello():
<<<<<<< HEAD
=======
    #return f'Avg F1-score:{f1_score(Y_test, yhat, average='weighted')}'
>>>>>>> 8ebf39e7ac3bb8437b32b8ff7794d173e7365177
    if request.method == 'POST':
        sequence = request.form['proteinInput']

        features = {
            'Sequence' : sequence
        }

        Length = len(sequence)
        features['Length'] = Length

        protein_features = compute_protein_features(sequence)
        features.update(protein_features)
        features.update(amino_Sequence(sequence))

        print(features)
        del features['Sequence']
        print(features)
        print(type(features))
        print(len(features))

        test_data = pd.DataFrame([features])
        print(test_data)

        pred = lr.predict(test_data)
        print(pred)
        if(pred[0]==1):
            res = 'This protein sequence is related to cancer'
        else:
            res = 'This protein sequence is not related to cancer'

        result = {'model':"Logistic Regression", 'pred':res, 'fsc':fscore, 'jsc':jscore}
        return render_template('result.html', result=result)
        
    return render_template('form.html')


if __name__=='__main__':
    app.run(debug=True)
<<<<<<< HEAD


=======
>>>>>>> 8ebf39e7ac3bb8437b32b8ff7794d173e7365177
