import warnings
warnings.filterwarnings('ignore')
import pandas as pd
import requests
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix,ConfusionMatrixDisplay,f1_score,jaccard_score

df = pd.read_csv('Human_Cancer_Genes_Data.csv')


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

# get protein sequence from protein name 


def get_protein_sequence(protein_name):
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    query = f"query={protein_name}&fields=accession,sequence"
    headers = {"Accept": "application/json"}
    
    response = requests.get(f"{base_url}?{query}", headers=headers)
    
    if response.status_code == 200:
        data = response.json()
        results = data.get("results", [])
        if results:
            # Access the first result (or iterate for multiple hits)
            protein_data = results[0]
            accession = protein_data.get("primaryAccession", "Unknown")
            sequence = protein_data.get("sequence", {}).get("value", "No sequence found")
            return accession, sequence
        else:
            return None, "No matching protein found."
    else:
        return None, f"Failed to fetch data. Status code: {response.status_code}"

@app.route('/', methods=['GET','POST'])
def hello():
    if request.method == 'POST':
        protein_Name = request.form['proteinInput']
        accession, sequence = get_protein_sequence(protein_Name)
        if not accession:
            return "Errror!! Not a protein sequence."+"\nPlease enter a protein sequence"


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
            res = 'This protein is related to cancer'
        else:
            res = 'This protein is not related to cancer'

        result = {'model':"Logistic Regression", 'pred':res, 'fsc':fscore, 'jsc':jscore}
        return render_template('result.html', result=result)
        
    return render_template('form.html')


if __name__=='__main__':
    app.run(debug=True)

