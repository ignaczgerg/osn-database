from flask import Flask, render_template, url_for, request, redirect
from similarity_calc import DataLoader, Similarity

from rdkit import Chem
from rdkit.Chem import Draw

app = Flask(__name__)

@app.route('/', methods=['POST', 'GET'])
def index():
    if request.method == 'POST':
        user_smiles = request.form['content']
        comparing_dataset = DataLoader.csv_loader('data/osn_data_similarity.csv')
        try:
            similar = Similarity.rejection_similarity(comparing_dataset, reference=user_smiles)
        except:
            return 'Wrong SMILES, please try again'

        Draw.MolToFile(similar[1],'static/dataset_image.png') 
        retrieved_features = [user_smiles, 
                                str(similar[2].iloc[0]),
                                round(similar[0][1],3), 
                                similar[2].iloc[1],
                                similar[2].iloc[2],
                                similar[2].iloc[3],
                                similar[2].iloc[4],
                                similar[2].iloc[5],
                                similar[2].iloc[6]]
        try:
            return render_template("index.html", tasks = retrieved_features)
        except:
            return 'There was an issue adding your task'

    if request.method == "GET":
        return render_template('index.html')


if __name__ == "__main__":
    app.run()