from flask import Flask, render_template, url_for, request, redirect
from mordred._base.descriptor import Descriptor
from similarity_calc import DataLoader, Similarity
import predictor
from rdkit import Chem
from rdkit.Chem import Draw

app = Flask(__name__)


@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        return render_template('similarity.html')
    return render_template('index.html')

@app.route('/similarity', methods=['POST', 'GET'])
def similarity():
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
                                round(similar[0][1],2), 
                                int(similar[2].iloc[1]*100),
                                int(similar[2].iloc[2]*100),
                                int(similar[2].iloc[3]*100),
                                int(similar[2].iloc[4]*100),
                                int(similar[2].iloc[5]*100),
                                int(similar[2].iloc[6]*100)]
        try:
            return render_template("similarity.html", tasks = retrieved_features)
        except:
            return 'There was an error adding your task'

    if request.method == "GET":
        return render_template("similarity.html")


@app.route('/pls-prediction', methods=['POST', 'GET'])
def pls_prediction():
    if request.method == 'POST':
        user_smiles = request.form['content']
        #comparing_dataset = DataLoader.csv_loader('data/osn_data_similarity.csv')
        try:
        #    similar = Similarity.rejection_similarity(comparing_dataset, reference=user_smiles)
            user_descr = predictor.descripter(user_smiles)
            prediction = predictor.predictor(user_descr)
        except:
            return 'Wrong SMILES, please try again. Currently, molecules with more than 5 heavy atoms work'

        Draw.MolToFile(Chem.MolFromSmiles(user_smiles),'static/dataset_image.png') 
        retrieved_features = [user_smiles, 
                                round(prediction[0][0], 3)*100]
        try:
            return render_template("pls-prediction.html", tasks = retrieved_features)
        except:
            return 'There was an error adding your task'
    if request.method == "GET":
        return render_template("pls-prediction.html")


@app.route('/publications', methods=['GET', 'POST'])
def publications():
    if request.method == 'POST':
        return render_template('publications.html')
    return render_template('publications.html')


@app.route('/contact', methods=['GET', 'POST'])
def contact():
    if request.method == 'POST':
        return render_template('contact.html')
    return render_template('contact.html')


if __name__ == "__main__":
    app.run(debug=False)